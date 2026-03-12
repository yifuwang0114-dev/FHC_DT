import os
import subprocess
import pandas as pd
from pathlib import Path
import time
import re
import logging
import psutil
from concurrent.futures import ThreadPoolExecutor, as_completed

# Debug mode parameter: If not 0, only the debug_mode prefix will be tested.
debug_mode = 0  # Set to 0 to run all tasks, set to n to run only the n prefix tasks

# obj_dir = Path('./nointersection_meshes_obj')
obj_dir = Path('../examples/Thingi10K_Valid')
# mesh_dir = Path('../examples/Thingi10K_Valid_mesh')  # Still used for tetgen_Y

# Get the file prefix of obj_dir
common_prefixes = sorted([f.stem for f in obj_dir.glob('*.obj')])

# If debug mode is enabled, only the first debug_mode tasks will be retrieved.
if debug_mode > 0:
    common_prefixes = common_prefixes[:debug_mode]

commands = [
    # CDT 
    # {
    #     'type': 'cdt',
    #     'exe': './cdt.exe',
    #     'command': lambda prefix: ['./cdt.exe', str(obj_dir / f"{prefix}.obj"), '-vwrf'],
    #     'success_indicator': 'ST:',
    #     'time_pattern': r"Time:\s+([\d\.]+)", 
    #     'boundary_pattern': r"(\d+) nononono",
    #     'volume_pattern': r"ST:\s+(\d+)", 
    #     'memory_pattern': r"Mem:\s+([\d\.]+)"
    # },
    # {
    #     'type': 'cdt_float',
    #     'exe': './cdt.exe',
    #     'command': lambda prefix: ['./cdt.exe', str(obj_dir / f"{prefix}.obj"), '-vwfr'],
    #     'success_indicator': 'ST:',
    #     'time_pattern': r"Time:\s+([\d\.]+)",  
    #     'boundary_pattern': r"(\d+) nononono",
    #     'volume_pattern': r"ST:\s+(\d+)", 
    #     'memory_pattern': r"Mem:\s+([\d\.]+)"
    # },
    # Tetgen 
    # {
    #     'type': 'tetgen',
    #     'exe': './tetgen.exe',
    #     'command': lambda prefix: ['./tetgen.exe', str(mesh_dir / f"{prefix}.mesh"), '-pYNEFk'],
    #     'success_indicator': 'Statistics:',
    #     'time_pattern': r"Total running seconds:\s+([\d\.]+)",
    #     'volume_pattern': r"Steiner points inside domain:\s+(\d+)",
    #     'boundary_pattern': r"Steiner points on input segments:\s+(\d+)|Steiner points on input facets:\s+(\d+)",
    #     'memory_pattern': r"Memory cost\s+:\s+([\d\.]+)\s+MB"
    # },
    # dt 
    {
        'type': 'dt',
        'exe': 'dt.exe',
        'command': lambda prefix: [
            './dt.exe', '--input', str(obj_dir / f"{prefix}.obj"), '--refine', '0'
        ],
        'success_indicator': 'Tetrahedralize success!',
        'volume_pattern': r"(\d+) steiner in volume",
        'boundary_pattern': r"(\d+) steiner in Boundary",
        'memory_pattern': r"Memory cost\s+:\s+([\d\.]+)\s+MB",
        'time_pattern': r"Total cost\s+:\s+([\d\.]+)\s+s" 
    },
]

logging.basicConfig(
    level=logging.INFO,
    format='%(message)s',
    handlers=[logging.FileHandler("ans.txt", mode='w', encoding='utf-8'), logging.StreamHandler()]
)

def extract_steiner_points(output, volume_pattern, boundary_pattern, memory_pattern, time_pattern=None):

    volume_match = re.search(volume_pattern, output)
    boundary_matches = re.findall(boundary_pattern, output)
    memory_matches = re.findall(memory_pattern, output)

    # Obtain Steiner points within the body
    volume_points = int(volume_match.group(1)) if volume_match else 0
    # Get the boundary Steiner point, if it is 0
    boundary_points = sum(int(num) for match in boundary_matches for num in match if num.isdigit()) if boundary_matches else 0
    memory_cost = float(memory_matches[0]) if memory_matches else 0  # Use memory_matches[0] to get the first matching value

    time_value = 0
    if time_pattern:
        time_match = re.search(time_pattern, output)
        time_value = float(time_match.group(1)) if time_match else 0

    return time_value, volume_points, boundary_points, memory_cost

def check_available_memory(min_memory_gb=1):
    """
    Check system available memory; if it is below the specified value, wait.
    :param min_memory_gb: Minimum available memory, in GB
    """
    while True:
        # Get available memory (in bytes)
        available_memory = psutil.virtual_memory().available
        available_memory_gb = available_memory / (1024 ** 3)  # Convert to GB
        logging.info(f"Current available memory: {available_memory_gb:.2f} GB")

        if available_memory_gb >= min_memory_gb:
            break
        else:
            logging.info(f"Insufficient memory, waiting for {min_memory_gb} GB available memory...")
            time.sleep(5)  # Wait for 5 seconds before checking again

def execute_command(prefix, cmd_info, timeout=7200, memory_limit_gb=20):
    cmd = cmd_info['command'](prefix)
    start_time = time.perf_counter()
    output = ""
    success = False

    try:
       # Start child process
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding='utf-8',
            errors='ignore'
        )

        output, _ = proc.communicate(timeout=timeout)
        success = cmd_info['success_indicator'] in output

    except subprocess.TimeoutExpired:
        logging.error(f"{prefix} | {cmd_info['type']} Timeout, terminate process")
        proc.kill()
        output, _ = proc.communicate()
        return False, 9999999, 0, 0, 0

    except Exception as e:
        logging.error(f"{prefix} | {cmd_info['type']} Execution failed: {e}")
        return False, -1, 0, 0, 0

    # Extract metrics
    time_value, volume_points, boundary_points, memory_cost = extract_steiner_points(
        output, cmd_info['volume_pattern'], cmd_info['boundary_pattern'], 
        cmd_info['memory_pattern'], cmd_info.get('time_pattern')
    )

    return success, time_value, volume_points, boundary_points, memory_cost

# Parallel execution of tasks (modifying the result record section)
max_workers = min(1, os.cpu_count())
logging.info(f"Start execution {len(common_prefixes) * len(commands)} ,use {max_workers} threads...\n")

results_dict = {prefix: {} for prefix in common_prefixes}

with ThreadPoolExecutor(max_workers=max_workers) as executor:
    futures = {
        executor.submit(execute_command, prefix, cmd): (prefix, cmd['type'])
        for prefix in common_prefixes for cmd in commands
    }

    for i, future in enumerate(as_completed(futures), start=1):
        prefix, cmd_type = futures[future]

        check_available_memory(min_memory_gb=1)


        success, time_value, volume_points, boundary_points, max_memory = future.result()

        print(future.result())
        results_dict[prefix][f"{cmd_type}_ans"] = "success" if success else "failure"
        results_dict[prefix][f"{cmd_type}_time"] = time_value
        results_dict[prefix][f"{cmd_type}_volume"] = volume_points
        results_dict[prefix][f"{cmd_type}_boundary"] = boundary_points
        results_dict[prefix][f"{cmd_type}_memory"] = max_memory

        logging.info(
            f"{prefix} | {cmd_type} | {'success' if success else 'failure'} | "
            f"{time_value}s | Steiner_volume: {volume_points} | "
            f"Steiner_boundary: {boundary_points} | memory max: {max_memory:.2f} MB"
        )
        logging.info(f"schedule: {i}/{len(futures)}\n")

df = pd.DataFrame.from_dict(results_dict, orient='index').reset_index()
df.rename(columns={'index': 'prefix'}, inplace=True)
excel_file = './ans.xlsx'
df.to_excel(excel_file, index=False)
logging.info(f"File saved to {excel_file}\n")