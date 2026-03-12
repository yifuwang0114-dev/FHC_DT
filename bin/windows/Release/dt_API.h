#ifndef dt_API_h
#define dt_API_h

#ifdef __linux__
#define DECL_VOLTET
#elif _WIN32
#define DECL_VOLTET __declspec(dllexport)
#else
#error "Unsupported platform"
#endif

#include <functional>
#include <vector>
#include <array>
#include <string>
#include <set>
#include <map>
#include <unordered_map>

namespace dt {
	struct Mesh {//for io easily
		std::vector<std::array<double, 3>> V;
		std::vector<std::array<int, 3>> S;
		std::vector<std::array<int, 4>> F;
		std::vector<std::array<int, 5>> T;
	};

    struct Args {
        int infolevel = 1;          // Verbosity level of output messages (1每3)
        int refine = 0;             // Enable mesh refinement and optimization (0 = off, 1 = on)
        int optloop = 5;            // Number of optimization iterations
        double size = 1;            // Relative mesh size (smaller value produces denser meshes)
        double growsize = 1.07;     // Mesh size growth rate between neighboring elements
        double optangle = 10;       // Target minimum dihedral angle for optimization (0每70 degrees)
        double optratio = 0.1;      // Target volume每edge ratio for mesh quality (0每1)
        std::string filename;      // Input geometry file name
        std::vector<double> hole;  // Hole coordinates: x1,y1,z1, x2,y2,z2, ...
        std::function<double(const double&, const double&, const double&)> sizingFunc = nullptr;
        // User-defined sizing field function: size = f(x, y, z)
    };
}

#ifdef DT_LIBRARY
/*
 *  @brief	Main Tetrahedralize
 *
 *  @param[in]		mesh		input mesh
 *  @param[in]		args		input arguments
 *  @param[out]		mesh		output mesh
 *  @return
 *	1				Completely succeed
 *	otherwise		Fail. Error info. will be defined later
 */
DECL_VOLTET int API_Tetrahedralize(
	dt::Mesh& mesh, dt::Args args
);

/*
 *  @brief	Create Mesh from vtk file
 *
 *  @param[in]		filename		input filename
 *  @param[out]		mesh			output mesh
 *  @return
 *	1				Completely succeed
 *	otherwise		Fail. Error info. will be defined later
 */
DECL_VOLTET int API_ReadMesh_File(
	std::string in_filename, dt::Mesh& mesh
);

/*
 *  @brief	Create Mesh from vtk file
 *
 *  @param[in]		filename		output filename
 *  @param[in]		mesh			output mesh
 */
DECL_VOLTET int API_WriteMesh(
	std::string out_filename, dt::Mesh& mesh, bool addSurTri = false
);
#endif
#endif
