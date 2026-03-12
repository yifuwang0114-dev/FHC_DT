#include "dt_API.h"
using namespace std;

int main(int argc, char *argv[])
{
	std::string in_filename;
	std::string out_filename;

	dt::Args args;
	dt::Mesh mesh;
	dt::DT d;
	int ret;

	if (dt::readMesh(in_filename, mesh))
	{
		if (d.tetrahedralize(mesh, args))
			dt::writeVTK(out_filename, mesh, false);
		else
			printf("tetrahedralize failed!\n");
	}
	else
		printf("readVtk failed!\n");

	return 0;
}
