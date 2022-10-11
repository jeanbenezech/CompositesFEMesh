#include <CTData.h>
#include <gridMod.h>
#include <gridTransformation.h>
#include <mesh.h>
#include <parameters.h>
#include <utils.h>
#include <cstring>
#include <string>
#include <sstream>

// #include <filesystem>

using namespace Eigen;

int main(int argc, const char *argv[]) {

	// std::filesystem::create_directories("DUNE");
	
	std::string dataname = argv[1];
	CT_DATA ct;
	ct.readMeasurements(dataname);

	// std::cout << ct.data[0].Coord(2) << std::endl;

	Parameters param;
	param.read("parameters.txt");
	// param.add_wrinkles = true;

	// read & initialisation du mesh
	Mesh m;
	m.mesh_type="gmsh";
	size_t lastindex = param.meshname.find_last_of(".");
	std::string mesh_name = param.meshname.substr(0, lastindex);
	m.isShell = param.isShell;
	m.read_msh(mesh_name+".msh", true, param.cz_id);

	m.initialise(param);

	m.read_points("input.txt");
	m.exportDir = true;


	globalCoorSyst(m, param);


	CT_GeometricTransformation(m, ct, param);
	// StackingSequence(m, param);


	m.write_vtk(mesh_name);



	return 0;
}
