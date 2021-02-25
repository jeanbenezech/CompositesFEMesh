#include <mesh.h>
#include <parameters.h>
#include <utils.h>
#include <cstring>
#include <string>
#include <sstream>

using namespace Eigen;

int main(int argc, const char *argv[]) {
	Parameters param;
	param.read("parameters.txt");

	// lecture & initialisation du maillage
	Mesh m;
	m.mesh_type="gmsh";
	size_t lastindex = param.meshname.find_last_of(".");
	std::string mesh_name = param.meshname.substr(0, lastindex);
	m.read_msh(mesh_name+".msh", true, param.cz_id);

	m.read_elem_fields(param.path_to_abaqus_result+param.AbaqusOdbName+".txt");
	m.read_node_fields(param.path_to_abaqus_result+param.AbaqusOdbName+"_U.txt");

	m.exportAbaqusFields = true;

	m.write_vtk("Abaqus_result");

	return 0;
}
