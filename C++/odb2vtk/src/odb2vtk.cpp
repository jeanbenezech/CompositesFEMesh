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
	// m.read_msh(param.path_to_abaqus_result+mesh_name+".msh", true, param.cz_id);
	m.read_msh(mesh_name+".msh", true, param.cz_id);


	m.read_elem_fields_stresses_strains_only(param.path_to_abaqus_result+param.AbaqusOdbName+".txt");
	m.read_node_fields(param.path_to_abaqus_result+param.AbaqusOdbName+"_U.txt");

	m.exportAbaqusFields = true;
	m.exportAbaqusDisplacement = true;
	m.write_vtk(param.path_to_abaqus_result+param.AbaqusOdbName);

	// for(int i=1; i<42;i++){
	// 	m.read_elem_fields(param.path_to_abaqus_result+param.AbaqusOdbName+"_frame_"+std::to_string(i)+".txt");
	// 	m.read_node_fields(param.path_to_abaqus_result+param.AbaqusOdbName+"_frame_"+std::to_string(i)+"_U.txt");

	// 	m.exportAbaqusFields = true;
	// 	m.exportAbaqusDisplacement = true;

	// 	m.write_vtk(param.path_to_abaqus_result+param.AbaqusOdbName+"_frame_"+std::to_string(i));
	// }

	

	return 0;
}
