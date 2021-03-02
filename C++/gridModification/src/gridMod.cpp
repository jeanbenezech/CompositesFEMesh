#include <gridMod.h>
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

	// read & initialisation du mesh
	Mesh m;
	m.mesh_type="gmsh";
	size_t lastindex = param.meshname.find_last_of(".");
	std::string mesh_name = param.meshname.substr(0, lastindex);
	m.read_msh(mesh_name+".msh", true, param.cz_id);

	m.initialise(param.nb_plies);
	m.read_points("input.txt");
	m.exportDir = true;

	if(param.Shape==0)
		localCoorSyst(m);
	else
		globalCoorSyst(m);

	GeometricTransformation(m, param);
	StackingSequence(m);

	if (param.Abaqus_output){ // ABAQUS
		m.write_ori_inp("Abaqus/"+mesh_name);
		m.write_inp("Abaqus/"+mesh_name);
		m.write_abaqus_cae_input("Abaqus/"+mesh_name);
	}

	if (param.Dune_output){ // DUNE
		m.write_ori_txt("DUNE/"+mesh_name);
		m.write_msh("DUNE/"+mesh_name);
	}

	m.write_vtk(mesh_name);

	return 0;
}
