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

	m.initialise(param);
	m.read_points("input.txt");
	m.exportDir = true;

	if(param.Shape==0)
		localCoorSyst(m, param);
	else
		globalCoorSyst(m);

	GeometricTransformation(m, param);
	StackingSequence(m, param);

	if (param.Abaqus_output){ // ABAQUS
		m.write_ori_inp("Abaqus/"+mesh_name);
		m.write_inp("Abaqus/"+mesh_name);
		m.write_abaqus_cae_input("Abaqus/"+mesh_name, param);
	}

	if (param.Dune_output){ // DUNE
		if (param.add_wrinkles>=1){
			m.write_ori_txt("DUNE/"+param.WID+"_"+mesh_name);
			m.write_msh("DUNE/"+param.WID+"_"+mesh_name);
		} else {
			m.write_ori_txt("DUNE/"+mesh_name);
			m.write_msh("DUNE/"+mesh_name);
		}
	}

	m.write_vtk(mesh_name);

	return 0;
}
