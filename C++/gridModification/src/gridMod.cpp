#include <gridMod.h>
#include <mesh.h>
#include <parameters.h>
#include <utils.h>
#include <cstring>
#include <string>
#include <sstream>

// #include <experimental/filesystem> // TODO:: No working with ARM os
// namespace fs = std::experimental::filesystem; // TODO:: No working with ARM os

using namespace Eigen;

int main(int argc, const char *argv[]) {
	// for (int i=0; i<101; i++){

	Parameters param;
	param.read("parameters");
	// param.read("parameters_wrinkles.txt");
	// param.read("parameters_"+std::to_string(i)+".txt");

	// read & initialisation du mesh
	Mesh m;
	m.mesh_type="gmsh";
	size_t lastindex = param.meshname.find_last_of(".");
	std::string mesh_name = param.meshname.substr(0, lastindex);
	m.isShell = param.isShell;
	m.initialise(param);
	m.read_msh(mesh_name+".msh", true, param);

	m.read_points("input.txt");
	m.exportDir = true;

	if(param.Shape==0 || param.Shape==2){ // Cspar
		localCoorSyst(m, param);}
	else { // laminate
		globalCoorSyst(m, param);
	}

	if (param.Dune_output){ // DUNE
		attribute_weight(m, param);
	}

	// if (param.rotateRVE && param.Abaqus_output) {
	// 	// m.extract_AbaqusSets(); // To be done before rotation of the RVE to keep the numbering of master nodes consistent
	// }

	// if (param.rotateRVE && param.Abaqus_output) {
	// 	// m.extract_AbaqusSets(); // To be done before rotation of the RVE to keep the numbering of master nodes consistent
	// 	m.elSets_delam(); // To be done before rotation of the RVE to keep the numbering of master nodes consistent
	// }
	// m.elSets_delam(); // To be done before rotation of the RVE to keep the numbering of master nodes consistent

	GeometricTransformation(m, param);
	StackingSequence(m, param);


	if(param.RigidBoundary){
		Rigid_Boundary(m, param);
	}

	if (param.Abaqus_output){ // ABAQUS
		// fs::create_directories("Abaqus"); // TODO:: No working with ARM os
		// m.write_msh("Abaqus/results/"+mesh_name); // for visualisation: mesh with in it wrinkle
		m.write_ori_inp("Abaqus/"+mesh_name, param);
		m.write_inp("Abaqus/"+mesh_name, param);
		// m.write_abaqus_cae_input("Abaqus/"+mesh_name, param);

		if (param.writeTransformedMSH)
			m.write_msh("tranformed_"+mesh_name);
	}

	if (param.Dune_output){ // DUNE
		// fs::create_directories("DUNE"); // TODO:: No working with ARM os
		if (param.add_wrinkles>=1){
			m.write_ori_txt("DUNE/"+param.WID+"_"+mesh_name, param);
			m.write_msh("DUNE/"+param.WID+"_"+mesh_name);
		} else {
			m.write_ori_txt("DUNE/"+mesh_name, param);
			m.write_msh("DUNE/"+mesh_name);
		}
	}

	if (param.do_test_delam){
		TESTinitialiseDamage(m, param);
	}

	m.write_vtk(mesh_name, param);
	// m.write_vtk(mesh_name+"_"+std::to_string(i));
	// }

	return 0;
}
