#ifndef ELEMENT_H
#define ELEMENT_H

#include "utils.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include "vertice.h"

using namespace Eigen;

class Random_field_parameters {
	public:
	std::vector<Vector3d> nodes;
	std::vector<double> values;
	double average;
};

class Element {
public:

	bool isShell = false;
	std::string type;
	int nb;
	int msh_type;
	int vtk_type;
	int size;
	std::string inria_mesh_name;
	std::string abaqus_type;

	// PRINCIPAL VARIABLES
	std::vector<int> global_indices;
	Matrix<int, Dynamic, Dynamic> Nodes;
	Matrix<int, Dynamic, Dynamic> Markers;
	Matrix<int, Dynamic, Dynamic> DD_weight;

	Matrix<double, Dynamic, Dynamic> Initial_barycenter;

	std::vector<Random_field_parameters> Random_Field;
	std::vector<Random_field_parameters> CT_Field;

	Matrix<double, Dynamic, Dynamic> U;
	Matrix<double, Dynamic, Dynamic> V;
	Matrix<double, Dynamic, Dynamic> W;

	// ABAQUS FIELDS
	Matrix<double, Dynamic, Dynamic> S;
	Matrix<double, Dynamic, Dynamic> E;
	Matrix<double, Dynamic, Dynamic> SDEG;
	Matrix<double, Dynamic, Dynamic> QUADSCRT;

	// FUNCTION
	void initialise(std::string type, int nb_elem, int nb_marker, bool isShell_);
	void fill(std::vector<std::vector<int>> l, int& nb_elset, std::vector<int>& material_types);
	Vector3d center(std::vector<Vertice> vertices, int id);

	//~Element();
private:
};

// ~~~~~~~~~~~~~~ INIT ~~~~~~~~~~~~~~

Vector3d Element::center(std::vector<Vertice> vertices, int id){
	Vector3d out = {0.,0.,0.};
	// std::cout << this->Nodes.rows() << std::endl;
	for(int node_iter =0; node_iter<this->Nodes.rows(); node_iter++){
		out += vertices[this->Nodes(node_iter,id)].coord;
	}
	return out /= this->Nodes.rows();
}

void Element::initialise(std::string key, int nb_elem, int nb_marker, bool isShell_=false) {
	type = key;
	nb = nb_elem;
	isShell = isShell_;

	if(type=="Triangle"){
		size = 3;
		msh_type = 2;
		vtk_type = 5;
		inria_mesh_name = "Triangles";
		abaqus_type = "None";
	} else if(type=="Quadrilateral"){
		size = 4;
		msh_type = 3;
		vtk_type = 9;
		inria_mesh_name = "Quadrilaterals";
		abaqus_type = "None";
	} else if(type=="Tetrahedron"){
		size = 4;
		msh_type = 4;
		vtk_type = 10;
		inria_mesh_name = "Tetrahedra";
		abaqus_type = "C3D4";
	} else if(type=="Hexahedron"){
		size = 8;
		msh_type = 5;
		vtk_type = 12;
		inria_mesh_name = "Hexahedra";
		abaqus_type = "C3D8";
		if(isShell){
			abaqus_type = "SC8R";
		}
	} else if(type=="Cohesive"){
		size = 8;
		msh_type = 5;
		vtk_type = 12;
		inria_mesh_name = "Hexahedra";
		abaqus_type = "COH3D8";
	} else if(type=="Prism"){
		size = 6;
		msh_type = 6;
		vtk_type = 13;
		inria_mesh_name = "None";
		abaqus_type = "C3D6";
	} else if(type=="Pyramid"){
		size = 5;
		msh_type = 7;
		vtk_type = 14;
		inria_mesh_name = "None";
		abaqus_type = "C3D5";
	}
	Nodes.resize(size, nb);
	Markers.resize(nb_marker + 1, nb); // nb_marker=1 for the physical group +1 for the material type
	DD_weight.resize(1, nb);
	global_indices.resize(nb);

	Initial_barycenter.resize(3, nb);

	Random_Field.resize(nb);
	CT_Field.resize(nb);

	U.resize(3, nb);
	V.resize(3, nb);
	W.resize(3, nb);

	S.resize(6, nb);
	E.resize(6, nb);
	SDEG.resize(1, nb);
	QUADSCRT.resize(1, nb);
}

void Element::fill(std::vector<std::vector<int>> l, int& nb_set, std::vector<int>& material_types) {

	for(int i = 0; i < nb; i++) {
		int beg_node = 2+l[i][2]+1;
		for (int j = 0; j < size; j++)
			Nodes(j,i) = l[i][beg_node+j]-1;

		global_indices[i] = l[i][0];
		Markers(0,i) = l[i][3]; // Marker[0] = indice ply in incremental order from bottom to top surface
		Markers(1,i) = material_types[Markers(0,i)-1]; // Maker[1] = type of material behaviour
		if (Markers(0,i)>nb_set) {nb_set=Markers(0,i);}
	}

}


#endif /* end of include guard: ELEMENT_H */
