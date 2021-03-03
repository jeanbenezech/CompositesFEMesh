#ifndef ELEMENT_H
#define ELEMENT_H

#include "utils.h"
#include <cmath>
#include <cstdlib>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>

using namespace Eigen;

class Element {
public:

	std::string type;
	int nb;
	int msh_type;
	int vtk_type;
	int size;
	std::string inria_mesh_name;
	std::string abaqus_type;

	// PRINCIPAL VARIABLES
	VectorXi global_indices;
	Matrix<int, Dynamic, Dynamic> Nodes;
	Matrix<int, Dynamic, Dynamic> Markers;

	Matrix<float, Dynamic, Dynamic> U;
	Matrix<float, Dynamic, Dynamic> V;
	Matrix<float, Dynamic, Dynamic> W;

	// FUNCTION
	void initialise(std::string type, int nb_elem, int nb_marker);
	void fill(std::vector<std::vector<int>> l, int& nb_elset);

	//~Element();
private:
};

// ~~~~~~~~~~~~~~ INIT ~~~~~~~~~~~~~~

void Element::initialise(std::string key, int nb_elem, int nb_marker) {
	type = key;
	nb = nb_elem;
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
	Markers.resize(nb_marker, nb);
	global_indices.resize(nb);

	U.resize(3, nb);
	V.resize(3, nb);
	W.resize(3, nb);
}

void Element::fill(std::vector<std::vector<int>> l, int& nb_set) {

	for(int i = 0; i < nb; i++) {
		int beg_node = 2+l[i][2]+1;
		for (int j = 0; j < size; j++)
			Nodes(j,i) = l[i][beg_node+j]-1;

		global_indices[i] = l[i][0];
		Markers(0,i) = l[i][3]; // Maker[0] = physical group
		if (Markers(0,i)>nb_set) {nb_set=Markers(0,i);}
	}

}


#endif /* end of include guard: ELEMENT_H */
