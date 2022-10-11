#ifndef VERTICE_H
#define VERTICE_H

#include "utils.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>

using namespace Eigen;

class Vertice {
public:


	Vector3d coord;
	Vector3d normal;
	// std::vector<Vector2i> element_ids;

	int top_surf_indice;
	int weight;

	double thickness;

	Vector3d X;

	// FUNCTION
	// void initialise(std::string type, int nb_elem, int nb_marker, bool isShell_);
	// void fill(std::vector<std::vector<int>> l, int& nb_elset);

	//~Vertice();
private:
};

// ~~~~~~~~~~~~~~ INIT ~~~~~~~~~~~~~~

// void Vertice::initialise() {

// }

// void Vertice::fill() {

// }


#endif /* end of include guard: VERTICE_H */
