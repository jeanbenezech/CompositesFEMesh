#ifndef MAILLAGE_H
#define MAILLAGE_H

#include "utils.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <time.h>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Dense>

using namespace Eigen;

class Maillage {
public:
	// PRINCIPAL VARIABLES
	Matrix<float, Dynamic, Dynamic> vertices;
	Matrix<int, Dynamic, Dynamic> Triangles;
	Matrix<int, Dynamic, Dynamic> Quadrilaterals;
	Matrix<int, Dynamic, Dynamic> Tetrahedra;
	Matrix<int, Dynamic, Dynamic> Hexahedra;
	Matrix<int, Dynamic, Dynamic> Cohesive;
	Matrix<int, Dynamic, Dynamic> Wedges;
	Matrix<int, Dynamic, Dynamic> Pyramids;

	Matrix<double, Dynamic, Dynamic> S;
	Matrix<double, Dynamic, Dynamic> E;
	Matrix<double, Dynamic, Dynamic> U;

	Matrix<double, Dynamic, Dynamic> SDEG;
	Matrix<double, Dynamic, Dynamic> QUADSCRT;

	std::string mesh_type;

	int nb_nset=0;
	int nb_elset=0;

	Matrix<float, Dynamic, Dynamic> sortie;
	std::vector<float> stacking_sequence;
	Vector2f P1, P2, P3;

	// FUNCTION
	void lecture_mesh(const std::string& filename);
	void lecture_msh(const std::string& filename, bool only_3D, int cz_id);
	void lecture_points(const std::string& filename);
	void lecture_fields(const std::string& filename);
	void lecture_node_fields(const std::string& filename);

	void ecriture_mesh(const std::string& filename);
	void ecriture_msh(const std::string& filename);
	void ecriture_vtk(const std::string& filename);
	void ecriture_fields_vtk(const std::string& filename);
	void ecriture_inp(const std::string& filename);

	void ecriture_ori_txt(const std::string& filename);
	void ecriture_ori_inp(const std::string& filename);

	void initialise();
	void print_hex(int number);
	void copy(Maillage& m);
	void duplicates();

	int Nb_vertices(){return nb_vertices_;}
	int Nb_Tetrahedra(){return nb_Tetrahedra_;}
	int Nb_Triangles(){return nb_Triangles_;}
	int Nb_Quadrilaterals(){return nb_Quadrilaterals_;}
	int Nb_Hexahedra(){return nb_Hexahedra_;}
	int Nb_Wedges(){return nb_Wedges_;}
	int Nb_Pyramids(){return nb_Pyramids_;}
	int Nb_Tot_cell(){return Tot_cells_;}

	int Nb_plies(){return nb_plies_;}

	//~Maillage();
private:
	int nb_vertices_;
	int nb_Tetrahedra_;
	int nb_Triangles_;
	int nb_Quadrilaterals_;
	int nb_Hexahedra_;
	int nb_Cohesive_;
	int nb_Wedges_;
	int nb_Pyramids_;
	int Tot_cells_;

	int nb_plies_;

};

// ~~~~~~~~~~~~~~ INIT ~~~~~~~~~~~~~~

void Maillage::initialise() {
}

// ~~~~~~~~~~~~~~ PRINT ~~~~~~~~~~~~~~

void Maillage::print_hex(int number) {
	for(int n = 0; n < Hexahedra.rows()-1; n++) {
		std::cout << n << " : ";
		std::cout << vertices(0, Hexahedra(n, number)) << ", ";
		std::cout << vertices(1, Hexahedra(n, number)) << ", ";
		std::cout << vertices(2, Hexahedra(n, number)) << std::endl;
	}
}

// ~~~~~~~~~~~~~~ READ ~~~~~~~~~~~~~~

void Maillage::lecture_mesh(const std::string& filename) {
	std::ifstream input;
	input.open(filename, std::ios::in);
	if (!input.is_open()) {
	std::cout << "Error: Cannot open file" << std::endl;
	}
	std::cout << "Reading " << filename << std::endl;

	int en_tete_size = 4;
	if(mesh_type=="cgal") {en_tete_size = 3;}
	//En-tete
	for(int i=0; i<en_tete_size;++i) {
		std::string ligne {};
		std::getline(input, ligne);
	}

	// Vertices
	input >> nb_vertices_;
	vertices.resize(3, nb_vertices_);
	//std::cout << nb_vertices_ << std::endl;
	for(int i=0; i< nb_vertices_; ++i) {
		float x, y, z, d;
		input >> x >> y >> z >> d;
		vertices(0,i)=x;
		vertices(1,i)=y;
		vertices(2,i)=z;
		//~ if (i<10){std::cout << vertices(0,i) << ", " << vertices(1,i) << ", " << vertices(2,i) << std::endl;}
	}
	std::string ligne {};
	// on met deux lectures de ligne, car le cpp lit la fin de la derniere ligne des vertices
	std::getline(input, ligne);
	std::getline(input, ligne);

	//~ std::cout << ligne << std::endl;

	if (ligne.find("Triangles") != std::string::npos) {
		input >> nb_Triangles_;
		Triangles.resize(4, nb_Triangles_);

		for(int i=0; i< nb_Triangles_; ++i) {
			int a, b, c, marker;
			input >> a >> b >> c >> marker;
			Triangles(0,i)=a-1;
			Triangles(1,i)=b-1;
			Triangles(2,i)=c-1;
			Triangles(3,i)=marker;
		}
		std::getline(input, ligne);
		std::getline(input, ligne);
	} else {
		nb_Triangles_=0;
		Triangles.resize(4, 1);
	}

	//~ std::cout << ligne << std::endl;

	if (ligne.find("Quadrilaterals") != std::string::npos) {
		input >> nb_Quadrilaterals_;
		Quadrilaterals.resize(5, nb_Quadrilaterals_);

		for(int i=0; i< nb_Quadrilaterals_; ++i) {
			int a, b, c, d, marker;
			input >> a >> b >> c >> d >> marker;
			Quadrilaterals(0,i)=a-1;
			Quadrilaterals(1,i)=b-1;
			Quadrilaterals(2,i)=c-1;
			Quadrilaterals(3,i)=d-1;
			Quadrilaterals(4,i)=marker;
		}
		std::getline(input, ligne);
		std::getline(input, ligne);
	} else {
		nb_Quadrilaterals_=0;
		Quadrilaterals.resize(4, 1);
	}

	//~ std::cout << ligne << std::endl;

	if (ligne.find("Tetrahedra") != std::string::npos) {
		input >> nb_Tetrahedra_;
		Tetrahedra.resize(5, nb_Tetrahedra_);

		for(int i=0; i < nb_Tetrahedra_; ++i) {
			int a, b, c, d, type;
			input >> a >> b >> c >> d >> type;
			Tetrahedra(0,i)=a-1;
			Tetrahedra(1,i)=b-1;
			Tetrahedra(2,i)=c-1;
			Tetrahedra(3,i)=d-1;
			Tetrahedra(4,i)=type;
			//std::cout << Tetrahedra(0,i) << ", " << Tetrahedra(1,i) << ", " << Tetrahedra(2,i) << ", " << Tetrahedra(3,i) << ", " << Tetrahedra(4,i) << std::endl;
		}
		std::getline(input, ligne);
		std::getline(input, ligne);
	} else {
		nb_Tetrahedra_=0;
		Tetrahedra.resize(5, 1);
	}

	//~ std::cout << ligne << std::endl;

	if (ligne.find("Hexahedra") != std::string::npos) {
		input >> nb_Hexahedra_;
		Hexahedra.resize(9, nb_Hexahedra_);

		for(int i=0; i < nb_Hexahedra_; ++i) {
			int a, b, c, d, e, f, g, h, type;
			input >> a >> b >> c >> d >> e >> f >> g >> h >> type;
			Hexahedra(0,i)=a-1;
			Hexahedra(1,i)=b-1;
			Hexahedra(2,i)=c-1;
			Hexahedra(3,i)=d-1;
			Hexahedra(4,i)=e-1;
			Hexahedra(5,i)=f-1;
			Hexahedra(6,i)=g-1;
			Hexahedra(7,i)=h-1;
			Hexahedra(8,i)=type;
		}
		std::getline(input, ligne);
		std::getline(input, ligne);
	} else {
		nb_Hexahedra_=0;
		Hexahedra.resize(9, 1);
	}
	//~ std::cout << ligne << std::endl;

}

void Maillage::lecture_msh(const std::string& filename, bool only_3D, int cz_id) {
	std::ifstream input;
	input.open(filename, std::ios::in);
	if (!input.is_open()) {
	std::cout << "Error: Cannot open file" << std::endl;
	}
	std::cout << "Reading " << filename << std::endl;

	std::string ligne {};
	std::getline(input, ligne); // $MeshFormat

	float version;
	input >> version;
	//~ std::cout << version << std::endl;

	if (version >= 2.3 || version <= 2.1) {
		std::cout << "Error: Wrong format of the msh file, 2.2 is needed" << std::endl;
		return;
	}
	std::getline(input, ligne); // End of the format ligne
	std::getline(input, ligne); // $EndMeshFormat

	std::getline(input, ligne); // $Nodes
	input >> nb_vertices_;

	vertices.resize(3, nb_vertices_);
	for(int i=0; i< nb_vertices_; ++i) {
		float x, y, z, d;
		input >> d >> x >> y >> z;
		vertices(0,i)=x;
		vertices(1,i)=y;
		vertices(2,i)=z;
		//~ if (i<10){std::cout << vertices(0,i) << ", " << vertices(1,i) << ", " << vertices(2,i) << std::endl;}
	}
	std::getline(input, ligne); // End last node line
	std::getline(input, ligne); // $EndNodes

	std::getline(input, ligne); // $Elements
	input >> Tot_cells_;

	std::getline(input, ligne); // End tot cell line
	std::getline(input, ligne);

	std::vector<std::vector<int>> ltri, lquad, ltet, lhex, lwed, lpyr, lcoh;

	while(ligne!="$EndElements") {
		std::vector<int> line;
		line = splitngetINT(ligne);

		if (line[1]==2) {ltri.push_back(line);}
		if (line[1]==3) {lquad.push_back(line);}
		if (line[1]==4) {ltet.push_back(line);}
		if (line[1]==5) {
			if (line[3]==cz_id && cz_id>=0) {lcoh.push_back(line);}
			else {lhex.push_back(line);}
		}
		if (line[1]==6) {lwed.push_back(line);}
		if (line[1]==7) {lpyr.push_back(line);}

		std::getline(input, ligne);
	}

	if (only_3D==false) {
		nb_Triangles_=ltri.size();
		Triangles.resize(4, nb_Triangles_);

		for(int i = 0; i < nb_Triangles_; i++) {
			int beg_node = 2+ltri[i][2]+1;
			Triangles(0,i) = ltri[i][beg_node]-1;
			Triangles(1,i) = ltri[i][beg_node+1]-1;
			Triangles(2,i) = ltri[i][beg_node+2]-1;
			Triangles(3,i) = ltri[i][3];

			if (Triangles(3,i)>nb_nset) {nb_nset=Triangles(3,i);}
		}

		nb_Quadrilaterals_=lquad.size();
		Quadrilaterals.resize(5, nb_Quadrilaterals_);
		for(int i = 0; i < nb_Quadrilaterals_; i++) {
			int beg_node = 2+lquad[i][2]+1;
			Quadrilaterals(0,i) = lquad[i][beg_node]-1;
			Quadrilaterals(1,i) = lquad[i][beg_node+1]-1;
			Quadrilaterals(2,i) = lquad[i][beg_node+2]-1;
			Quadrilaterals(3,i) = lquad[i][beg_node+3]-1;
			Quadrilaterals(4,i) = lquad[i][3];

			if (Quadrilaterals(4,i)>nb_nset) {nb_nset=Quadrilaterals(4,i);}
		}
	} else {
		Tot_cells_-=(ltri.size()+lquad.size());
		nb_Triangles_=0;
		nb_Quadrilaterals_=0;
	}

	nb_Tetrahedra_=ltet.size();
	Tetrahedra.resize(5, nb_Tetrahedra_);
	for(int i = 0; i < nb_Tetrahedra_; i++) {
		int beg_node = 2+ltet[i][2]+1;
		Tetrahedra(0,i) = ltet[i][beg_node]-1;
		Tetrahedra(1,i) = ltet[i][beg_node+1]-1;
		Tetrahedra(2,i) = ltet[i][beg_node+2]-1;
		Tetrahedra(3,i) = ltet[i][beg_node+3]-1;
		Tetrahedra(4,i) = ltet[i][3];

		if (Tetrahedra(4,i)>nb_elset) {nb_elset=Tetrahedra(4,i);}
	}

	nb_Wedges_=lwed.size();
	Wedges.resize(7, nb_Wedges_);
	for(int i = 0; i < nb_Wedges_; i++) {
		int beg_node = 2+lwed[i][2]+1;
		Wedges(0,i) = lwed[i][beg_node]-1;
		Wedges(1,i) = lwed[i][beg_node+1]-1;
		Wedges(2,i) = lwed[i][beg_node+2]-1;
		Wedges(3,i) = lwed[i][beg_node+3]-1;
		Wedges(4,i) = lwed[i][beg_node+4]-1;
		Wedges(5,i) = lwed[i][beg_node+5]-1;
		Wedges(6,i) = lwed[i][3];

		if (Wedges(6,i)>nb_elset) {nb_elset=Wedges(6,i);}
	}

	nb_Hexahedra_=lhex.size();
	Hexahedra.resize(9, nb_Hexahedra_);
	for(int i = 0; i < nb_Hexahedra_; i++) {
		int beg_node = 2+lhex[i][2]+1;
		Hexahedra(0,i) = lhex[i][beg_node]-1;
		Hexahedra(1,i) = lhex[i][beg_node+1]-1;
		Hexahedra(2,i) = lhex[i][beg_node+2]-1;
		Hexahedra(3,i) = lhex[i][beg_node+3]-1;
		Hexahedra(4,i) = lhex[i][beg_node+4]-1;
		Hexahedra(5,i) = lhex[i][beg_node+5]-1;
		Hexahedra(6,i) = lhex[i][beg_node+6]-1;
		Hexahedra(7,i) = lhex[i][beg_node+7]-1;
		Hexahedra(8,i) = lhex[i][3];

		if (Hexahedra(8,i)>nb_elset) {nb_elset=Hexahedra(8,i);}
	}

	nb_Cohesive_=lcoh.size();
	Cohesive.resize(9, nb_Cohesive_);
	for(int i = 0; i < nb_Cohesive_; i++) {
		int beg_node = 2+lcoh[i][2]+1;
		Cohesive(0,i) = lcoh[i][beg_node]-1;
		Cohesive(1,i) = lcoh[i][beg_node+1]-1;
		Cohesive(2,i) = lcoh[i][beg_node+2]-1;
		Cohesive(3,i) = lcoh[i][beg_node+3]-1;
		Cohesive(4,i) = lcoh[i][beg_node+4]-1;
		Cohesive(5,i) = lcoh[i][beg_node+5]-1;
		Cohesive(6,i) = lcoh[i][beg_node+6]-1;
		Cohesive(7,i) = lcoh[i][beg_node+7]-1;
		Cohesive(8,i) = lcoh[i][3];

		if (Cohesive(8,i)>nb_elset) {nb_elset=Cohesive(8,i);}
	}

	nb_Pyramids_=lpyr.size();
	Pyramids.resize(6, nb_Pyramids_);
	for(int i = 0; i < nb_Pyramids_; i++) {
		int beg_node = 2+lpyr[i][2]+1;
		Pyramids(0,i) = lpyr[i][beg_node]-1;
		Pyramids(1,i) = lpyr[i][beg_node+1]-1;
		Pyramids(2,i) = lpyr[i][beg_node+2]-1;
		Pyramids(3,i) = lpyr[i][beg_node+3]-1;
		Pyramids(4,i) = lpyr[i][beg_node+4]-1;
		Pyramids(5,i) = lpyr[i][3];

		if (Pyramids(5,i)>nb_elset) {nb_elset=Pyramids(5,i);}
	}

	nb_nset=std::max(nb_nset-nb_elset, 0);

	std::cout << nb_vertices_       << " nodes" << std::endl;
	std::cout << nb_Triangles_      << " triangles" << std::endl;
	std::cout << nb_Quadrilaterals_ << " quadrilaterals"  << std::endl;
	std::cout << nb_Tetrahedra_     << " tetrahedra"  << std::endl;
	std::cout << nb_Hexahedra_      << " hexahedra"  << std::endl;
	std::cout << nb_Cohesive_       << " cohesive hexahedra"  << std::endl;
	std::cout << nb_Wedges_         << " wedges"  << std::endl;
	std::cout << nb_Pyramids_       << " pyramids"  << std::endl;
	std::cout << nb_nset            << " NSETs"  << std::endl;
	std::cout << nb_elset           << " ELSETs"  << std::endl;
	std::cout << std::endl;
}

void Maillage::lecture_points(const std::string& filename) {
	std::ifstream input;
	input.open(filename, std::ios::in);
	if (!input.is_open()) {
	std::cout << "Error: Cannot open file" << std::endl;
	}
	std::cout << "Reading " << filename << std::endl;

	input >> P1(0) >> P1(1);
	input >> P2(0) >> P2(1);
	input >> P3(0) >> P3(1);

}

void Maillage::lecture_fields(const std::string& filename) {
	std::ifstream input;
	input.open(filename, std::ios::in);
	if (!input.is_open()) {
	std::cout << "Error: Cannot open file" << std::endl;
	}
	std::cout << "Reading " << filename << std::endl;

	std::string ligne {};
	std::getline(input, ligne); // Commentaires

	int Tot_3D_elem = Tot_cells_-nb_Triangles_-nb_Quadrilaterals_;// Only 3D elements

	E.resize(6,Tot_3D_elem);
	S.resize(6,Tot_3D_elem);

	SDEG.resize(1,Tot_3D_elem);
	QUADSCRT.resize(1,Tot_3D_elem);

	for(int i = 0; i < Tot_3D_elem; i++) {
		int num;
		input >> num >> S(0, i) >> S(1, i) >> S(2, i) >> S(3, i) >> S(4, i) >> S(5, i)\
		>> E(0, i) >> E(1, i) >> E(2, i) >> E(3, i) >> E(4, i) >> E(5, i) >> SDEG(0,i)\
		>> QUADSCRT(0,i);
	}

}

void Maillage::lecture_node_fields(const std::string& filename) {
	std::ifstream input;
	input.open(filename, std::ios::in);
	if (!input.is_open()) {
	std::cout << "Error: Cannot open file" << std::endl;
	}
	std::cout << "Reading " << filename << std::endl;

	std::string ligne {};
	std::getline(input, ligne); // Commentaires

	U.resize(3,nb_vertices_);

	for(int i = 0; i < nb_vertices_; i++) {
		int num;
		input >> num >> U(0, i) >> U(1, i) >> U(2, i);
	}

}

// ~~~~~~~~~~~~~~ COPY | REMOVE | RESTORE ~~~~~~~~~~~~~~

void Maillage::copy(Maillage& m) {

	nb_vertices_= m.Nb_vertices();
	nb_Tetrahedra_ = m.Nb_Tetrahedra();
	nb_Triangles_ = m.Nb_Triangles();
	nb_Quadrilaterals_ = m.Nb_Quadrilaterals();
	nb_Hexahedra_ = m.Nb_Hexahedra();
	nb_Wedges_ = m.Nb_Wedges();
	nb_Pyramids_ = m.Nb_Pyramids();

	Tot_cells_ = m.Nb_Tot_cell();

	vertices.resize(3, nb_vertices_);

	for(int i = 0; i < nb_vertices_; i++) {
		for (int j = 0; j < 3; j++) {
			vertices(j,i) = m.vertices(j,i);
		}
	}

	Triangles.resize(4, nb_Triangles_);

	for(int i = 0; i < nb_Triangles_; i++) {
		for (int j = 0; j < 4; j++) {
			Triangles(j,i) = m.Triangles(j,i);
		}
	}

	Tetrahedra.resize(5, nb_Tetrahedra_);

	for(int i = 0; i < nb_Tetrahedra_; i++) {
		for (int j = 0; j < 5; j++) {
			Tetrahedra(j,i) = m.Tetrahedra(j,i);
		}
	}

	Quadrilaterals.resize(5, nb_Quadrilaterals_);

	for(int i = 0; i < nb_Quadrilaterals_; i++) {
		for (int j = 0; j < 5; j++) {
			Quadrilaterals(j,i) = m.Quadrilaterals(j,i);
		}
	}

	Wedges.resize(7, nb_Wedges_);

	for(int i = 0; i < nb_Wedges_; i++) {
		for (int j = 0; j < 7; j++) {
			Wedges(j,i) = m.Wedges(j,i);
		}
	}

	Hexahedra.resize(9, nb_Hexahedra_);

	for(int i = 0; i < nb_Hexahedra_; i++) {
		for (int j = 0; j < 9; j++) {
			Hexahedra(j,i) = m.Hexahedra(j,i);
		}
	}

	Pyramids.resize(6, nb_Pyramids_);

	for(int i = 0; i < nb_Pyramids_; i++) {
		for (int j = 0; j < 6; j++) {
			Pyramids(j,i) = m.Pyramids(j,i);
		}
	}
}

void Maillage::duplicates() {

	std::vector<int> dupl(nb_vertices_);
	dupl[0]=-1;

	int count=0;

	for(int i=1;i<nb_vertices_;++i) {
		dupl[i]=-1;
		for(int j=0;j<i;++j) {

			Vector3f u, v;
			u=vertices.col(i);
			v=vertices.col(j);

			if (issame(u,v)){
				dupl[i] = j;
				count+=1;
				break;
			}
		}
	}

	std::cout << "There are " << count << " duplicated nodes." << std::endl;
	std::cout << std::endl;

	for(int i = 0; i < nb_Wedges_; i++) {
		for (int j = 0; j < 6; j++) {
			if (dupl[Wedges(j,i)]>0) {Wedges(j,i) = dupl[Wedges(j,i)];}
		}
	}
}


// ~~~~~~~~~~~~~~ DUMP ~~~~~~~~~~~~~~

void Maillage::ecriture_mesh(const std::string& filename) {
	std::ofstream output;
	output.open(filename+".mesh", std::ios::out);
	if (!output.is_open()) {
	std::cout << "Error: Cannot open file" << std::endl;
	}
	std::cout << "Writting " << filename+".mesh" << std::endl;

	// En-tete
	output << "MeshVersionFormatted 2" << std::endl;
	output << "Dimension" << std::endl;
	output << "3" << std::endl;

	// VERTICES
	output << "Vertices" << std::endl;
	output << nb_vertices_ << std::endl;
	for(int i=0; i< nb_vertices_; ++i) {
		output << vertices(0,i) << " " << vertices(1,i) << " " << vertices(2,i) << " " << vertices(3,i) << std::endl;
	}

	if (nb_Triangles_>0) {
		// Triangles
		output << "Triangles" << std::endl;
		output << nb_Triangles_ << std::endl;
		for(int i=0; i< nb_Triangles_; ++i) {
			output << Triangles(0,i)+1 << " " << Triangles(1,i)+1 << " " << Triangles(2,i)+1 << " ";
			output << Triangles(3,i) << std::endl;
		}
	}

	if (nb_Quadrilaterals_>0) {
		// Quadrilaterals
		output << "Quadrilaterals" << std::endl;
		output << nb_Quadrilaterals_ << std::endl;
		for(int i=0; i< nb_Quadrilaterals_; ++i) {
			output << Quadrilaterals(0,i)+1 << " " << Quadrilaterals(1,i)+1 << " " << Quadrilaterals(2,i)+1 << " " << Quadrilaterals(3,i)+1 << " ";
			output << Quadrilaterals(4,i) << std::endl;
		}
	}

	if (nb_Tetrahedra_>0) {
		// Tetrahedra
		output << "Tetrahedra" << std::endl;
		output << nb_Tetrahedra_ << std::endl;
		for(int i=0; i< nb_Tetrahedra_; ++i) {
			output << Tetrahedra(0,i)+1 << " " << Tetrahedra(1,i)+1 << " " << Tetrahedra(2,i)+1 << " " << Tetrahedra(3,i)+1 << " ";
			output << Tetrahedra(4,i) << std::endl;
		}
	}

	if (nb_Hexahedra_>0) {
		// Hexahedra
		output << "Hexahedra" << std::endl;
		output << nb_Hexahedra_ << std::endl;
		for(int i=0; i< nb_Hexahedra_; ++i) {
			output << Hexahedra(0,i)+1 << " " << Hexahedra(1,i)+1 << " " << Hexahedra(2,i)+1 << " " << Hexahedra(3,i)+1 << " ";
			output << Hexahedra(4,i)+1 << " " << Hexahedra(5,i)+1 << " " << Hexahedra(6,i)+1 << " " << Hexahedra(7,i)+1 << " ";
			output << Hexahedra(8,i) << std::endl;
		}
	}
	output << "End" << std::endl;
}

void Maillage::ecriture_msh(const std::string& filename) {
	std::ofstream output;
	output.open(filename+".msh", std::ios::out);
	if (!output.is_open()) {
	std::cout << "Error: Cannot open file" << std::endl;
	}
	std::cout << "Writting " << filename+".msh" << std::endl;

	std::cout << nb_vertices_       << " nodes" << std::endl;
	std::cout << nb_Triangles_      << " triangles" << std::endl;
	std::cout << nb_Quadrilaterals_ << " quadrilaterals"  << std::endl;
	std::cout << nb_Tetrahedra_     << " tetrahedra"  << std::endl;
	std::cout << nb_Hexahedra_      << " hexahedra"  << std::endl;
	std::cout << nb_Wedges_         << " wedges"  << std::endl;
	std::cout << nb_Pyramids_       << " pyramids"  << std::endl;

	// En-tete
	output << "$MeshFormat" << std::endl;
	output << "2.2 0 8" << std::endl;
	output << "$EndMeshFormat" << std::endl;

	// VERTICES
	output << "$Nodes" << std::endl;
	output << nb_vertices_ << std::endl;
	for(int i=0; i< nb_vertices_; ++i) {
		output << i+1 << " " << vertices(0,i) << " " << vertices(1,i) << " " << vertices(2,i) << std::endl;
	}

	output << "$EndNodes" << std::endl;
	output << "$Elements" << std::endl;
	output << Tot_cells_ << std::endl;
	int count = 1;
	for (int i=0; i< nb_Triangles_; i++) {
		output << count << " "  << 2 << " "  << 1 << " " << Triangles(3,i) << " " ;
		output << Triangles(0,i)+1 << " " << Triangles(1,i)+1 << " " << Triangles(2,i)+1 << std::endl;
		count+=1;
	}
	for (int i=0; i< nb_Quadrilaterals_; i++) {
		output << count << " "  << 3 << " "  << 1 << " " << Quadrilaterals(4,i) << " " ;
		output << Quadrilaterals(0,i)+1 << " " << Quadrilaterals(1,i)+1 << " " << Quadrilaterals(2,i)+1 << " "  << Quadrilaterals(3,i)+1 << std::endl;
		count+=1;
	}
	for (int i=0; i< nb_Tetrahedra_; i++) {
		output << count << " "  << 4 << " "  << 1 << " " << Tetrahedra(4,i) << " " ;
		output << Tetrahedra(0,i)+1 << " " << Tetrahedra(1,i)+1 << " " << Tetrahedra(2,i)+1 << " "  << Tetrahedra(3,i)+1 << std::endl;
		count+=1;
	}
	for (int i=0; i< nb_Hexahedra_; i++) {
		output << count << " "  << 5 << " "  << 1 << " " << Hexahedra(8,i) << " " ;
		//~ output << count << " "  << 5 << " "  << 1 << " " << vol_marker << " " ;
		output << Hexahedra(0,i)+1 << " " << Hexahedra(1,i)+1 << " " << Hexahedra(2,i)+1 << " " << Hexahedra(3,i)+1 << " ";
		output << Hexahedra(4,i)+1 << " " << Hexahedra(5,i)+1 << " " << Hexahedra(6,i)+1 << " " << Hexahedra(7,i)+1 << std::endl;
		count+=1;
	}
	for (int i=0; i< nb_Wedges_; i++) {
		output << count << " "  << 6 << " "  << 1 << " " << Wedges(6,i) << " " ;
		//~ output << count << " "  << 6 << " "  << 1 << " " << vol_marker << " " ;
		output << Wedges(0,i)+1 << " " << Wedges(1,i)+1 << " " << Wedges(2,i)+1 << " "  << Wedges(3,i)+1 << " " << Wedges(4,i)+1 << " "  << Wedges(5,i)+1 << std::endl;
		count+=1;
	}
	for (int i=0; i< nb_Pyramids_; i++) {
		output << count << " "  << 7 << " "  << 1 << " " << Pyramids(5,i) << " " ;
		output << Pyramids(0,i)+1 << " " << Pyramids(1,i)+1 << " " << Pyramids(2,i)+1 << " "  << Pyramids(3,i)+1 << " "  << Pyramids(4,i)+1 << std::endl;
		count+=1;
	}
	output << "$EndElements" << std::endl;
}

void Maillage::ecriture_vtk(const std::string& filename) {

	std::string cell_type_tri = "triangle";
	int cell_type_tri_id = 5;

	std::string cell_type_quad = "quadrilater";
	int cell_type_quad_id = 9;

	std::string cell_type_tet = "tetrahedron";
	int cell_type_tet_id = 10;

	std::string cell_type_hex = "hexahedron";
	int cell_type_hex_id = 12;

	std::string cell_type_wed = "wedge";
	int cell_type_wed_id = 13;

	std::string cell_type_pyr = "pyramid";
	int cell_type_pyr_id = 14;

	std::string outfile(filename+".vtk");
	std::ofstream output;
	output.open(outfile.c_str());
	std::cout << "Writting " << outfile << std::endl;

	output << "# vtk DataFile Version 3.0\n";
	output << "Generated by Jean Bénézech\n";
	output << "ASCII\n";
	output << "DATASET UNSTRUCTURED_GRID\n";

	output << "POINTS " << nb_vertices_ << " float" << std::endl;
	for (int i=0; i< nb_vertices_; i++) {
		output << vertices(0,i) << " " << vertices(1,i) << " " << vertices(2,i) << std::endl;
	}

	//~ Tot_cells_ = nb_Hexahedra_;

	int second_number = (3+1)*nb_Triangles_ + (4+1)*nb_Quadrilaterals_ + (4+1)*nb_Tetrahedra_ + (5+1)*nb_Pyramids_ + (6+1)*nb_Wedges_ + (8+1)*nb_Hexahedra_;
	//~ int second_number = (8+1)*nb_Hexahedra_;

	output << "CELLS " << Tot_cells_ << " " << second_number << std::endl;
	for (int i=0; i< nb_Triangles_; i++) {
		output << 3 << " " << Triangles(0,i) << " " << Triangles(1,i) << " " << Triangles(2,i) << std::endl;
	}
	for (int i=0; i< nb_Quadrilaterals_; i++) {
		output << 4 << " " << Quadrilaterals(0,i) << " " << Quadrilaterals(1,i) << " " << Quadrilaterals(2,i) << " "  << Quadrilaterals(3,i) << std::endl;
	}
	for (int i=0; i< nb_Tetrahedra_; i++) {
		output << 4 << " " << Tetrahedra(0,i) << " " << Tetrahedra(1,i) << " " << Tetrahedra(2,i) << " "  << Tetrahedra(3,i) << std::endl;
	}
	for (int i=0; i< nb_Pyramids_; i++) {
		output << 5 << " " << Pyramids(0,i) << " " << Pyramids(1,i) << " " << Pyramids(2,i) << " "  << Pyramids(3,i) << " "  << Pyramids(4,i) << std::endl;
	}
	for (int i=0; i< nb_Wedges_; i++) {
		output << 6 << " " << Wedges(0,i) << " " << Wedges(1,i) << " " << Wedges(2,i) << " "  << Wedges(3,i) << " ";
		output << Wedges(4,i) << " "  << Wedges(5,i) << std::endl;
	}
	for (int i=0; i< nb_Hexahedra_; i++) {
		output << 8 << " " << Hexahedra(0,i) << " " << Hexahedra(1,i) << " " << Hexahedra(2,i) << " " << Hexahedra(3,i) << " ";
		output << Hexahedra(4,i) << " " << Hexahedra(5,i) << " " << Hexahedra(6,i) << " " << Hexahedra(7,i) << std::endl;
	}

	output << "CELL_TYPES " << Tot_cells_ << std::endl;
	for (int i=0; i< nb_Triangles_; i++) {
		output << cell_type_tri_id << std::endl;
	}
	for (int i=0; i< nb_Quadrilaterals_; i++) {
		output << cell_type_quad_id << std::endl;
	}
	for (int i=0; i< nb_Tetrahedra_; i++) {
		output << cell_type_tet_id << std::endl;
	}
	for (int i=0; i< nb_Pyramids_; i++) {
		output << cell_type_pyr_id << std::endl;
	}
	for (int i=0; i< nb_Wedges_; i++) {
		output << cell_type_wed_id << std::endl;
	}
	for (int i=0; i< nb_Hexahedra_; i++) {
		output << cell_type_hex_id << std::endl;
	}

	output << "CELL_DATA " << Tot_cells_ << std::endl;
	output << "SCALARS marker int 1"<< std::endl;
	output << "LOOKUP_TABLE default"<< std::endl;
	for (int i=0; i< nb_Triangles_; i++) {
		output << Triangles(3,i) << std::endl;
	}
	for (int i=0; i< nb_Quadrilaterals_; i++) {
		output << Quadrilaterals(4,i) << std::endl;
	}
	for (int i=0; i< nb_Tetrahedra_; i++) {
		output << Tetrahedra(4,i) << std::endl;
	}
	for (int i=0; i< nb_Pyramids_; i++) {
		output << Pyramids(5,i) << std::endl;
	}
	for (int i=0; i< nb_Wedges_; i++) {
		output << Wedges(6,i) << std::endl;
	}
	for (int i=0; i< nb_Hexahedra_; i++) {
		output << Hexahedra(8,i) << std::endl;
	}

	output << "SCALARS dir float 3\n";
	output << "LOOKUP_TABLE default\n";
	for (int i=0; i< Tot_cells_; i++) {
		output << sortie(1, i) << " " << sortie(2, i) << " " << sortie(3, i) << std::endl;
	}
	output << "SCALARS dir2 float 3\n";
	output << "LOOKUP_TABLE default\n";
	for (int i=0; i< Tot_cells_; i++) {
		output << sortie(4, i) << " " << sortie(5, i) << " " << sortie(6, i) << std::endl;
	}
	output.close();
}

void Maillage::ecriture_fields_vtk(const std::string& filename) {

	std::string cell_type_tri = "triangle";
	int cell_type_tri_id = 5;

	std::string cell_type_quad = "quadrilater";
	int cell_type_quad_id = 9;

	std::string cell_type_tet = "tetrahedron";
	int cell_type_tet_id = 10;

	std::string cell_type_hex = "hexahedron";
	int cell_type_hex_id = 12;

	std::string cell_type_coh = "cohesive";
	int cell_type_coh_id = 12;

	std::string cell_type_wed = "wedge";
	int cell_type_wed_id = 13;

	std::string cell_type_pyr = "pyramid";
	int cell_type_pyr_id = 14;

	std::string outfile(filename+".vtk");
	std::ofstream output;
	output.open(outfile.c_str());
	std::cout << "Writting " << outfile << std::endl;

	output << "# vtk DataFile Version 3.0\n";
	output << "Generated by Jean Bénézech\n";
	output << "ASCII\n";
	output << "DATASET UNSTRUCTURED_GRID\n";

	output << "POINTS " << nb_vertices_ << " float" << std::endl;
	for (int i=0; i< nb_vertices_; i++) {
		output << vertices(0,i) << " " << vertices(1,i) << " " << vertices(2,i) << std::endl;
	}

	//~ Tot_cells_ = nb_Hexahedra_;
	int end_2D = nb_Triangles_+nb_Quadrilaterals_;
	int Tot_3D_elem = Tot_cells_-end_2D;// Only 3D elements


	//~ int second_number = (3+1)*nb_Triangles_ + (4+1)*nb_Quadrilaterals_ + (4+1)*nb_Tetrahedra_ + (5+1)*nb_Pyramids_ + (6+1)*nb_Wedges_ + (8+1)*nb_Hexahedra_;
	int second_number = (4+1)*nb_Tetrahedra_ + (5+1)*nb_Pyramids_ + (6+1)*nb_Wedges_ + (8+1)*(nb_Hexahedra_+nb_Cohesive_);

	output << "CELLS " << Tot_3D_elem << " " << second_number << std::endl;
	//~ for (int i=0; i< nb_Triangles_; i++) {
		//~ output << 3 << " " << Triangles(0,i) << " " << Triangles(1,i) << " " << Triangles(2,i) << std::endl;
	//~ }
	//~ for (int i=0; i< nb_Quadrilaterals_; i++) {
		//~ output << 4 << " " << Quadrilaterals(0,i) << " " << Quadrilaterals(1,i) << " " << Quadrilaterals(2,i) << " "  << Quadrilaterals(3,i) << std::endl;
	//~ }
	for (int i=0; i< nb_Tetrahedra_; i++) {
		output << 4 << " " << Tetrahedra(0,i) << " " << Tetrahedra(1,i) << " " << Tetrahedra(2,i) << " "  << Tetrahedra(3,i) << std::endl;
	}
	for (int i=0; i< nb_Hexahedra_; i++) {
		output << 8 << " " << Hexahedra(0,i) << " " << Hexahedra(1,i) << " " << Hexahedra(2,i) << " " << Hexahedra(3,i) << " ";
		output << Hexahedra(4,i) << " " << Hexahedra(5,i) << " " << Hexahedra(6,i) << " " << Hexahedra(7,i) << std::endl;
	}
	for (int i=0; i< nb_Cohesive_; i++) {
		output << 8 << " " << Cohesive(0,i) << " " << Cohesive(1,i) << " " << Cohesive(2,i) << " " << Cohesive(3,i) << " ";
		output << Cohesive(4,i) << " " << Cohesive(5,i) << " " << Cohesive(6,i) << " " << Cohesive(7,i) << std::endl;
	}
	for (int i=0; i< nb_Pyramids_; i++) {
		output << 5 << " " << Pyramids(0,i) << " " << Pyramids(1,i) << " " << Pyramids(2,i) << " "  << Pyramids(3,i) << " "  << Pyramids(4,i) << std::endl;
	}
	for (int i=0; i< nb_Wedges_; i++) {
		output << 6 << " " << Wedges(0,i) << " " << Wedges(1,i) << " " << Wedges(2,i) << " "  << Wedges(3,i) << " ";
		output << Wedges(4,i) << " "  << Wedges(5,i) << std::endl;
	}



	output << "CELL_TYPES " << Tot_3D_elem << std::endl;
	//~ for (int i=0; i< nb_Triangles_; i++) {
		//~ output << cell_type_tri_id << std::endl;
	//~ }
	//~ for (int i=0; i< nb_Quadrilaterals_; i++) {
		//~ output << cell_type_quad_id << std::endl;
	//~ }
	for (int i=0; i< nb_Tetrahedra_; i++) {
		output << cell_type_tet_id << std::endl;
	}
	for (int i=0; i< nb_Hexahedra_; i++) {
		output << cell_type_hex_id << std::endl;
	}
	for (int i=0; i< nb_Cohesive_; i++) {
		output << cell_type_coh_id << std::endl;
	}
	for (int i=0; i< nb_Pyramids_; i++) {
		output << cell_type_pyr_id << std::endl;
	}
	for (int i=0; i< nb_Wedges_; i++) {
		output << cell_type_wed_id << std::endl;
	}


	output << "CELL_DATA " << Tot_3D_elem << std::endl;
	output << "SCALARS marker int 1"<< std::endl;
	output << "LOOKUP_TABLE default"<< std::endl;
	//~ for (int i=0; i< nb_Triangles_; i++) {
		//~ output << Triangles(3,i) << std::endl;
	//~ }
	//~ for (int i=0; i< nb_Quadrilaterals_; i++) {
		//~ output << Quadrilaterals(4,i) << std::endl;
	//~ }
	for (int i=0; i< nb_Tetrahedra_; i++) {
		output << Tetrahedra(4,i) << std::endl;
	}
	for (int i=0; i< nb_Hexahedra_; i++) {
		output << Hexahedra(8,i) << std::endl;
	}
	for (int i=0; i< nb_Cohesive_; i++) {
		output << Cohesive(8,i) << std::endl;
	}
	for (int i=0; i< nb_Pyramids_; i++) {
		output << Pyramids(5,i) << std::endl;
	}
	for (int i=0; i< nb_Wedges_; i++) {
		output << Wedges(6,i) << std::endl;
	}


	output << "SCALARS S float 6\n";
	output << "LOOKUP_TABLE default\n";
	for (int i=0; i< Tot_3D_elem; i++) {
		output << S(0, i) << " " << S(1, i) << " " << S(2, i) << " " <<  S(3, i) << " " << S(4, i) << " " << S(5, i) << std::endl;
	}
	output << "SCALARS E float 6\n";
	output << "LOOKUP_TABLE default\n";
	for (int i=0; i< Tot_3D_elem; i++) {
		output << E(0, i) << " " << E(1, i) << " " << E(2, i) << " " <<  E(3, i) << " " << E(4, i) << " " << E(5, i) << std::endl;
	}

	output << "SCALARS SDEG float 1\n";
	output << "LOOKUP_TABLE default\n";
	for (int i=0; i< Tot_3D_elem; i++) {
		output << SDEG(0, i) << std::endl;
	}
	output << "SCALARS QUADSCRT float 1\n";
	output << "LOOKUP_TABLE default\n";
	for (int i=0; i< Tot_3D_elem; i++) {
		output << QUADSCRT(0, i) << std::endl;
	}


	output << "POINT_DATA " << nb_vertices_ << std::endl;
	output << "SCALARS U float 3\n";
	output << "LOOKUP_TABLE default\n";
	for (int i=0; i< nb_vertices_; i++) {
		output << U(0, i) << " " << U(1, i) << " " << U(2, i) << std::endl;
	}

	output.close();
}

void Maillage::ecriture_ori_txt(const std::string& filename) {
	std::ofstream output;
	output.open(filename, std::ios::out);
	if (!output.is_open()) {
	std::cout << "Error: Cannot open file" << std::endl;
	}
	std::cout << "Writting " << filename << std::endl;

	// En-tete
	output << Tot_cells_-(Nb_Quadrilaterals()+Nb_Triangles()) << std::endl;

	for(int i=(Nb_Quadrilaterals()+Nb_Triangles()); i< sortie.cols(); ++i) {
		output << i+1 << " ";
		for(int j=1; j< 12; ++j) {
			output << sortie(j,i) << " ";
		}
		output << std::endl;
	}
}

void Maillage::ecriture_ori_inp(const std::string& filename) {

	std::ofstream output;
	output.open(filename+"_ori.inp", std::ios::out);
	if (!output.is_open()) {
	std::cout << "Error: Cannot open file" << std::endl;
	}
	std::cout << "Writting " << filename+"_ori.inp" << std::endl;


	output << "*ORIENTATION, NAME=ori_glob" << std::endl;
	output << "1., 0., 0., 0., 1., 0." << std::endl;
	output << "1, 0." << std::endl;
	output << "*DISTRIBUTION, NAME=ori_loc_dist, LOCATION=ELEMENT, TABLE=ori_tab" << std::endl;
	output << ", 1.0,  0.0,  0.0,  0.0,  1.0, 0.0" << std::endl;

	int a = (Nb_Quadrilaterals()+Nb_Triangles());

	for(int i=a; i< sortie.cols(); ++i) {
		output << i+1-a ;
		for(int j=3; j< 9; ++j) {
			output << ", " << sortie(j,i);
		}
		output << std::endl;
	}
	output << "*ORIENTATION, NAME=ori_loc" << std::endl;
	output << "ori_loc_dist" << std::endl;
}

void Maillage::ecriture_inp(const std::string& filename) {
	std::ofstream output;
	output.open(filename+"_mesh.inp", std::ios::out);
	if (!output.is_open()) {
	std::cout << "Error: Cannot open file" << std::endl;
	} else {
		std::cout << "Writting " << filename+"_mesh.inp" << std::endl;
	}

	// ~~~~~~~~~ NSET and ElSET preparation ~~~~~~~~~ fonction des labels donnes sur gmsh
	std::vector<std::vector<int>> nsets, _tet, _quad, elsets;
	std::vector<Vector3f> masterNodes;
	nsets.resize(nb_nset);
	masterNodes.resize(nb_nset);
	_tet.resize(nb_nset);
	_quad.resize(nb_nset);


	//NSET
	for(int i=0; i< nb_Triangles_; ++i) {
		_tet[Triangles(3, i)-nb_elset-1].push_back(i);
	}
	for(int i=0; i< nb_Quadrilaterals_; ++i) {
		_quad[Quadrilaterals(4, i)-nb_elset-1].push_back(i);
	}

	for(int j=0; j< nb_nset;j++) {
		for(int i=0; i<_tet[j].size(); i++){
			for(int k=0; k<3; k++) {
				nsets[j].push_back(Triangles(k,_tet[j][i]));
			}
		}
		for(int i=0; i<_quad[j].size(); i++){
			for(int k=0; k<4; k++) {
				nsets[j].push_back(Quadrilaterals(k,_quad[j][i]));
			}
		}
		std::sort(nsets[j].begin(), nsets[j].end());
		auto last = std::unique(nsets[j].begin(), nsets[j].end());
		nsets[j].erase(last, nsets[j].end());

		masterNodes[j](0)=0.0;
		masterNodes[j](1)=0.0;
		masterNodes[j](2)=0.0;
		for(int i=0; i<nsets[j].size(); i++) {
			masterNodes[j](0) += vertices(0,nsets[j][i]);
			masterNodes[j](1) += vertices(1,nsets[j][i]);
			masterNodes[j](2) += vertices(2,nsets[j][i]);
		}
		masterNodes[j](0)/=nsets[j].size();
		masterNodes[j](1)/=nsets[j].size();
		masterNodes[j](2)/=nsets[j].size();

		//~ std::cout << "MasterNode " << j << " : " << masterNodes[j](0) << ", " << masterNodes[j](1) << ", " << masterNodes[j](2) << std::endl;
	}

	// ~~~~~~~~~ N O D E S ~~~~~~~~~

	output << "*NODE" << std::endl;
	output << nb_vertices_+4 << std::endl;
	for(int i=0; i< nb_vertices_; ++i) {
		output << i+1 << ", " << vertices(0,i) << ", " << vertices(1,i) << ", " << vertices(2,i) << std::endl;
	}
	// Ajout de master nodes
	for(int j=0; j< nb_nset;j++) {
		output << nb_vertices_+1+j << ", " <<  masterNodes[j](0) << ", " << masterNodes[j](1) << ", " << masterNodes[j](2) << std::endl;
	}

	// ~~~~~~~~~ E L E M E N T S ~~~~~~~~~

	int incr_tet, incr_hex, incr_wed;
	int cnt=1;
	output << "******* E L E M E N T S *************" << std::endl;
	// Tetrahedra
	incr_tet = cnt;
	if (nb_Tetrahedra_>0) {
		output << "*ELEMENT, type=C3D4, ELSET=Tetrahedra" << std::endl;
		for(int i=0; i< nb_Tetrahedra_; ++i) {
			output << cnt << ", " << Tetrahedra(0,i)+1 << ", " << Tetrahedra(1,i)+1 << ", " << Tetrahedra(2,i)+1 << ", " << Tetrahedra(3,i)+1 << std::endl;
			cnt+=1;
		}
	}
	// Hexahedra
	incr_hex = cnt;
	if (nb_Hexahedra_>0) {
		output << "*ELEMENT, type=C3D8, ELSET=Hexahedra" << std::endl;
		for(int i=0; i< nb_Hexahedra_; ++i) {
			output << cnt << ", " << Hexahedra(0,i)+1 << ", " << Hexahedra(1,i)+1 << ", " << Hexahedra(2,i)+1 << ", " << Hexahedra(3,i)+1 << ", ";
			output << Hexahedra(4,i)+1 << ", " << Hexahedra(5,i)+1 << ", " << Hexahedra(6,i)+1 << ", " << Hexahedra(7,i)+1 << std::endl;
			cnt+=1;
		}
	}
	// Wedges
	incr_wed = cnt;
	if (nb_Wedges_>0) {
		output << "*ELEMENT, type=C3D6, ELSET=Wedges" << std::endl;
		for(int i=0; i< nb_Wedges_; ++i) {
			output << cnt << ", " << Wedges(0,i)+1 << ", " << Wedges(1,i)+1 << ", " << Wedges(2,i)+1 << ", " << Wedges(3,i)+1 << ", ";
			output << Wedges(4,i)+1 << ", " << Wedges(5,i)+1 << std::endl;
			cnt+=1;
		}
	}

	//ELSET
	elsets.resize(nb_elset);
	for(int i=0; i< nb_Tetrahedra_; ++i) {
		elsets[Tetrahedra(4, i)-1].push_back(i+incr_tet);
	}
	for(int i=0; i< nb_Hexahedra_; ++i) {
		elsets[Hexahedra(4, i)-1].push_back(i+incr_hex);
	}
	for(int i=0; i< nb_Wedges_; ++i) {
		elsets[Wedges(4, i)-1].push_back(i+incr_wed);
	}

	// ~~~~~~~~~ E L S E T S ~~~~~~~~~
	output << "******* E L E M E N T S   S E T S *************" << std::endl;
	output << "*ELSET,ELSET=All_elements" << std::endl;
	for(int i=0; i< nb_Tetrahedra_; ++i) {
		output << i+incr_tet;
		if ((i+1)%10 == 0) { output << "," << std::endl; }
		else { output << ", "; }
	}
	for(int i=0; i< nb_Hexahedra_; ++i) {
		output << i+incr_hex;
		if ((i+1)%10 == 0) { output << "," << std::endl; }
		else { output << ", "; }
	}
	for(int i=0; i< nb_Wedges_; ++i) {
		output << i+incr_wed;
		if ((i+1)%10 == 0) { output << "," << std::endl; }
		else { output << ", "; }
	}
	output << std::endl;
	for(int j=0; j< nb_elset;j++) {
		output << "*ELSET,ELSET=elset" << j << std::endl;
		for(int i=0; i< elsets[j].size(); ++i) {
			output << elsets[j][i];
			if (i==elsets[j].size()-1) { output << "," << std::endl; }
			else if ((i+1)%10 == 0) { output << "," << std::endl; }
			else { output << ", "; }
		}
	}

	// ~~~~~~~~~ NSET ~~~~~~~~~
	output << "******* N O D E S   S E T S *************" << std::endl;
	for(int j=0; j< nb_nset;j++) {
		if (nsets[j].size()>0) {
			output << "*NSET,NSET=nset" << j << std::endl;
			for(int i=0; i< nsets[j].size(); ++i) {
				output << nsets[j][i]+1;
				if (i==nsets[j].size()-1) { output << "," << std::endl; }
				else if ((i+1)%10 == 0) { output << "," << std::endl; }
				else { output << ", "; }
			}
		}
	}

	for(int j=0; j< nb_nset;j++) {
		output << "*NSET,NSET=MasterNode" << j << std::endl;
		output << nb_vertices_+1+j << ", " << std::endl;
	}

}

#endif /* end of include guard: MAILLAGE_H */
