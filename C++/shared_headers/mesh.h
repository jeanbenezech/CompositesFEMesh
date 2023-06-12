#ifndef MESH_H
#define MESH_H

#include "utils.h"
#include "element.h"
#include "vertice.h"
#include "parameters.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <Eigen/Eigenvalues>
#include <Eigen/Dense>

using namespace Eigen;

class Mesh {
public:
	// PRINCIPAL VARIABLES
	std::vector<Vertice> Vertices;
	std::vector<Element> Elements;

	std::string mesh_type;

	int nb_nset=0;
	int nb_elset=0;

	bool isShell=false;

	std::vector<double> stacking_sequence;
	std::vector<double> ply_type;
	int nb_corner, type;
	std::vector<Vector2f> P;

	// Abaqus specific
	std::vector<std::vector<int>> nsets, elsets, Coh_elsets;
	std::vector<Vector3d> masterNodes;
	bool isNsetDone=false;


	// FUNCTION
	void read_mesh(const std::string& filename);
	void read_msh(const std::string& filename, bool only_3D, int cz_id);
	void read_points(const std::string& filename);
	void read_elem_fields(const std::string& filename);
	void read_elem_fields_stresses_strains_only(const std::string& filename);
	void read_node_fields(const std::string& filename);

	void write_mesh(const std::string& filename);
	void write_msh(const std::string& filename, int verbosity);
	bool exportDir=false;
	bool exportAbaqusFields=false;
	bool exportAbaqusDisplacement=false;
	bool exportRandomFieldOnElement=false;
	void write_vtk(const std::string& filename, int verbosity);
	void write_point_cloud_vtk(const std::string& filename, int verbosity);
	void write_inp(const std::string& filename);
	void write_ori_txt(const std::string& filename);
	void write_ori_inp(const std::string& filename);
	void write_abaqus_cae_input(const std::string& filename, Parameters& param);

	void extract_AbaqusSets();
	void elSets_delam(double locX, double locZ, double radius, double radius_offset_X,  double radius_offset_Z);

	void initialise(Parameters& param);
	void print(std::string type, int number);

	int Nb_vertices(){return nb_vertices_;}
	int Nb_Tot_cell(){return Tot_cells_;}

	void Set_vertices(int& value){nb_vertices_=value;}

	int Nb_plies(){return nb_plies_;}

	//~Mesh();
private:
	int nb_vertices_;
	int Tot_cells_;

	int nb_plies_;

};

// ~~~~~~~~~~~~~~ INIT ~~~~~~~~~~~~~~

void Mesh::initialise(Parameters& param) {

	nb_plies_=param.nb_plies;

	stacking_sequence.resize(nb_plies_);

	for(int i=0; i<nb_plies_; i++)
		stacking_sequence[i] = param.StackSeq[i](0);

}

// ~~~~~~~~~~~~~~ PRINT ~~~~~~~~~~~~~~

void Mesh::print(std::string type, int number) {
	for(auto& elem : Elements){
		if(elem.type != "Triangle" || elem.type != "Quadrilateral"){
			for(int n = 0; n < elem.Nodes.rows(); n++) {
				std::cout << n << " : ";
				std::cout << Vertices[elem.Nodes(n, number)].coord(0) << ", ";
				std::cout << Vertices[elem.Nodes(n, number)].coord(1) << ", ";
				std::cout << Vertices[elem.Nodes(n, number)].coord(2) << std::endl;
			}
		}
	}
}
// ~~~~~~~~~~~~~~ READ ~~~~~~~~~~~~~~

void Mesh::read_mesh(const std::string& filename) {
	std::ifstream input;
	input.open(filename, std::ios::in);
	if (!input.is_open()) {
	std::cout << "Error: Cannot open file " << std::endl;
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
	Vertices.resize(nb_vertices_);
	
	//std::cout << nb_vertices_ << std::endl;
	for(int i=0; i< nb_vertices_; ++i) {
		double x, y, z, d;
		input >> x >> y >> z >> d;
		Vertices[i].coord(0)=x;
		Vertices[i].coord(1)=y;
		Vertices[i].coord(2)=z;
	}
	std::string ligne {};

	std::getline(input, ligne);
	std::getline(input, ligne);

	if (ligne.find("Triangles") != std::string::npos) {
		Element tmp;
		int nb_elem;
		input >> nb_elem;
		tmp.initialise("Triangle",nb_elem,1);

		for(int i=0; i< nb_elem; ++i) {
			for(int j=0; j<tmp.size; j++)
				input >> tmp.Nodes(j,i);
			input >> tmp.Markers(0,i);
		}
		std::getline(input, ligne);
		std::getline(input, ligne);
		Elements.push_back(tmp);
	}

	if (ligne.find("Quadrilaterals") != std::string::npos) {
		Element tmp;
		int nb_elem;
		input >> nb_elem;
		tmp.initialise("Quadrilateral",nb_elem,1);

		for(int i=0; i< nb_elem; ++i) {
			for(int j=0; j<tmp.size; j++)
				input >> tmp.Nodes(j,i);
			input >> tmp.Markers(0,i);
		}
		std::getline(input, ligne);
		std::getline(input, ligne);
		Elements.push_back(tmp);
	}

	if (ligne.find("Tetrahedra") != std::string::npos) {
		Element tmp;
		int nb_elem;
		input >> nb_elem;
		tmp.initialise("Tetrahedron",nb_elem,1);

		for(int i=0; i< nb_elem; ++i) {
			for(int j=0; j<tmp.size; j++)
				input >> tmp.Nodes(j,i);
			input >> tmp.Markers(0,i);
		}
		std::getline(input, ligne);
		std::getline(input, ligne);
		Elements.push_back(tmp);
	}

	if (ligne.find("Hexahedra") != std::string::npos) {
		Element tmp;
		int nb_elem;
		input >> nb_elem;
		tmp.initialise("Hexahedron",nb_elem,1);

		for(int i=0; i< nb_elem; ++i) {
			for(int j=0; j<tmp.size; j++)
				input >> tmp.Nodes(j,i);
			input >> tmp.Markers(0,i);
		}
		std::getline(input, ligne);
		std::getline(input, ligne);
		Elements.push_back(tmp);
	}

	Tot_cells_ = 0;
	for(auto& elem : Elements)
		Tot_cells_+=elem.nb;

}

void Mesh::read_msh(const std::string& filename, bool only_3D, int cz_id) {
	std::ifstream input;
	input.open(filename, std::ios::in);
	if (!input.is_open()) {
	std::cout << "Error: Cannot open file " << filename  << std::endl;
	}
	std::cout << "Reading " << filename << std::endl;

	std::string ligne {};
	std::getline(input, ligne); // $MeshFormat

	double version;
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

	// vertice_parents.resize(nb_vertices_);
	// for(int i=0; i< nb_vertices_; ++i) {
	// 	vertice_parents[i].resize(0);
	// }

	Vertices.resize(nb_vertices_);

	for(int i=0; i< nb_vertices_; ++i) {
		double x, y, z, d;
		input >> d >> x >> y >> z;
		Vertices[d-1].coord(0)=x;
		Vertices[d-1].coord(1)=y;
		Vertices[d-1].coord(2)=z;
		//~ if (i<10){std::cout << Vertices[i].coord(0) << ", " << Vertices[i].coord(1) << ", " << Vertices[i].coord(2) << std::endl;}

		// initialize the normal to vertices at zero
		Vertices[d-1].normal = Vector3d::Zero();

	}
	std::getline(input, ligne); // End last node line
	std::getline(input, ligne); // $EndNodes

	std::getline(input, ligne); // $Elements
	input >> Tot_cells_;

	std::getline(input, ligne); // End tot cell line
	std::getline(input, ligne);

	std::vector<std::vector<int>> ltri, lquad, ltet, lhex, lwed, lpyr, lcoh;

	int nb_diff_type = 0;
	bool istri=false , isquad=false , istet=false , ishex=false , iswed=false , ispyr=false , iscoh=false;

	while(ligne!="$EndElements") {
		std::vector<int> line;
		line = splitngetINT(ligne);

		if (line[1]==2) {
			ltri.push_back(line);
			if(!istri){istri = true; nb_diff_type+=1;}
		}
		if (line[1]==3) {
			lquad.push_back(line);
			if(!isquad){isquad = true; nb_diff_type+=1;}
		}
		if (line[1]==4) {
			ltet.push_back(line);
			if(!istet){istet = true; nb_diff_type+=1;}
		}
		if (line[1]==5) {
			if (line[3]==cz_id && cz_id>=0) {
				lcoh.push_back(line);
				if(!iscoh){iscoh = true; nb_diff_type+=1;}
			} else {
				lhex.push_back(line);
				if(!ishex){ishex = true; nb_diff_type+=1;}
			}
		}
		if (line[1]==6) {
			lwed.push_back(line);
			if(!ispyr){ispyr = true; nb_diff_type+=1;}
		}
		if (line[1]==7) {
			lpyr.push_back(line);
			if(!ispyr){ispyr = true; nb_diff_type+=1;}
		}

		std::getline(input, ligne);
	}

	if (only_3D==true) {
		if(istri){nb_diff_type-=1;}
		if(isquad){nb_diff_type-=1;}
		Tot_cells_-=(ltri.size()+lquad.size());
	}


	Elements.resize(nb_diff_type);

	// std::cout << "Elements.size() : " << Elements.size() << std::endl;

	int incr_elem=0;
	if (only_3D==false) {
		if(istri){
			Elements[incr_elem].initialise("Triangle", ltri.size(), ltri[0][2]-1);
			Elements[incr_elem].fill(ltri, nb_nset);
			incr_elem += 1;
		}
		if(isquad){
			Elements[incr_elem].initialise("Quadrilateral", lquad.size(), lquad[0][2]-1);
			Elements[incr_elem].fill(lquad, nb_nset);
			incr_elem += 1;
		}
	}

	if(istet){
		Elements[incr_elem].initialise("Tetrahedron", ltet.size(), 1);//ltet[0][2]-1);
		Elements[incr_elem].fill(ltet, nb_elset);
		incr_elem += 1;
	}
	if(ishex){
		Elements[incr_elem].initialise("Hexahedron", lhex.size(), 1, isShell);//lhex[0][2]-1);
		Elements[incr_elem].fill(lhex, nb_elset);
		incr_elem += 1;
	}
	if(iscoh){
		Elements[incr_elem].initialise("Cohesive", lcoh.size(), 1);//lcoh[0][2]-1);
		Elements[incr_elem].fill(lcoh, nb_elset);
		incr_elem += 1;
	}
	if(iswed){
		Elements[incr_elem].initialise("Prism", lwed.size(), 1);//lwed[0][2]-1);
		Elements[incr_elem].fill(lwed, nb_elset);
		incr_elem += 1;
	}
	if(ispyr){
		Elements[incr_elem].initialise("Pyramid", lpyr.size(), 1);//lpyr[0][2]-1);
		Elements[incr_elem].fill(lpyr, nb_elset);
		incr_elem += 1;
	}

	// Initial barycenter
	for(auto& elem : Elements){
		for(int l=0; l < elem.nb; ++l) {
			elem.Initial_barycenter.col(l) = Vector3d::Zero();
			for(int i=0; i<elem.size; i++)
				elem.Initial_barycenter.col(l) += Vertices[elem.Nodes(i,l)].coord;
			elem.Initial_barycenter.col(l) *= 1.0 / (elem.size+0.0);
		}
	}

	// Vectices element_ids
	// for(auto& elem : Elements){
	// 	for(int l=0; l < elem.nb; ++l) {
	// 		for(int i=0; i<elem.size; i++){
	// 			Vector2i tmp;
	// 			tmp(0)=i;
	// 			tmp(1)=l;
	// 			Vertices[elem.Nodes(i,l)].element_ids.push_back(tmp);
	// 		}
	// 	}
	// }

}

void Mesh::read_points(const std::string& filename) {
	std::ifstream input;
	input.open(filename, std::ios::in);
	if (!input.is_open()) {
		std::cout << "Error: Cannot open file " << filename << std::endl;
	} else{
		std::cout << "Reading " << filename << std::endl;
	}

	input >> nb_corner;

	P.resize(nb_corner);
	for (int i=0; i<P.size();i++) {
		input >> P[i](0) >> P[i](1);
	}
}

void Mesh::read_elem_fields(const std::string& filename) {
	std::ifstream input;
	input.open(filename, std::ios::in);
	if (!input.is_open()) {
	std::cout << "Error: Cannot open file " << filename << std::endl;
	}
	std::cout << "Reading " << filename << std::endl;

	std::string ligne {};
	std::getline(input, ligne);

	Matrix<double, Dynamic, Dynamic> tmp;
	tmp.resize(14,Tot_cells_);

	for(int i = 0; i < Tot_cells_; i++) {
		int num;
		input >> num >> tmp(0, i) >> tmp(1, i) >> tmp(2, i) >> tmp(3, i) >> tmp(4, i) >> tmp(5, i)\
		>> tmp(6, i) >> tmp(7, i) >> tmp(8, i) >> tmp(9, i) >> tmp(10, i) >> tmp(11, i) >> tmp(12,i)\
		>> tmp(13,i);
	}

	for(auto& elem : Elements){
		for (int i=0; i< elem.nb; i++) {
			for(int j=0; j<6; j++){
				elem.S(j,i) = tmp(j,elem.global_indices[i]);
				elem.E(j,i) = tmp(j+6,elem.global_indices[i]);
			}
			elem.SDEG(0,i) = tmp(12,elem.global_indices[i]);
			elem.QUADSCRT(0,i) = tmp(13,elem.global_indices[i]);
		}
	}

}

void Mesh::read_elem_fields_stresses_strains_only(const std::string& filename) {
	std::ifstream input;
	input.open(filename, std::ios::in);
	if (!input.is_open()) {
	std::cout << "Error: Cannot open file " << filename << std::endl;
	}
	std::cout << "Reading " << filename << std::endl;

	std::string ligne {};
	std::getline(input, ligne);

	Matrix<double, Dynamic, Dynamic> tmp;
	tmp.resize(12,Tot_cells_);

	for(int i = 0; i < Tot_cells_; i++) {
		int num;
		input >> num >> tmp(0, i) >> tmp(1, i) >> tmp(2, i) >> tmp(3, i) >> tmp(4, i) >> tmp(5, i)\
		>> tmp(6, i) >> tmp(7, i) >> tmp(8, i) >> tmp(9, i) >> tmp(10, i) >> tmp(11, i);
	}

	for(auto& elem : Elements){
		for (int i=0; i< elem.nb; i++) {
			// std::cout << elem.S.rows() << "; " << elem.S.cols() << std::endl;
			// std::cout << elem.E.rows() << "; " << elem.E.cols() << std::endl;
			// std::cout <<  "i: " << elem.global_indices[-1] << std::endl;
			for(int j=0; j<6; j++){
				elem.S(j,i) = tmp(j,i);
				elem.E(j,i) = tmp(j+6,i);
				// elem.S(j,i) = tmp(j,elem.global_indices[i]-1);
				// elem.E(j,i) = tmp(j+6,elem.global_indices[i]-1);
			}
		}
	}

	// std::cout << "end of read element fields." << std::endl;

}

void Mesh::read_node_fields(const std::string& filename) {
	std::ifstream input;
	input.open(filename, std::ios::in);
	if (!input.is_open()) {
	std::cout << "Error: Cannot open file " << filename << std::endl;
	}
	std::cout << "Reading " << filename << std::endl;

	std::string ligne {};
	std::getline(input, ligne); // Commentaires

	for(int i = 0; i < nb_vertices_; i++) {
		int num;
		double ux, uy, uz;
		input >> num >> ux >> uy >> uz;
		if(abs(ux)>1e-8){Vertices[num-1].X(0) = ux;}else{Vertices[num-1].X(0) = 0.;}
		if(abs(uy)>1e-8){Vertices[num-1].X(1) = uy;}else{Vertices[num-1].X(1) = 0.;}
		if(abs(uz)>1e-8){Vertices[num-1].X(2) = uz;}else{Vertices[num-1].X(2) = 0.;}
	}

}

// ~~~~~~~~~~~~~~ DUMP ~~~~~~~~~~~~~~

void Mesh::write_mesh(const std::string& filename) {
	std::ofstream output;
	output.open(filename+".mesh", std::ios::out);
	if (!output.is_open()) {
	std::cout << "Error: Cannot open file " << filename << std::endl;
	}
	std::cout << "writing " << filename+".mesh" << std::endl;

	// En-tete
	output << "MeshVersionFormatted 2" << std::endl;
	output << "Dimension" << std::endl;
	output << "3" << std::endl;

	// VERTICES
	output << "Vertices" << std::endl;
	output << nb_vertices_ << std::endl;
	for(int i=0; i< nb_vertices_; i++) {
		output << Vertices[i].coord(0) << " " << Vertices[i].coord(1) << " " << Vertices[i].coord(2) << " " << i << std::endl;
	}

	for(auto& elem : Elements){
		if(elem.inria_mesh_name!="None"){
			output << elem.inria_mesh_name << std::endl;
			output << elem.nb << std::endl;
			for (int i=0; i< elem.nb; i++) {
				for (int j =0; j < elem.Nodes.rows(); j++)
					output << elem.Nodes(j,i)+1 << " ";
				for (int j =0; j < elem.Markers.rows()-1; j++)
					output << elem.Markers(j,i) << " ";
				output << elem.Markers(elem.Markers.rows()-1,i)+1;
				output << std::endl;
			}
		}
	}
	output << "End" << std::endl;
}

void Mesh::write_msh(const std::string& filename, int verbosity=1) {
	std::ofstream output;
	output.open(filename+".msh", std::ios::out);
	if (!output.is_open()) {
	std::cout << "Error: Cannot open file " << filename << std::endl;
	}
	std::cout << "writing " << filename+".msh" << std::endl;

	if(verbosity){
		std::cout << nb_vertices_ << " nodes" << std::endl;
		for(auto& elem : Elements){
			std::cout << elem.nb << " " << elem.type << std::endl;
		}
	}
	// En-tete
	output << "$MeshFormat" << std::endl;
	output << "2.2 0 8" << std::endl;
	output << "$EndMeshFormat" << std::endl;

	// VERTICES
	output << "$Nodes" << std::endl;
	output << nb_vertices_ << std::endl;
	for(int i=0; i< nb_vertices_; ++i) {
		output << std::setprecision(15) << i+1 << " " << Vertices[i].coord(0) << " " << Vertices[i].coord(1) << " " << Vertices[i].coord(2) << std::endl;
	}
	output << "$EndNodes" << std::endl;
	output << "$Elements" << std::endl;
	output << Tot_cells_ << std::endl;

	for(auto& elem : Elements){
		for (int i=0; i< elem.nb; i++) {
			output << elem.global_indices[i] << " " << elem.msh_type << " ";
			output << elem.Markers.rows() << " ";
			for (int j =0; j < elem.Markers.rows(); j++)
				output << elem.Markers(j,i) << " ";
			for (int j =0; j < elem.Nodes.rows()-1; j++)
				output << elem.Nodes(j,i)+1 << " ";
			output << elem.Nodes(elem.Nodes.rows()-1,i)+1;
			output << std::endl;
		}
	}
	output << "$EndElements" << std::endl;
}

void Mesh::write_vtk(const std::string& filename, int verbosity=0) {

	if(verbosity){
		std::cout << nb_vertices_ << " nodes" << std::endl;
		for(auto& elem : Elements){
			std::cout << elem.nb << " " << elem.type << std::endl;
		}
	}

	std::string outfile(filename+".vtk");
	std::ofstream output;
	output.open(outfile.c_str());
	std::cout << "writing " << outfile << std::endl;

	output << "# vtk DataFile Version 3.0\n";
	output << "Generated by Jean Bénézech\n";
	output << "ASCII\n";
	output << "DATASET UNSTRUCTURED_GRID\n";

	output << "POINTS " << nb_vertices_ << " double" << std::endl;
	for (int i=0; i< nb_vertices_; i++) {
		output << Vertices[i].coord(0) << " " << Vertices[i].coord(1) << " " << Vertices[i].coord(2) << std::endl;
	}

	int second_number = 0;
	for(auto& elem : Elements)
		second_number += (elem.size+1)*elem.nb;

	output << "CELLS " << Tot_cells_ << " " << second_number << std::endl;
	for(auto& elem : Elements){
		for (int i=0; i< elem.nb; i++) {
			output << elem.size << " ";
			for (int j =0; j < elem.Nodes.rows()-1; j++)
				output << elem.Nodes(j,i) << " ";
			output << elem.Nodes(elem.Nodes.rows()-1,i) << std::endl;
		}
	}

	output << "CELL_TYPES " << Tot_cells_ << std::endl;
	for(auto& elem : Elements){
		for (int i=0; i< elem.nb; i++) {
			output << elem.vtk_type << std::endl;
		}
	}

	output << "CELL_DATA " << Tot_cells_ << std::endl;
	output << "SCALARS marker int 1"<< std::endl;
	output << "LOOKUP_TABLE default"<< std::endl;

	for(auto& elem : Elements){
		for (int i=0; i< elem.nb; i++) {
			for (int j =0; j < elem.Markers.rows()-1; j++)
				output << elem.Markers(j,i) << " ";
			output << elem.Markers(elem.Markers.rows()-1,i) << std::endl;
		}
	}

	// output << "SCALARS DD_weight int 1"<< std::endl;
	// output << "LOOKUP_TABLE default"<< std::endl;

	// for(auto& elem : Elements){
	// 	for (int i=0; i< elem.nb; i++) {
	// 		for (int j =0; j < elem.DD_weight.rows()-1; j++)
	// 			output << elem.DD_weight(j,i) << " ";
	// 		output << elem.DD_weight(elem.DD_weight.rows()-1,i) << std::endl;
	// 	}
	// }

	if (exportRandomFieldOnElement){
		output << "SCALARS Random_Field double 1\n";
		output << "LOOKUP_TABLE default\n";
		for(auto& elem : Elements){
			for (int i=0; i< elem.nb; i++) {
				output << elem.Random_Field[i].average << std::endl;
			}
		}
	}

	if (exportDir){
		output << "SCALARS U double 3\n";
		output << "LOOKUP_TABLE default\n";
		for(auto& elem : Elements){
			for (int i=0; i< elem.nb; i++) {
				output << elem.U(0, i) << " " << elem.U(1, i) << " " << elem.U(2, i) << std::endl;
			}
		}
		output << "SCALARS V double 3\n";
		output << "LOOKUP_TABLE default\n";
		for(auto& elem : Elements){
			for (int i=0; i< elem.nb; i++) {
				output << elem.V(0, i) << " " << elem.V(1, i) << " " << elem.V(2, i) << std::endl;
			}
		}
		output << "SCALARS W double 3\n";
		output << "LOOKUP_TABLE default\n";
		for(auto& elem : Elements){
			for (int i=0; i< elem.nb; i++) {
				output << elem.W(0, i) << " " << elem.W(1, i) << " " << elem.W(2, i) << std::endl;
			}
		}
	}

	if (exportAbaqusFields){
		output << "SCALARS S double 6\n";
		output << "LOOKUP_TABLE default\n";
		for(auto& elem : Elements){
			for (int i=0; i< elem.nb; i++) {
				output << elem.S(0, i) << " " << elem.S(1, i) << " " << elem.S(2, i) << " " <<  elem.S(3, i) << " " << elem.S(4, i) << " " << elem.S(5, i) << std::endl;
			}
		}

		output << "SCALARS E double 6\n";
		output << "LOOKUP_TABLE default\n";
		for(auto& elem : Elements){
			for (int i=0; i< elem.nb; i++) {
				output << elem.E(0, i) << " " << elem.E(1, i) << " " << elem.E(2, i) << " " <<  elem.E(3, i) << " " << elem.E(4, i) << " " << elem.E(5, i) << std::endl;
			}
		}

		// output << "SCALARS SDEG double 1\n";
		// output << "LOOKUP_TABLE default\n";
		// for(auto& elem : Elements){
		// 	for (int i=0; i< elem.nb; i++) {
		// 		output << elem.SDEG(0, i) << std::endl;
		// 	}
		// }
		// output << "SCALARS QUADSCRT double 1\n";
		// output << "LOOKUP_TABLE default\n";
		// for(auto& elem : Elements){
		// 	for (int i=0; i< elem.nb; i++) {
		// 		output << elem.QUADSCRT(0, i) << std::endl;
		// 	}
		// }
	}

	if (exportAbaqusDisplacement){
		output << "POINT_DATA " << nb_vertices_ << std::endl;
		output << "SCALARS U double 3\n";
		output << "LOOKUP_TABLE default\n";
		for (int i=0; i< nb_vertices_; i++) {
			output << Vertices[i].X(0) << " " << Vertices[i].X(1) << " " << Vertices[i].X(2) << std::endl;
		}
	}

	// output << "POINT_DATA " << nb_vertices_ << std::endl;
	// output << "SCALARS extra double 2\n";
	// output << "LOOKUP_TABLE default\n";
	// for (int i=0; i< nb_vertices_; i++) {
	// 	output << Vertices[i].top_surf_indice << " " << Vertices[i].weight << std::endl;
	// }

	output.close();
}

void Mesh::write_point_cloud_vtk(const std::string& filename, int verbosity=0) {

	if(verbosity){
		std::cout << nb_vertices_ << " nodes" << std::endl;
		for(auto& elem : Elements){
			std::cout << elem.nb << " " << elem.type << std::endl;
		}
	}

	std::string outfile(filename+".vtk");
	std::ofstream output;
	output.open(outfile.c_str());
	std::cout << "writing " << outfile << std::endl;

	output << "# vtk DataFile Version 3.0\n";
	output << "Generated by Jean Bénézech\n";
	output << "ASCII\n";
	output << "DATASET UNSTRUCTURED_GRID\n";

	output << "POINTS " << nb_vertices_ << " double" << std::endl;
	for (int i=0; i< nb_vertices_; i++) {
		output << Vertices[i].coord(0) << " " << Vertices[i].coord(1) << " " << Vertices[i].coord(2) << std::endl;
	}

	int second_number = 0;
	for(auto& elem : Elements)
		second_number += (elem.size+1)*elem.nb;


	// if (exportDir){
	// 	output << "SCALARS U double 3\n";
	// 	output << "LOOKUP_TABLE default\n";
	// 	for(auto& elem : Elements){
	// 		for (int i=0; i< elem.nb; i++) {
	// 			output << elem.U(0, i) << " " << elem.U(1, i) << " " << elem.U(2, i) << std::endl;
	// 		}
	// 	}
	// 	output << "SCALARS V double 3\n";
	// 	output << "LOOKUP_TABLE default\n";
	// 	for(auto& elem : Elements){
	// 		for (int i=0; i< elem.nb; i++) {
	// 			output << elem.V(0, i) << " " << elem.V(1, i) << " " << elem.V(2, i) << std::endl;
	// 		}
	// 	}
	// 	output << "SCALARS W double 3\n";
	// 	output << "LOOKUP_TABLE default\n";
	// 	for(auto& elem : Elements){
	// 		for (int i=0; i< elem.nb; i++) {
	// 			output << elem.W(0, i) << " " << elem.W(1, i) << " " << elem.W(2, i) << std::endl;
	// 		}
	// 	}
	// }

	output << "POINT_DATA " << nb_vertices_ << std::endl;
	output << "SCALARS N double 3\n";
	output << "LOOKUP_TABLE default\n";
	for (int i=0; i< nb_vertices_; i++) {
		output << Vertices[i].X(0) << " " << Vertices[i].X(1) << " " << Vertices[i].X(2) << std::endl;
	}

	// output << "POINT_DATA " << nb_vertices_ << std::endl;
	output << "SCALARS thickness double 1\n";
	output << "LOOKUP_TABLE default\n";
	for (int i=0; i< nb_vertices_; i++) {
		output << Vertices[i].thickness << std::endl;
	}

	output.close();
}

void Mesh::write_ori_txt(const std::string& filename) {
	std::ofstream output;
	output.open(filename+"_ori.txt", std::ios::out);
	if (!output.is_open()) {
		std::cout << "Error: Cannot open file " << filename << std::endl;
	}
	std::cout << "writing " << filename+"_ori.txt" << std::endl;

	// En-tete
	output << Tot_cells_ << std::endl;
	// std::cout << "nb_plies_: " << nb_plies_ << std::endl;

	for(auto& elem : Elements){
		for (int i=0; i< elem.nb; i++) {
			output << elem.global_indices[i] << " ";
			
			// output << elem.DD_weight(0, i) << " ";
			
			if (elem.Markers(0, i) <= nb_plies_){
				output << elem.Markers(0, i) << " ";
				output << stacking_sequence[elem.Markers(0, i) - 1] << " ";
			}
			else{
				output << nb_plies_-elem.Markers(0, i) << " 0 "; // "-1 0"
			}
			output << std::setprecision(15) <<  elem.U(0, i) << " " <<  elem.U(1, i) << " " <<  elem.U(2, i) << " ";
			output << std::setprecision(15) << -elem.W(0, i) << " " << -elem.W(1, i) << " " << -elem.W(2, i) << " ";
			output << std::setprecision(15) <<  elem.V(0, i) << " " <<  elem.V(1, i) << " " <<  elem.V(2, i) << std::endl;
			
			// output << "-1 0 1 0 0 -0 -0 1 0 -1 0" << std::endl;
		}
	}
}

void Mesh::write_ori_inp(const std::string& filename) {

	std::ofstream output;
	output.open(filename+"_ori.inp", std::ios::out);
	if (!output.is_open()) {
	std::cout << "Error: Cannot open file " << filename << std::endl;
	}
	std::cout << "writing " << filename+"_ori.inp" << std::endl;


	output << "*ORIENTATION, NAME=ori_glob" << std::endl;
	output << "1., 0., 0., 0., 1., 0." << std::endl;
	output << "1, 0." << std::endl;
	output << "*DISTRIBUTION, NAME=ori_loc_dist, LOCATION=ELEMENT, TABLE=ori_tab" << std::endl;
	output << ", 1.0,  0.0,  0.0,  0.0,  1.0, 0.0" << std::endl;

	for(auto& elem : Elements){
		for (int i=0; i< elem.nb; i++) {
			output << elem.global_indices[i] << ", ";
			if(!isShell){
				output << std::setprecision(15) <<  elem.U(0, i) << ", " <<  elem.U(1, i) << ", " <<  elem.U(2, i) << ", ";
				output << std::setprecision(15) << -elem.W(0, i) << ", " << -elem.W(1, i) << ", " << -elem.W(2, i) << std::endl;
			} else {
				output << std::setprecision(15) <<  elem.V(0, i) << ", " <<  elem.V(1, i) << ", " <<  elem.V(2, i) << ", ";
				output << std::setprecision(15) <<  elem.W(0, i) << ", " <<  elem.W(1, i) << ", " <<  elem.W(2, i) << std::endl;
			}
		}
	}

	output << "*ORIENTATION, NAME=ori_loc" << std::endl;
	output << "ori_loc_dist" << std::endl;
}

void Mesh::extract_AbaqusSets(){

	if (isNsetDone)
		return;

	// ~~~~~~~~~ NSET ~~~~~~~~~
	nb_nset=6;
	nsets.resize(nb_nset);
	masterNodes.resize(nb_nset);

	double xmin, xmax, ymin, ymax, zmin, zmax;
	xmin=ymin=zmin=  10000;
	xmax=ymax=zmax= -10000;
	for (int i=0; i<nb_vertices_;i++) {
		if (xmin > Vertices[i].coord(0)) {xmin = Vertices[i].coord(0);}
		if (xmax < Vertices[i].coord(0)) {xmax = Vertices[i].coord(0);}
		if (ymin > Vertices[i].coord(1)) {ymin = Vertices[i].coord(1);}
		if (ymax < Vertices[i].coord(1)) {ymax = Vertices[i].coord(1);}
		if (zmin > Vertices[i].coord(2)) {zmin = Vertices[i].coord(2);}
		if (zmax < Vertices[i].coord(2)) {zmax = Vertices[i].coord(2);}
	}

	double epsi=0.01;

	for (int i=0; i<nb_vertices_;i++) {
		if (Vertices[i].coord(0)-xmin<epsi) {nsets[0].push_back(i);}
		if (xmax-Vertices[i].coord(0)<epsi) {nsets[1].push_back(i);}
		if (Vertices[i].coord(1)-ymin<epsi) {nsets[2].push_back(i);}
		if (ymax-Vertices[i].coord(1)<epsi) {nsets[3].push_back(i);}
		if (Vertices[i].coord(2)-zmin<epsi) {nsets[4].push_back(i);}
		if (zmax-Vertices[i].coord(2)<epsi) {nsets[5].push_back(i);}
	}

	for(int j=0; j< nb_nset;j++) {
		masterNodes[j](0)=0.0;
		masterNodes[j](1)=0.0;
		masterNodes[j](2)=0.0;
		for(int i=0; i<nsets[j].size(); i++) {
			masterNodes[j](0) += Vertices[nsets[j][i]].coord(0);
			masterNodes[j](1) += Vertices[nsets[j][i]].coord(1);
			masterNodes[j](2) += Vertices[nsets[j][i]].coord(2);
		}
		masterNodes[j](0)/=nsets[j].size();
		masterNodes[j](1)/=nsets[j].size();
		masterNodes[j](2)/=nsets[j].size();
	}

	isNsetDone = true;
	return;
}

void Mesh::elSets_delam(double locX=40., double locZ=160., double radius=10., double radius_offset_X=0.,  double radius_offset_Z=0.){

	if (isNsetDone)
		return;

	// ~~~~~~~~~ ELSET ~~~~~~~~~

	// std::vector<int> new_nset;
	// new_nset.resize(0);
	// for(int nset : nsets)
	// 	new_nsets
	nb_elset+=1;

	for(auto& elem : Elements){
		for(int l=0; l < elem.nb; ++l) {

			Vector3d elem_center = elem.center(Vertices, l);
			// Elliptic region with big axes in a direction (X or Z) where opennig is critical
			double dist= pow(locX-elem_center[0],2)/pow(radius+radius_offset_X,2) +\
						pow(locZ-elem_center[2],2)/pow(radius+radius_offset_Z,2);

			if (dist <= 1 && elem.Markers(0,l)==5){ // 1 for interlayer :: specific to an example
				elem.Markers(0,l)=nb_elset;
			}
		}
	}




	// ~~~~~~~~~ NSET ~~~~~~~~~
	nb_nset=6;
	nsets.resize(nb_nset);
	masterNodes.resize(nb_nset);

	double xmin, xmax, ymin, ymax, zmin, zmax;
	xmin=ymin=zmin=  10000;
	xmax=ymax=zmax= -10000;
	for (int i=0; i<nb_vertices_;i++) {
		if (xmin > Vertices[i].coord(0)) {xmin = Vertices[i].coord(0);}
		if (xmax < Vertices[i].coord(0)) {xmax = Vertices[i].coord(0);}
		if (ymin > Vertices[i].coord(1)) {ymin = Vertices[i].coord(1);}
		if (ymax < Vertices[i].coord(1)) {ymax = Vertices[i].coord(1);}
		if (zmin > Vertices[i].coord(2)) {zmin = Vertices[i].coord(2);}
		if (zmax < Vertices[i].coord(2)) {zmax = Vertices[i].coord(2);}
	}

	double epsi=0.01;

	for (int i=0; i<nb_vertices_;i++) {
		if (Vertices[i].coord(0)-xmin<epsi) {nsets[0].push_back(i);}
		if (xmax-Vertices[i].coord(0)<epsi) {nsets[1].push_back(i);}
		if (Vertices[i].coord(1)-ymin<epsi) {nsets[2].push_back(i);}
		if (ymax-Vertices[i].coord(1)<epsi) {nsets[3].push_back(i);}
		if (Vertices[i].coord(2)-zmin<epsi) {nsets[4].push_back(i);}
		if (zmax-Vertices[i].coord(2)<epsi) {nsets[5].push_back(i);}
	}

	for(int j=0; j< nb_nset;j++) {
		masterNodes[j](0)=0.0;
		masterNodes[j](1)=0.0;
		masterNodes[j](2)=0.0;
		for(int i=0; i<nsets[j].size(); i++) {
			masterNodes[j](0) += Vertices[nsets[j][i]].coord(0);
			masterNodes[j](1) += Vertices[nsets[j][i]].coord(1);
			masterNodes[j](2) += Vertices[nsets[j][i]].coord(2);
		}
		masterNodes[j](0)/=nsets[j].size();
		masterNodes[j](1)/=nsets[j].size();
		masterNodes[j](2)/=nsets[j].size();
	}

	isNsetDone = true;
	return;
}

void Mesh::write_inp(const std::string& filename) {
	std::ofstream output;
	output.open(filename+"_mesh.inp", std::ios::out);
	if (!output.is_open()) {
	std::cout << "Error: Cannot open file " << filename << std::endl;
	} else {
		std::cout << "writing " << filename+"_mesh.inp" << std::endl;
	}

	extract_AbaqusSets();

	std::cout << nb_nset            << " NSETs"  << std::endl;
	std::cout << nb_elset           << " ELSETs"  << std::endl;

	// ~~~~~~~~~ N O D E S ~~~~~~~~~

	output << "*NODE" << std::endl;
	output << nb_vertices_+nb_nset << std::endl; // Plus the number of master nodes
	for(int i=0; i< nb_vertices_; ++i) {
		output << i+1 << ", " << Vertices[i].coord(0) << ", " << Vertices[i].coord(1) << ", " << Vertices[i].coord(2) << std::endl;
	}
	// Adding master nodes
	for(int j=0; j< nb_nset;j++) {
		output << nb_vertices_+1+j << ", " <<  masterNodes[j](0) << ", " << masterNodes[j](1) << ", " << masterNodes[j](2) << std::endl;
	}

	// ~~~~~~~~~ E L E M E N T S ~~~~~~~~~

	output << "******* E L E M E N T S *************" << std::endl;

	std::vector<int> cno = {2, 1, 5, 6, 3, 0, 4, 7};
	for(auto& elem : Elements){
	if(elem.type != "Triangle" || elem.type != "Quadrilateral"){ // We would ony use 3D elements in Abaqus (that could change later)
		output << "*ELEMENT, type=" << elem.abaqus_type << ", ELSET=" << elem.type << std::endl;
		for (int i=0; i< elem.nb; i++) {
			output << elem.global_indices[i] << ", ";
			if(elem.type == "Cohesive"){
				for (int j =0; j < elem.Nodes.rows()-1; j++)
					output << elem.Nodes(cno[j],i)+1 << ", ";
				output << elem.Nodes(cno[elem.Nodes.rows()-1],i)+1 << std::endl;
			}else{
				if(!isShell){
					for (int j =0; j < elem.Nodes.rows()-1; j++)
						output << elem.Nodes(j,i)+1 << ", ";
					output << elem.Nodes(elem.Nodes.rows()-1,i)+1 << std::endl;
				} else {
					for (int j =0; j < elem.Nodes.rows()-1; j++)
						output << elem.Nodes(cno[j],i)+1 << ", ";
					output << elem.Nodes(cno[elem.Nodes.rows()-1],i)+1 << std::endl;
				}
			}
		}
	}
	}

	// for(auto& elem : Elements){
	// if(elem.type != "Triangle" || elem.type != "Quadrilateral"){ // We would ony use 3D elements in Abaqus (that could change later)
	// 	output << "*ELEMENT, type=" << elem.abaqus_type << ", ELSET=" << elem.type << std::endl;
	// 	for (int i=0; i< elem.nb; i++) {
	// 		output << elem.global_indices[i] << ", ";
	// 		for (int j =0; j < elem.Nodes.rows()-1; j++)
	// 			output << elem.Nodes(j,i)+1 << ", ";
	// 		output << elem.Nodes(elem.Nodes.rows()-1,i)+1 << std::endl;
	// 	}
	// }
	// }

	//ELSET
	elsets.resize(nb_elset);
	for(auto& elem : Elements){
		if(elem.type != "Triangle" || elem.type != "Quadrilateral"){ // We would ony use 3D elements in Abaqus (that could change later)
			for (int i=0; i< elem.nb; i++)
				elsets[elem.Markers(0,i)-1].push_back(elem.global_indices[i]);
		}
	}


	// Create an elset for each cohesive layer
	bool is_coh=false;
	for(auto& elem : Elements)
		if(elem.type == "Cohesive")
			is_coh=true;

	if (is_coh){
		Coh_elsets.resize(nb_elset-2); // -1 because {nb_ply-1 = nb_coh} and -1 because {nb_elset = nb_ply +1} to be adjust if some resin elements are in the model
		for(auto& elem : Elements){
			if(elem.type == "Cohesive"){
				int cnt = 0;
				Coh_elsets[cnt].push_back(elem.global_indices[0]);
				for (int i=1; i< elem.nb; i++){
					if(elem.global_indices[i]-elem.global_indices[i-1]>1){
						cnt+=1;
					}
					Coh_elsets[cnt].push_back(elem.global_indices[i]);
				}
			}
		}
	}


	// ~~~~~~~~~ E L S E T S ~~~~~~~~~
	output << "******* E L E M E N T S   S E T S *************" << std::endl;
	output << "*ELSET,ELSET=All_elements" << std::endl;
	for(auto& elem : Elements){
		if(elem.type != "Triangle" || elem.type != "Quadrilateral"){ // We would ony use 3D elements in Abaqus (that could change later)
			for (int i=0; i< elem.nb; i++){
				output << elem.global_indices[i];
				if ((i+1)%10 == 0) { output << "," << std::endl; }
				else { output << ", "; }
			}
		}
	}
	output << std::endl;
	for(int j=0; j< nb_elset;j++) {
		output << "*ELSET,ELSET=elset" << j+1 << std::endl;
		for(int i=0; i< elsets[j].size(); ++i) {
			output << elsets[j][i];
			if (i==elsets[j].size()-1) { output << "," << std::endl; }
			else if ((i+1)%10 == 0) { output << "," << std::endl; }
			else { output << ", "; }
		}
	}

	if (is_coh){
		for(int j=0; j< Coh_elsets.size();j++) {
			output << "*ELSET,ELSET=cohesive_elset" << j+1 << std::endl;
			for(int i=0; i< Coh_elsets[j].size(); ++i) {
				output << Coh_elsets[j][i];
				if (i==Coh_elsets[j].size()-1) { output << "," << std::endl; }
				else if ((i+1)%10 == 0) { output << "," << std::endl; }
				else { output << ", "; }
			}
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

#endif /* end of include guard: MESH_H */
