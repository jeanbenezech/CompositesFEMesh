#ifndef MESH_H
#define MESH_H

#include "utils.h"
#include "element.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Dense>

using namespace Eigen;

class Mesh {
public:
	// PRINCIPAL VARIABLES
	Matrix<float, Dynamic, Dynamic> vertices;
	std::vector<Element> Elements;

	std::string mesh_type;

	int nb_nset=0;
	int nb_elset=0;

	// Matrix<float, Dynamic, Dynamic> sortie;
	std::vector<float> stacking_sequence;
	int nb_corner, type;
	std::vector<Vector2f> P;
	float theta;

	// FUNCTION
	void read_mesh(const std::string& filename);
	void read_msh(const std::string& filename, bool only_3D, int cz_id);
	void read_points(const std::string& filename);
	void write_mesh(const std::string& filename);
	void write_msh(const std::string& filename, int verbosity);
	bool exportDir=false;
	void write_vtk(const std::string& filename, int verbosity);
	void write_inp(const std::string& filename);
	void write_ori_txt(const std::string& filename);
	void write_ori_inp(const std::string& filename);

	void initialise(const int nb_plies);
	void print(std::string type, int number);

	int Nb_vertices(){return nb_vertices_;}
	int Nb_Tot_cell(){return Tot_cells_;}

	int Nb_plies(){return nb_plies_;}

	//~Mesh();
private:
	int nb_vertices_;
	int Tot_cells_;

	int nb_plies_;

};

// ~~~~~~~~~~~~~~ INIT ~~~~~~~~~~~~~~

void Mesh::initialise(const int nb_plies) {

	stacking_sequence.resize(nb_plies);
	
	stacking_sequence[0] =  45.0;
	stacking_sequence[1] = -45.0;
	stacking_sequence[2] =  90.0;
	stacking_sequence[3] =   0.0;
	stacking_sequence[4] = -45.0;
	stacking_sequence[5] =  45.0;
	stacking_sequence[6] =  45.0;
	stacking_sequence[7] = -45.0;
	stacking_sequence[8] =  90.0;
	stacking_sequence[9] =   0.0;
	stacking_sequence[10]= -45.0;
	stacking_sequence[11]=  45.0;
	for(int i=0; i< 12;i++){
		stacking_sequence[12+i] = stacking_sequence[11-i];
	}

}

// ~~~~~~~~~~~~~~ PRINT ~~~~~~~~~~~~~~

void Mesh::print(std::string type, int number) {
	for(auto& elem : Elements){
		if(elem.type != "Triangle" || elem.type != "Quadrilateral"){
			for(int n = 0; n < elem.Nodes.rows(); n++) {
				std::cout << n << " : ";
				std::cout << vertices(0, elem.Nodes(n, number)) << ", ";
				std::cout << vertices(1, elem.Nodes(n, number)) << ", ";
				std::cout << vertices(2, elem.Nodes(n, number)) << std::endl;
			}
		}
	}
}

// ~~~~~~~~~~~~~~ READ ~~~~~~~~~~~~~~

void Mesh::read_mesh(const std::string& filename) {
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
		vertices(0,d-1)=x;
		vertices(1,d-1)=y;
		vertices(2,d-1)=z;
		//~ if (i<10){std::cout << vertices(0,i) << ", " << vertices(1,i) << ", " << vertices(2,i) << std::endl;}
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
		Elements[incr_elem].initialise("Tetrahedron", ltet.size(), ltet[0][2]-1);
		Elements[incr_elem].fill(ltet, nb_elset);
		incr_elem += 1;
	}
	if(ishex){
		Elements[incr_elem].initialise("Hexahedron", lhex.size(), lhex[0][2]-1);
		Elements[incr_elem].fill(lhex, nb_elset);
		incr_elem += 1;
	}
	if(iscoh){
		Elements[incr_elem].initialise("Cohesive", lcoh.size(), lcoh[0][2]-1);
		Elements[incr_elem].fill(lcoh, nb_elset);
		incr_elem += 1;
	}
	if(iswed){
		Elements[incr_elem].initialise("Prism", lwed.size(), lwed[0][2]-1);
		Elements[incr_elem].fill(lwed, nb_elset);
		incr_elem += 1;
	}
	if(ispyr){
		Elements[incr_elem].initialise("Pyramid", lpyr.size(), lpyr[0][2]-1);
		Elements[incr_elem].fill(lpyr, nb_elset);
		incr_elem += 1;
	}

}

void Mesh::read_points(const std::string& filename) {
	std::ifstream input;
	input.open(filename, std::ios::in);
	if (!input.is_open()) {
		std::cout << "Error: Cannot open file" << std::endl;
	} else{
		std::cout << "Reading " << filename << std::endl;
	}

	input >> nb_corner;

	// if (type<4) theta = 90;

	P.resize(nb_corner);
	for (int i=0; i<P.size();i++) {
		input >> P[i](0) >> P[i](1);
	}
}

// ~~~~~~~~~~~~~~ DUMP ~~~~~~~~~~~~~~

void Mesh::write_mesh(const std::string& filename) {
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
	for(int i=0; i< nb_vertices_; i++) {
		output << vertices(0,i) << " " << vertices(1,i) << " " << vertices(2,i) << " " << i << std::endl;
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
	std::cout << "Error: Cannot open file" << std::endl;
	}
	std::cout << "Writting " << filename+".msh" << std::endl;

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
		output << i+1 << " " << vertices(0,i) << " " << vertices(1,i) << " " << vertices(2,i) << std::endl;
	}
	output << "$EndNodes" << std::endl;
	output << "$Elements" << std::endl;
	output << Tot_cells_ << std::endl;

	for(auto& elem : Elements){
		for (int i=0; i< elem.nb; i++) {
			output << elem.global_indices(i) << " " << elem.msh_type << " ";
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
	std::cout << "Writting " << outfile << std::endl;

	output << "# vtk DataFile Version 3.0\n";
	output << "Generated by Jean Bénézech\n";
	output << "ASCII\n";
	output << "DATASET UNSTRUCTURED_GRID\n";

	output << "POINTS " << nb_vertices_ << " float" << std::endl;
	for (int i=0; i< nb_vertices_; i++) {
		output << vertices(0,i) << " " << vertices(1,i) << " " << vertices(2,i) << std::endl;
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

	if (exportDir){
		output << "SCALARS U float 3\n";
		output << "LOOKUP_TABLE default\n";
		for(auto& elem : Elements){
			for (int i=0; i< elem.nb; i++) {
				output << elem.U(0, i) << " " << elem.U(1, i) << " " << elem.U(2, i) << std::endl;
			}
		}
		output << "SCALARS V float 3\n";
		output << "LOOKUP_TABLE default\n";
		for(auto& elem : Elements){
			for (int i=0; i< elem.nb; i++) {
				output << elem.V(0, i) << " " << elem.V(1, i) << " " << elem.V(2, i) << std::endl;
			}
		}
		output << "SCALARS W float 3\n";
		output << "LOOKUP_TABLE default\n";
		for(auto& elem : Elements){
			for (int i=0; i< elem.nb; i++) {
				output << elem.W(0, i) << " " << elem.W(1, i) << " " << elem.W(2, i) << std::endl;
			}
		}
	}

	output.close();
}

void Mesh::write_ori_txt(const std::string& filename) {
	std::ofstream output;
	output.open(filename+"_ori.txt", std::ios::out);
	if (!output.is_open()) {
	std::cout << "Error: Cannot open file" << std::endl;
	}
	std::cout << "Writting " << filename+"_ori.txt" << std::endl;

	// En-tete
	output << Tot_cells_ << std::endl;

	for(auto& elem : Elements){
		for (int i=0; i< elem.nb; i++) {
			output << elem.global_indices(i) << " ";
			output << elem.Markers(0, i) << " ";
			output << stacking_sequence[elem.Markers(0, i) - 1] << " ";
			output << elem.U(0, i) << " " << elem.U(1, i) << " " << elem.U(2, i) << " ";
			output << elem.V(0, i) << " " << elem.V(1, i) << " " << elem.V(2, i) << " ";
			output << elem.W(0, i) << " " << elem.W(1, i) << " " << elem.W(2, i) << std::endl;
		}
	}
}

void Mesh::write_ori_inp(const std::string& filename) {

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

	for(auto& elem : Elements){
		for (int i=0; i< elem.nb; i++) {
			output << elem.global_indices(i) << ", ";
			output << elem.U(0, i) << ", " << elem.U(1, i) << ", " << elem.U(2, i) << ", ";
			output << elem.V(0, i) << ", " << elem.V(1, i) << ", " << elem.V(2, i) << std::endl;
		}
	}

	output << "*ORIENTATION, NAME=ori_loc" << std::endl;
	output << "ori_loc_dist" << std::endl;
}

void Mesh::write_inp(const std::string& filename) {
	std::ofstream output;
	output.open(filename+"_mesh.inp", std::ios::out);
	if (!output.is_open()) {
	std::cout << "Error: Cannot open file" << std::endl;
	} else {
		std::cout << "Writting " << filename+"_mesh.inp" << std::endl;
	}

	// ~~~~~~~~~ NSET and ElSET preparation ~~~~~~~~~
	std::vector<std::vector<int>> nsets, elsets;
	std::vector<Vector3f> masterNodes;

	// ~~~~~~~~~ NSET ~~~~~~~~~
	nb_nset=6;
	nsets.resize(nb_nset);
	masterNodes.resize(nb_nset);

	float xmin, xmax, ymin, ymax, zmin, zmax;
	xmin=ymin=zmin=  10000;
	xmax=ymax=zmax= -10000;
	for (int i=0; i<nb_vertices_;i++) {
		if (xmin > vertices(0, i)) {xmin = vertices(0, i);}
		if (xmax < vertices(0, i)) {xmax = vertices(0, i);}
		if (ymin > vertices(1, i)) {ymin = vertices(1, i);}
		if (ymax < vertices(1, i)) {ymax = vertices(1, i);}
		if (zmin > vertices(2, i)) {zmin = vertices(2, i);}
		if (zmax < vertices(2, i)) {zmax = vertices(2, i);}
	}

	float epsi=0.01;

	for (int i=0; i<nb_vertices_;i++) {
		if (vertices(0, i)-xmin<epsi) {nsets[0].push_back(i);}
		if (xmax-vertices(0, i)<epsi) {nsets[1].push_back(i);}
		if (vertices(1, i)-ymin<epsi) {nsets[2].push_back(i);}
		if (ymax-vertices(1, i)<epsi) {nsets[3].push_back(i);}
		if (vertices(2, i)-zmin<epsi) {nsets[4].push_back(i);}
		if (zmax-vertices(2, i)<epsi) {nsets[5].push_back(i);}
	}

	for(int j=0; j< nb_nset;j++) {
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
	}


	std::cout << nb_nset            << " NSETs"  << std::endl;
	std::cout << nb_elset           << " ELSETs"  << std::endl;

	// ~~~~~~~~~ N O D E S ~~~~~~~~~

	output << "*NODE" << std::endl;
	output << nb_vertices_+nb_nset << std::endl; // Plus the number of master nodes
	for(int i=0; i< nb_vertices_; ++i) {
		output << i+1 << ", " << vertices(0,i) << ", " << vertices(1,i) << ", " << vertices(2,i) << std::endl;
	}
	// Adding master nodes
	for(int j=0; j< nb_nset;j++) {
		output << nb_vertices_+1+j << ", " <<  masterNodes[j](0) << ", " << masterNodes[j](1) << ", " << masterNodes[j](2) << std::endl;
	}

	// ~~~~~~~~~ E L E M E N T S ~~~~~~~~~

	output << "******* E L E M E N T S *************" << std::endl;

	for(auto& elem : Elements){
	if(elem.type != "Triangle" || elem.type != "Quadrilateral"){ // We would ony use 3D elements in Abaqus (that could change later)
		output << "*ELEMENT, type=" << elem.abaqus_type << ", ELSET=" << elem.type << std::endl;
		for (int i=0; i< elem.nb; i++) {
			output << elem.global_indices(i) << ", ";
			for (int j =0; j < elem.Nodes.rows()-1; j++)
				output << elem.Nodes(j,i)+1 << ", ";
			output << elem.Nodes(elem.Nodes.rows()-1,i)+1 << std::endl;
		}
	}
	}

	//ELSET
	elsets.resize(nb_elset);
	for(auto& elem : Elements){
		if(elem.type != "Triangle" || elem.type != "Quadrilateral"){ // We would ony use 3D elements in Abaqus (that could change later)
			for (int i=0; i< elem.nb; i++)
				elsets[elem.Markers(0,i)-1].push_back(elem.global_indices(i));
		}
	}

	// ~~~~~~~~~ E L S E T S ~~~~~~~~~
	output << "******* E L E M E N T S   S E T S *************" << std::endl;
	output << "*ELSET,ELSET=All_elements" << std::endl;
	for(auto& elem : Elements){
		if(elem.type != "Triangle" || elem.type != "Quadrilateral"){ // We would ony use 3D elements in Abaqus (that could change later)
			for (int i=0; i< elem.nb; i++){
				output << elem.global_indices(i);
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
