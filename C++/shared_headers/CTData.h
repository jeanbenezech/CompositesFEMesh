#pragma once

#include "mesh.h"
#include "element.h"
#include "parameters.h"
#include "vector_tools.h"
#include "gridTransformation.h"

#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <Eigen/LU>
#include <Eigen/Geometry>

#define PI 3.14159265

using namespace Eigen;

class Measurement{
	public:
	Vector3d N, Coord;
	double thickness;
	int ID;

	void splitnget(std::string& line, std::string& delim) {
		auto start = 0U;
		auto end = line.find(delim);

		thickness = std::stof(line.substr(start, end - start));
		start = end + 1U;
		end = line.find(delim, start);

		N(0) = std::stof(line.substr(start, end - start));
		start = end + 1U;
		end = line.find(delim, start);
		N(1) = std::stof(line.substr(start, end - start));
		start = end + 1U;
		end = line.find(delim, start);
		N(2) = std::stof(line.substr(start, end - start));
		start = end + 1U;
		end = line.find(delim, start);

		ID = std::stoi(line.substr(start, end - start));
		start = end + 1U;
		end = line.find(delim, start);

		Coord(0) = std::stof(line.substr(start, end - start));
		start = end + 1U;
		end = line.find(delim, start);
		Coord(1) = std::stof(line.substr(start, end - start));
		start = end + 1U;
		end = line.find(delim, start);
		Coord(2) = std::stof(line.substr(start, end - start));
		start = end + 1U;
		end = line.find(delim, start);

	}

	void print(){
		std::cout << ID << " : " << thickness << std::endl;
	}

};

class CT_DATA{
	public:

	// Global gemoetry
	double zmin, zmax, ymin, ymax, xmin, xmax;
	double ymid, zmid;

	// Wrinkles parameters
	std::vector<Measurement> data;

	// General parameters
	int verbosity = 0;

	void readMeasurements(const std::string& filename);

	void init_mesh(Mesh& m, GridTransformation& GT, Parameters& param);

	private:
	std::string delim = ",";
};

void CT_DATA::readMeasurements(const std::string& filename) {

	std::ifstream input;
	input.open(filename, std::ios::in);
	if (!input.is_open()) {
		std::cout << "Error: Cannot open file " << filename << std::endl;
	}

	std::string line {};
	std::getline(input, line);
	std::getline(input, line);
	data.resize(0);

	zmax=0;
	zmin=10000;
	ymin=0;
	ymax=10000;
	xmin=0;
	xmax=10000;

	while(line!=""){
		Measurement tmp;
		tmp.splitnget(line, delim);

		if (tmp.Coord(2)>zmax){ zmax = tmp.Coord(2);}
		if (tmp.Coord(2)<zmin){ zmin = tmp.Coord(2);}
		if (tmp.Coord(1)>ymax){ ymax = tmp.Coord(1);}
		if (tmp.Coord(1)<ymin){ ymin = tmp.Coord(1);}
		if (tmp.Coord(0)>xmax){ xmax = tmp.Coord(0);}
		if (tmp.Coord(0)<xmin){ xmin = tmp.Coord(0);}

		data.push_back(tmp);
		std::getline(input, line);
	}


	ymid = (ymax-ymin)/2.0;
	zmid = (zmax-zmin)/2.0;

}

void CT_DATA::init_mesh(Mesh& m, GridTransformation& GT, Parameters& param){

	int size = data.size();
	double reduction = 10.0;

	int reduced_size = std::floor((size+0.0)/reduction)-1;
	m.Set_vertices(reduced_size);
	m.Vertices.resize(m.Nb_vertices());

	GT.Y.resize(GT.N.size());
	GT.Number_of_associated_point_to_a_node.resize(GT.N.size());
	for (int j=0; j< GT.N.size(); j++) {
		GT.Y[j]=0;
		GT.Number_of_associated_point_to_a_node[j]=0;
	}

	double xmin, xmax, ymin, ymax, zmin, zmax;

	xmin = 50.0;
	xmax = 550.0;
	ymin = 0.0;
	ymax = 150.0;
	zmin = 200;
	zmax = 800;

	for(int i=0; i< m.Nb_vertices(); ++i) { // Loop over all point in the cloud

		int j = reduction*i;
		if (data[j].Coord(1) < ymax && data[j].Coord(2) > zmin && data[j].Coord(2) < zmax){
			if (data[j].Coord(0) > xmin && data[j].Coord(0) < xmax){ // Reduction of the cloud domain


				m.Vertices[i].coord(0)=data[j].Coord(0)-xmin;
				m.Vertices[i].coord(1)=data[j].Coord(1)-ymin; // Translation to the sub cloud domain
				m.Vertices[i].coord(2)=data[j].Coord(2)-zmin;

				m.Vertices[i].thickness = 0.6*(data[j].thickness -  50.0);// / 0.114858; // Minus nominal thickness

				for (int k=0; k< GT.N.size(); k++){
					double dist = sqrt(pow(GT.N[k][0]-m.Vertices[i].coord(0),2) + \
										pow(GT.N[k][1]-m.Vertices[i].coord(2),2));
					if (dist<10){
						// std::cout << "point cloud " << i << " is associated to top surf node " << j << std::endl;
						GT.Number_of_associated_point_to_a_node[k]+=1.0;
						GT.Y(k) += m.Vertices[i].thickness;
					}
				}

				m.Vertices[i].X(0) = data[j].N(0);
				m.Vertices[i].X(1) = data[j].N(1);
				m.Vertices[i].X(2) = data[j].N(2);


			}
		}

	}


	for (int j=0; j< GT.N.size(); j++){
		if(GT.Number_of_associated_point_to_a_node[j]>0)
			GT.Y(j)/=GT.Number_of_associated_point_to_a_node[j];
		// std::cout << "Local thickness of pile " << j << " = " << GT.Y(j) << std::endl;
	}

}

