#pragma once

#include "mesh.h"
#include "element.h"
#include "parameters.h"
#include "vector_tools.h"

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

	private:
	std::string delim = ",";
	// double epsilon = 0.1;
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


		// double theta = -90.0;
		// double c = cos(theta*PI/180.0);
		// double s = sin(theta*PI/180.0);
		// rotZ = {{c,s,0},{-s,c,0},{0,0,1}};
		// rotY = {{c,0,-s},{0,1,0},{s,0,c}};


		if (data[i].Coord(2)>zmax){ zmax = data[i].Coord(2);}
		if (data[i].Coord(2)<zmin){ zmin = data[i].Coord(2);}
		if (data[i].Coord(1)>ymax){ ymax = data[i].Coord(1);}
		if (data[i].Coord(1)<ymin){ ymin = data[i].Coord(1);}
		if (data[i].Coord(0)>xmax){ xmax = data[i].Coord(0);}
		if (data[i].Coord(0)<xmin){ xmin = data[i].Coord(0);}

		data.push_back(tmp);
		std::getline(input, line);
	}

	ymid = (ymax-ymin)/2.0;
	zmid = (zmax-zmin)/2.0;

	// for(int node=0; node < m.Nb_vertices(); ++node) {

	// }

	// if (verbosity>0){
	// 	std::cout << "Xmin : " << xmin << std::endl;
	// 	std::cout << "Xmax : " << xmax << std::endl;
	// 	std::cout << "Ymin : " << ymin << std::endl;
	// 	std::cout << "Ymax : " << ymax << std::endl;
	// 	std::cout << "Zmin : " << zmin << std::endl;
	// 	std::cout << "Zmax : " << zmax << std::endl;
	// }



	// std::cout << "Ymid : " << ymid << std::endl;
}

