#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "utils.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#define PI 3.14159265

std::string extract(std::string line) {
	std::string out;
	size_t index = line.find_last_of(":");
	out = line.substr(index+2U, line.size());
	return out;
}

Vector3d vectorExtract(std::string line){
	Vector3d out;
	size_t begin = line.find_last_of(":")+2U;
	std::string myvector = line.substr(begin, line.size());
	size_t coma1 = myvector.find_first_of(",");
	size_t coma2 = myvector.find_last_of(",");
	out(0) = std::stod(myvector.substr(      0U, coma1));
	out(1) = std::stod(myvector.substr(coma1+1U, coma2));
	out(2) = std::stod(myvector.substr(coma2+1U, myvector.size()));
	return out;
}

Vector2d ExtractV2d(std::string line){
	Vector2d out;
	size_t begin = line.find_last_of(":")+2U;
	std::string myvector = line.substr(begin, line.size());
	size_t coma = myvector.find_first_of(",");
	out(0) = std::stod(myvector.substr(      0U, coma));
	out(1) = std::stod(myvector.substr(coma+1U, myvector.size()));
	return out;
}

class Parameters {
public:

	int Shape;
	std::string meshname;
	int add_wrinkles;
	bool add_ramp;
	bool Abaqus_output;
	bool Dune_output;
	int nb_plies;

	int cz_id, resin_id;
	bool isResin, isCZ, isShell, writeTransformedMSH;

	bool GaussianThickness, CornerThickness;
	double sigma, length, ThicknessVar;

	// std::vector<Vector2d> StackSeq;
	std::vector<Vector3d> StackSeq;

	// Wrinkle parameters
	std::string WID;
	std::vector<double> wrinkleSize;
	std::vector<Vector3d> wrinklePos;
	std::vector<Vector3d> wrinkleDamp;
	std::vector<double> wrinkleOri;

	// Ramp parameters
	double rampSize;
	double R;
	Vector3d StartEndinZdir;

	// RVE rotation parameters
	double angleRVE = 0.0;
	double AngleRotateFlangeR = 0.0;
	double AngleRotateFlangeL = 0.0;
	double interior_radius_RVE = 0.0;
	bool rotateRVE = false;
	bool RotateFlanges = false;
	std::string rotateAxis = "X";
	Vector2d rotateStartStop;
	Vector2d centerRot;

	// Geometric parameters
	double X, Y, Z, Height;

	// Path to results
	std::string path_to_abaqus_result;
	std::string AbaqusOdbName;

	void read(const std::string& filename);

	//~Parameters();
// private:
};

// ~~~~~~~~~~~~~~ READ ~~~~~~~~~~~~~~

void Parameters::read(const std::string& filename) {
	std::ifstream input;
	input.open(filename, std::ios::in);
	if (!input.is_open()) {
		std::cout << "Error: Cannot open file" << filename << std::endl;
	} else {
		//std::cout << "Reading " << filename << std::endl;
	}

	std::string line {};
	std::getline(input, line);

	StackSeq.resize(0);
	wrinkleSize.resize(0);
	wrinklePos.resize(0);
	wrinkleDamp.resize(0);
	wrinkleOri.resize(0);

	bool need_default_rotate_start_stop=true;

	while(!line.empty()){

		if (line.find("Shell")!=std::string::npos)
			isShell = std::stoi(extract(line));

		if (line.find("GaussianThickness")!=std::string::npos){
			GaussianThickness = std::stoi(extract(line));
		}
		if (line.find("CornerThickness")!=std::string::npos){
			CornerThickness = std::stoi(extract(line));
		}

		if (line.find("Sigma(d)")!=std::string::npos)
			sigma = std::stod(extract(line));
		if (line.find("Length(d)")!=std::string::npos)
			length = std::stod(extract(line));

		if (line.find("ThicknessVar(d)")!=std::string::npos)
			ThicknessVar = std::stod(extract(line));

		if (line.find("Shape")!=std::string::npos)
			Shape = std::stoi(extract(line));

		if (line.find("Auto_Resin_betw_plies")!=std::string::npos)
			isResin = std::stoi(extract(line));

		if (line.find("cohezive_elements")!=std::string::npos)
			isCZ = std::stoi(extract(line));

		if (line.find("np(i)")!=std::string::npos)
			nb_plies = std::stoi(extract(line));

		if (line.find("name(s)")!=std::string::npos)
			meshname = extract(line)+".msh";

		if (line.find("wrinkle")!=std::string::npos)
			add_wrinkles = std::stoi(extract(line));

		if (line.find("Ramp(b)")!=std::string::npos)
			add_ramp = std::stoi(extract(line));

		if (line.find("R(f)")!=std::string::npos){
			R = std::stod(extract(line));
		}

		if (line.find("X(f)")!=std::string::npos)
			X = std::stod(extract(line));

		if (line.find("Y(f)")!=std::string::npos)
			Y = std::stod(extract(line));
		
		if (line.find("Z(f)")!=std::string::npos)
			Z = std::stod(extract(line));

		if (line.find("Height(f)")!=std::string::npos)
			Height = std::stod(extract(line));

		if (line.find("StartEndinZdir(f)")!=std::string::npos)
			StartEndinZdir = vectorExtract(line);

		if (line.find("Abaqus_output")!=std::string::npos)
			Abaqus_output = std::stoi(extract(line));

		if (line.find("Dune_output")!=std::string::npos)
			Dune_output = std::stoi(extract(line));

		if (line.find("writeTransformedMSH(b)")!=std::string::npos)
			writeTransformedMSH = std::stoi(extract(line));

		if (line.find("Rsize(f)")!=std::string::npos)
			rampSize = std::stod(extract(line));

		if (line.find("WID(s)")!=std::string::npos)
			WID = extract(line);

		if (line.find("Wsize")!=std::string::npos){
			double tmp = std::stof(extract(line));
			wrinkleSize.push_back(tmp);
		}

		if (line.find("Wori")!=std::string::npos){
			double tmp = std::stof(extract(line));
			wrinkleOri.push_back(tmp);
		}

		if (line.find("Wpos")!=std::string::npos){
			Vector3d tmp = vectorExtract(line);
			wrinklePos.push_back(tmp);
		}

		if (line.find("Wdamp")!=std::string::npos){
			Vector3d tmp = vectorExtract(line);
			wrinkleDamp.push_back(tmp);
		}

		if (line.find("Path2result(s)")!=std::string::npos)
			path_to_abaqus_result = extract(line);

		if (line.find("AbaqusOdbName(s)")!=std::string::npos)
			AbaqusOdbName = extract(line);

		if (line.find("f,f,b")!=std::string::npos){
			Vector3d tmp = vectorExtract(line);
			StackSeq.push_back(tmp);
		}
		// if (line.find("f,f,b")!=std::sf,f,btring::npos){
		// 	Vector2d tmp = ExtractV2d(line);
		// 	StackSeq.push_back(tmp);
		// }

		if (line.find("RotateRVE(b)")!=std::string::npos){
			rotateRVE = std::stoi(extract(line));
		}
		if (line.find("AngleRotateRVE(f)")!=std::string::npos){
			angleRVE = std::stof(extract(line));
		}

		if (line.find("RotateAxis(f)")!=std::string::npos)
			rotateAxis = extract(line);

		if (line.find("Rotate_start_end(f)")!=std::string::npos){
			rotateStartStop = ExtractV2d(line);
			need_default_rotate_start_stop=false;
		}


		if (line.find("RotateFlanges(b)")!=std::string::npos){
			RotateFlanges = std::stoi(extract(line));
			// std::cout << "RotateFlanges: " << RotateFlanges << std::endl;
		}
		if (line.find("AngleRotateFlangeRi(f)")!=std::string::npos){
			AngleRotateFlangeR = std::stof(extract(line));
		}
		if (line.find("AngleRotateFlangeLe(f)")!=std::string::npos){
			AngleRotateFlangeL = std::stof(extract(line));
		}

		std::getline(input, line);
	}



	if (need_default_rotate_start_stop){
		// Default value
		if (rotateAxis=="X"){
			rotateStartStop = {0., Z };
		} else if (rotateAxis=="Z"){
			rotateStartStop = {0., X };
		}
	}

	if (rotateRVE){
		if (angleRVE>0){
			interior_radius_RVE = (rotateStartStop(1)-rotateStartStop(0)) / (angleRVE*PI / 180.);
		} else {
			std::cout << "the angle of rotate RVE must be positive" << std::endl; 
			exit( 0 );
		}
	}



	cz_id = -1;
	resin_id = -1;
	if (isResin==1 && isCZ==1){
		cz_id = nb_plies+2;
		resin_id = nb_plies+1;
	} else if (isResin==1 && !isCZ==1) {
		resin_id = nb_plies+1;
	} else if (!isResin==1 && isCZ==1) {
		cz_id = nb_plies+1;
	}

	if (Shape==0) { // Cspar: force this parameter to be false
		rotateRVE= false;
	}

}


#endif /* end of include guard: PARAMETERS_H */
