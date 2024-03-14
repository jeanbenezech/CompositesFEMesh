#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "utils.h"
#include "../mINI/src/mini/ini.h"
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

Vector3d vectorExtract(std::string line, bool isIni=false){
	Vector3d out;
	size_t begin = line.find_last_of(":")+2U;
	if (isIni == true)
		begin = 0U;
	std::string myvector = line.substr(begin, line.size());
	size_t coma1 = myvector.find_first_of(",");
	size_t coma2 = myvector.find_last_of(",");
	out(0) = std::stod(myvector.substr(      0U, coma1));
	out(1) = std::stod(myvector.substr(coma1+1U, coma2));
	out(2) = std::stod(myvector.substr(coma2+1U, myvector.size()));
	return out;
}

Vector2d ExtractV2d(std::string line, bool isIni=false){
	Vector2d out;
	size_t begin = line.find_last_of(":")+2U;
	if (isIni == true)
		begin = 0U;
	std::string myvector = line.substr(begin, line.size());
	size_t coma = myvector.find_first_of(",");
	out(0) = std::stod(myvector.substr(      0U, coma));
	out(1) = std::stod(myvector.substr(coma+1U, myvector.size()));
	return out;
}

std::vector<double> Extractvector(std::string line, int size, bool isIni=false){
	std::vector<double> out(size);
	size_t begin = line.find_last_of(":")+2U;
	if (isIni == true)
		begin = 0U;
	std::string myvector_init = line.substr(begin, line.size());
	std::string myvector = line.substr(begin, line.size());
	auto init_string = 0U;
	for (int i=0; i<size; i++){
		size_t coma = myvector.find_first_of(",");
		auto extracted = myvector.substr(0U, coma);
		out[i] = std::stod(extracted);
		// std::cout << "out["<<i<<"]= " << out[i] << std::endl;
		myvector = myvector_init.substr(coma + 1U, myvector_init.size());
		myvector_init = myvector;
	}
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
	

	// Rigid boundary element
	int RigidBoundary=0;
	int SizeIntervaldZ=2;
	std::vector<double> intervaldZ;

	void read(const std::string& filename);

	//~Parameters();
// private:
};

// ~~~~~~~~~~~~~~ READ ~~~~~~~~~~~~~~

void Parameters::read(const std::string& filename) {
	StackSeq.resize(0);
	wrinkleSize.resize(0);
	wrinklePos.resize(0);
	wrinkleDamp.resize(0);
	wrinkleOri.resize(0);


	// Deprecated
	// std::ifstream input;
	// input.open(filename, std::ios::in);
	// if (!input.is_open()) {
	// 	std::cout << "Error: Cannot open file" << filename << std::endl;
	// } else {
	// 	std::cout << "Reading " << filename << std::endl;
	// }
	// std::string line {};
	// std::getline(input, line);
	// while(!line.empty()){ //// Deprecated

	// 	if (line.find("Shell")!=std::string::npos)
	// 		isShell = std::stoi(extract(line));

	// 	if (line.find("GaussianThickness")!=std::string::npos){
	// 		GaussianThickness = std::stoi(extract(line));
	// 	}
	// 	if (line.find("CornerThickness")!=std::string::npos){
	// 		CornerThickness = std::stoi(extract(line));
	// 	}

	// 	if (line.find("Sigma(d)")!=std::string::npos)
	// 		sigma = std::stod(extract(line));
	// 	if (line.find("Length(d)")!=std::string::npos)
	// 		length = std::stod(extract(line));

	// 	if (line.find("ThicknessVar(d)")!=std::string::npos)
	// 		ThicknessVar = std::stod(extract(line));

	// 	if (line.find("Shape")!=std::string::npos)
	// 		Shape = std::stoi(extract(line));

	// 	if (line.find("Auto_Resin_betw_plies")!=std::string::npos)
	// 		isResin = std::stoi(extract(line));

	// 	if (line.find("cohezive_elements")!=std::string::npos)
	// 		isCZ = std::stoi(extract(line));

	// 	if (line.find("np(i)")!=std::string::npos)
	// 		nb_plies = std::stoi(extract(line));

	// 	if (line.find("name(s)")!=std::string::npos)
	// 		meshname = extract(line)+".msh";

	// 	if (line.find("wrinkle")!=std::string::npos)
	// 		add_wrinkles = std::stoi(extract(line));

	// 	if (line.find("Ramp(b)")!=std::string::npos)
	// 		add_ramp = std::stoi(extract(line));

	// 	if (line.find("R(f)")!=std::string::npos){
	// 		R = std::stod(extract(line));
	// 	}

	// 	if (line.find("X(f)")!=std::string::npos)
	// 		X = std::stod(extract(line));

	// 	if (line.find("Y(f)")!=std::string::npos)
	// 		Y = std::stod(extract(line));
		
	// 	if (line.find("Z(f)")!=std::string::npos)
	// 		Z = std::stod(extract(line));

	// 	if (line.find("Height(f)")!=std::string::npos)
	// 		Height = std::stod(extract(line));

	// 	if (line.find("StartEndinZdir(f)")!=std::string::npos)
	// 		StartEndinZdir = vectorExtract(line);

	// 	if (line.find("Abaqus_output")!=std::string::npos)
	// 		Abaqus_output = std::stoi(extract(line));

	// 	if (line.find("Dune_output")!=std::string::npos)
	// 		Dune_output = std::stoi(extract(line));

	// 	if (line.find("writeTransformedMSH(b)")!=std::string::npos)
	// 		writeTransformedMSH = std::stoi(extract(line));

	// 	if (line.find("Rsize(f)")!=std::string::npos)
	// 		rampSize = std::stod(extract(line));

	// 	if (line.find("WID(s)")!=std::string::npos)
	// 		WID = extract(line);

	// 	if (line.find("Wsize")!=std::string::npos){
	// 		double tmp = std::stof(extract(line));
	// 		wrinkleSize.push_back(tmp);
	// 	}

	// 	if (line.find("Wori")!=std::string::npos){
	// 		double tmp = std::stof(extract(line));
	// 		wrinkleOri.push_back(tmp);
	// 	}

	// 	if (line.find("Wpos")!=std::string::npos){
	// 		Vector3d tmp = vectorExtract(line);
	// 		wrinklePos.push_back(tmp);
	// 	}

	// 	if (line.find("Wdamp")!=std::string::npos){
	// 		Vector3d tmp = vectorExtract(line);
	// 		wrinkleDamp.push_back(tmp);
	// 	}

	// 	if (line.find("Path2result(s)")!=std::string::npos)
	// 		path_to_abaqus_result = extract(line);

	// 	if (line.find("AbaqusOdbName(s)")!=std::string::npos)
	// 		AbaqusOdbName = extract(line);

	// 	if (line.find("f,f,b")!=std::string::npos){
	// 		Vector3d tmp = vectorExtract(line);
	// 		StackSeq.push_back(tmp);
	// 	}
	// 	// if (line.find("f,f,b")!=std::sf,f,btring::npos){
	// 	// 	Vector2d tmp = ExtractV2d(line);
	// 	// 	StackSeq.push_back(tmp);
	// 	// }

	// 	if (line.find("RotateRVE(b)")!=std::string::npos){
	// 		rotateRVE = std::stoi(extract(line));
	// 	}
	// 	if (line.find("AngleRotateRVE(f)")!=std::string::npos){
	// 		angleRVE = std::stof(extract(line));
	// 	}

	// 	if (line.find("RotateAxis(f)")!=std::string::npos)
	// 		rotateAxis = extract(line);

	// 	if (line.find("Rotate_start_end(f)")!=std::string::npos){
	// 		rotateStartStop = ExtractV2d(line);
	// 		need_default_rotate_start_stop=false;
	// 	}


	// 	if (line.find("RotateFlanges(b)")!=std::string::npos){
	// 		RotateFlanges = std::stoi(extract(line));
	// 		// std::cout << "RotateFlanges: " << RotateFlanges << std::endl;
	// 	}
	// 	if (line.find("AngleRotateFlangeRi(f)")!=std::string::npos){
	// 		AngleRotateFlangeR = std::stof(extract(line));
	// 	}
	// 	if (line.find("AngleRotateFlangeLe(f)")!=std::string::npos){
	// 		AngleRotateFlangeL = std::stof(extract(line));
	// 	}

	// 	if (line.find("RigidBoundary(b)")!=std::string::npos){
	// 		RigidBoundary = std::stoi(extract(line));
	// 		// std::cout << "Is RigidBoundary? " << RigidBoundary << std::endl;
	// 	}
	
	// 	if (line.find("dZIntervalsize(i)")!=std::string::npos){
	// 		SizeIntervaldZ = std::stoi(extract(line));
	// 	}

	// 	if (line.find("dZInterval(f)")!=std::string::npos){
	// 		intervaldZ = Extractvector(line,SizeIntervaldZ);
	// 	}


	// 	std::getline(input, line);
	// }

	// INI reading
	mINI::INIFile file(filename+".ini");
	mINI::INIStructure config;
	file.read(config);

	// GENERAL parameters
	meshname = config["GENERAL"].get("name(s)")+".msh";
	if (config["GENERAL"].has("Shell(b)")==true)
		isShell = std::stoi(config["GENERAL"].get("Shell(b)"));
	if (config["GENERAL"].has("Shape(i)")==true)
		Shape = std::stoi(config["GENERAL"].get("Shape(i)"));
	if (config["GENERAL"].has("Auto_Resin_betw_plies(b)")==true)
		isResin = std::stoi(config["GENERAL"].get("Auto_Resin_betw_plies(b)"));
	if (config["GENERAL"].has("cohezive_elements(b)")==true)
		isCZ = std::stoi(config["GENERAL"].get("cohezive_elements(b)"));
	if (config["GENERAL"].has("Abaqus_output(b)")==true)
		Abaqus_output = std::stoi(config["GENERAL"].get("Abaqus_output(b)"));
	if (config["GENERAL"].has("Path2result(s)")==true)
		path_to_abaqus_result = config["GENERAL"].get("Path2result(s)");
	if (config["GENERAL"].has("AbaqusOdbName(s)")==true)
		AbaqusOdbName = config["GENERAL"].get("AbaqusOdbName(s)");
	if (config["GENERAL"].has("writeTransformedMSH(b)")==true)
		writeTransformedMSH = std::stoi(config["GENERAL"].get("writeTransformedMSH(b)"));
	if (config["GENERAL"].has("Dune_output(b)")==true)
		Dune_output = std::stoi(config["GENERAL"].get("Dune_output(b)"));


	// General Geometry
	if (config["GEOMETRY"].has("np(i)")==true)
		nb_plies = std::stoi(config["GEOMETRY"].get("np(i)"));
	if (config["GEOMETRY"].has("X(f)")==true)
		X = std::stoi(config["GEOMETRY"].get("X(f)"));
	if (config["GEOMETRY"].has("Y(f)")==true)
		Y = std::stoi(config["GEOMETRY"].get("Y(f)"));
	if (config["GEOMETRY"].has("Z(f)")==true)
		Z = std::stoi(config["GEOMETRY"].get("Z(f)"));
	if (config["GEOMETRY"].has("Height(f)")==true)
		Height = std::stoi(config["GEOMETRY"].get("Height(f)"));
	if (config["GEOMETRY"].has("R(f)")==true)
		R = std::stoi(config["GEOMETRY"].get("R(f)"));

	// StaqSeq
	StackSeq.resize(nb_plies);
	for(int iter_ply = 0; iter_ply<nb_plies; iter_ply++){
		if (config["STACKSEQ"].has("p"+std::to_string(iter_ply)+"(f,f,b)")==true)
			StackSeq[iter_ply] = vectorExtract(config["STACKSEQ"].get("p"+std::to_string(iter_ply)+"(f,f,b)"), true); // true : from config file, false: from txt older parameters file format
	}

	/// Grid transformation
	if (config["GEO-TRANSFORMATION"].has("Ramp(b)")==true)
		add_ramp = std::stoi(config["GEO-TRANSFORMATION"].get("Ramp(b)"));
	if (config["GEO-TRANSFORMATION"].has("StartEndinZdir(f)")==true)
		StartEndinZdir = vectorExtract(config["GEO-TRANSFORMATION"].get("StartEndinZdir(f)"), true);
	if (config["GEO-TRANSFORMATION"].has("Rsize(f)")==true)
		rampSize = std::stod(config["GEO-TRANSFORMATION"].get("Rsize(f)"));
	if (config["GEO-TRANSFORMATION"].has("RotateRVE(b)")==true)
		rotateRVE = std::stoi(config["GEO-TRANSFORMATION"].get("RotateRVE(b)"));
	if (config["GEO-TRANSFORMATION"].has("AngleRotateRVE(f)")==true)
		angleRVE = std::stof(config["GEO-TRANSFORMATION"].get("AngleRotateRVE(f)"));
	// if (config["GEO-TRANSFORMATION"].has("InteriorRadiusRVE(f)")==true)
	// 	interior_radius_RVE = std::stof(config["GEO-TRANSFORMATION"].get("InteriorRadiusRVE(f)"));
	rotateAxis = config["GEO-TRANSFORMATION"].get("RotateAxis(s)");

	bool need_default_rotate_start_stop=true;
	if (rotateAxis=="X"){
		rotateStartStop = {0., Z };
	} else if (rotateAxis=="Z"){
		rotateStartStop = {0., X };
	}
	if (config["GEO-TRANSFORMATION"].has("Rotate_start_end(f)")==true)
		rotateStartStop = ExtractV2d(config["GEO-TRANSFORMATION"].get("Rotate_start_end(f)"), true);
	if (config["GEO-TRANSFORMATION"].has("RotateFlanges(b)")==true)
		RotateFlanges = std::stoi(config["GEO-TRANSFORMATION"].get("RotateFlanges(b)"));
	if (config["GEO-TRANSFORMATION"].has("AngleRotateFlangeRi(f)")==true)
		AngleRotateFlangeR = std::stof(config["GEO-TRANSFORMATION"].get("AngleRotateFlangeRi(f)"));
	if (config["GEO-TRANSFORMATION"].has("AngleRotateFlangeLe(f)")==true)
		AngleRotateFlangeL = std::stof(config["GEO-TRANSFORMATION"].get("AngleRotateFlangeLe(f)"));

	if (config["GEO-TRANSFORMATION"].has("RigidBoundary(b)")==true)
		RigidBoundary = std::stoi(config["GEO-TRANSFORMATION"].get("RigidBoundary(b)"));
	if (config["GEO-TRANSFORMATION"].has("dZIntervalsize(i)")==true)
		SizeIntervaldZ = std::stoi(config["GEO-TRANSFORMATION"].get("dZIntervalsize(i)"));
	if (RigidBoundary==true)
		if (config["GEO-TRANSFORMATION"].has("dZInterval(f)")==true)
			intervaldZ = Extractvector(config["GEO-TRANSFORMATION"].get("dZInterval(f)"), SizeIntervaldZ, true);


	// GAUSSIAN parameters
	if (config["GENERAL"].has("GaussianThickness(b)")==true)
		GaussianThickness = std::stoi(config["GENERAL"].get("GaussianThickness(b)"));
	if (config["GAUSSIAN"].has("Sigma(d)")==true)
		sigma = std::stod(config["GAUSSIAN"].get("Sigma(d)"));
	if (config["GAUSSIAN"].has("Length(d)")==true)
		length = std::stod(config["GAUSSIAN"].get("Length(d)"));

	// std::cout << "length: " << length << std::endl;

	// CORNER THICKNESS parameters
	if (config["GENERAL"].has("CornerThickness(b)")==true)
		CornerThickness = std::stoi(config["GENERAL"].get("CornerThickness(b)"));
	if (config["GEO-TRANSFORMATION"].has("ThicknessVar(d)")==true)
		ThicknessVar = std::stod(config["GEO-TRANSFORMATION"].get("ThicknessVar(d)"));
	// std::cout << "CornerThickness: " << CornerThickness << std::endl;
	// std::cout << "ThicknessVar: " << ThicknessVar << std::endl;
	// WRINKLES
	if (config["WRINKLES"].has("wrinkle(i)")==true)
		add_wrinkles = std::stoi(config["WRINKLES"].get("wrinkle(i)"));
	WID = config["WRINKLES"].get("WID(s)");

	wrinkleSize.resize(add_wrinkles);
	wrinkleOri.resize(add_wrinkles);
	wrinklePos.resize(add_wrinkles);
	wrinkleDamp.resize(add_wrinkles);
	for(int iter_w = 0; iter_w<add_wrinkles; iter_w++){
		if (config["WRINKLES"].has("Wsize"+std::to_string(iter_w)+"(f)")==true)
			wrinkleSize[iter_w] = std::stof(config["WRINKLES"].get("Wsize"+std::to_string(iter_w)+"(f)"));
		if (config["WRINKLES"].has("Wori"+std::to_string(iter_w)+"(f)")==true)
			wrinkleOri[iter_w] = std::stof(config["WRINKLES"].get("Wori"+std::to_string(iter_w)+"(f)"));
		if (config["WRINKLES"].has("Wpos"+std::to_string(iter_w)+"(f)")==true)
			wrinklePos[iter_w] = vectorExtract(config["WRINKLES"].get("Wpos"+std::to_string(iter_w)+"(f)"), true);
		if (config["WRINKLES"].has("Wdamp"+std::to_string(iter_w)+"(f)")==true)
			wrinkleDamp[iter_w] = vectorExtract(config["WRINKLES"].get("Wdamp"+std::to_string(iter_w)+"(f)"), true);
	}

	/// Ajustment on parameters choices to avoid incompatibility
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

	if (Shape==0 || Shape==2) { // Cspar: force this parameter to be false
		rotateRVE= false;
	}

}


#endif /* end of include guard: PARAMETERS_H */
