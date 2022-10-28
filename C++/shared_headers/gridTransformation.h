#pragma once

#include "mesh.h"
#include "element.h"
#include "parameters.h"
#include "vector_tools.h"
#include "spline.h"

#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <Eigen/LU>
#include <Eigen/Geometry>

#include <omp.h>

#include <random>
#include <map> // for visualization

#include <algorithm>

#define PI 3.14159265

using namespace Eigen;

template<typename V>
double distance(V& p1, V& p2){
	double dist=0.0;
	for (int i = 0; i<p1.size(); i++)
		dist += pow(p1(i)-p2(i), 2);
	dist = sqrt(dist);
	return dist;
}

template<typename T>
std::vector<T> arange(T start, T stop, T step = 1) {
    std::vector<T> values;
    for (T value = start; value < stop; value += step)
        values.push_back(value);
    return values;
}

class GridTransformation{
	public:

	// Global gemoetry
	double z1, z2, z3, z4;
	double delta_max;
	double zmin, zmax, ymin, ymax, xmin, xmax;
	double ymid;

	double local_xmax;
	double local_ymin;
	double local_ymax;

	double threshold_avoiding_deformation_of_corners;

	// Wrinkles parameters
	Vector3d defect_location, damping, center_rot;
	double wrinkleOri;
	double ref_lenght;
	double defect_size;
	double interior_radius;


	// Gaussian random field
	MatrixXd K, L;
	std::vector<Vector2d> N;
	std::vector<Vector3d> Flat_map;
	std::vector<double> distFBS;
	VectorXd Z, Y;
	std::vector<int> indices_map_top_surf;
	std::vector<double> weight;
	std::vector<int> associated_top_surf_node;

	// Wrinkles from CT
	std::vector<int> Number_of_associated_point_to_a_node;

	// General parameters
	int verbosity = 0;

	void initialise(Mesh& m, Parameters& param);
	std::tuple<double,double,double> ramp(Vector3d& point);
	void Csection_wrinkles(Vector3d& point, Vector3d& normal, int number, Parameters& param, std::tuple<double,double,double>& ramp_param);
	void Csection_wrinkles_2(Vector3d& point, Vector3d& normal, int number, Parameters& param, std::tuple<double,double,double>& ramp_param);
	void Csection_wrinkles_splined(Vector3d& point, Vector3d& normal, int number, Parameters& param, std::tuple<double,double,double>& ramp_param);

	void flat_transformation(int& node, Vector3d& point, Parameters& param, std::tuple<double,double,double>& ramp_param);
	void prepare_GT_for_Flat_model(int& node, Vector3d& point, Parameters& param, std::tuple<double,double,double>& ramp_param);
	void AssociatedtoTopSurfaceNode(Mesh& m, Parameters& param);
	void Gaussian_random_field_K_initialisation(Parameters& param);
	void Transfo_CT(Mesh& CTmesh);
	void Apply_Gaussian_random_field(Vector3d& point, Vector3d& normal, double& value);

	void Corner_thickness(Vector3d& point, Vector3d& normal, Parameters& param, std::tuple<double,double,double>& ramp_param, double CTV);

	void wrinkles(Vector3d& point);

	private:
	// double epsilon = 0.1;
};

void GridTransformation::initialise(Mesh& m, Parameters& param) {
	zmax=0;
	zmin=m.Vertices[0].coord(2);
	ymin=0;
	ymax=m.Vertices[0].coord(1);
	xmin=0;
	xmax=m.Vertices[0].coord(0);

	for(int node=0; node < m.Nb_vertices(); ++node) {
		if (m.Vertices[node].coord(2)>zmax){ zmax = m.Vertices[node].coord(2);}
		if (m.Vertices[node].coord(2)<zmin){ zmin = m.Vertices[node].coord(2);}
		if (m.Vertices[node].coord(1)>ymax){ ymax = m.Vertices[node].coord(1);}
		if (m.Vertices[node].coord(1)<ymin){ ymin = m.Vertices[node].coord(1);}
		if (m.Vertices[node].coord(0)>xmax){ xmax = m.Vertices[node].coord(0);}
		if (m.Vertices[node].coord(0)<xmin){ xmin = m.Vertices[node].coord(0);}
	}

	if (verbosity>0){
		std::cout << "Xmin : " << xmin << std::endl;
		std::cout << "Xmax : " << xmax << std::endl;
		std::cout << "Ymin : " << ymin << std::endl;
		std::cout << "Ymax : " << ymax << std::endl;
		std::cout << "Zmin : " << zmin << std::endl;
		std::cout << "Zmax : " << zmax << std::endl;
	}

	// double deltaz = (zmax-zmin) / 16.0;
	// z1 = zmin + 2.0*deltaz;
	// z2 = zmin + 7.0*deltaz;
	// z3 = zmin + 9.0*deltaz;
	// z4 = zmin + 14.0*deltaz;
	// double deltaz = (zmax-zmin) / 16.0;
	z1 = zmin + param.StartEndinZdir(0);
	z2 = z1 + param.StartEndinZdir(1);
	z3 = z2 + param.StartEndinZdir(2);
	z4 = z3 + param.StartEndinZdir(1);

	local_xmax = param.Height;
	local_ymin = param.R;
	local_ymax = param.X+local_ymin;

	delta_max = param.rampSize;
	threshold_avoiding_deformation_of_corners = param.R+param.Y; // Exterior radius

	interior_radius = param.R;

	ymid = (ymax-ymin)/2.0;


	size_t lastindex = param.meshname.find_last_of(".");
	std::string mesh_name = param.meshname.substr(0, lastindex);
	/* Write C-spar parameters to be able to use them in DUNE */
	std::string filename= "DUNE/"+mesh_name+"_ParamForWrinklesRepresentationInCspar.txt";
	std::ofstream output;
	output.open(filename, std::ios::out);
	if (!output.is_open()) {
		std::cout << "Error: Cannot open file " << filename << std::endl;
	}
	output << "Global boundaries of the C-spar (xmin, xmax, ymin, ymax, zmin, zmax):" << std::endl;
	output << xmin << std::endl;
	output << xmax << std::endl;
	output << ymin << std::endl;
	output << ymax << std::endl;
	output << zmin << std::endl;
	output << zmax << std::endl;

	output << std::endl;
	output << "Ramp parameters of the C-spar (delta_max, z1, z2, z3, z4):" << std::endl;
	output << delta_max << std::endl;
	output << z1 << std::endl;
	output << z2 << std::endl;
	output << z3 << std::endl;
	output << z4 << std::endl;

	output << std::endl;
	output << "Wrinkles creation specific parameters of the C-spar (localXmax, localYmin, localYmax, threshold_avoiding_deformation_of_corners, interior_radius, big_radius, Cspar_width):" << std::endl;
	output << local_xmax << std::endl;
	output << local_ymin << std::endl;
	output << local_ymax << std::endl;
	output << threshold_avoiding_deformation_of_corners << std::endl;
	output << interior_radius << std::endl;
	output << param.R << std::endl;
	output << param.X << std::endl;

	// std::cout << "Ymid : " << ymid << std::endl;
}

void GridTransformation::flat_transformation(int& node, Vector3d& point, Parameters& param, std::tuple<double,double,double>& ramp_param){
	double local_ymid = (local_ymin+local_ymax)/2.0;

	Vector3d ramp;
	ramp[2] = 0.0;
	ramp[0] = -std::get<0>(ramp_param);
	if(point[1]>=local_ymid)
		ramp[1] = -std::get<1>(ramp_param);
	else
		ramp[1] = std::get<2>(ramp_param);

	Vector3d init = point-ramp;

	Vector3d ref, moved;
	double dist_from_bottom_surf=0.0;
	moved = init;
	if (init[1]>=local_ymax){
		if(init[0]<=local_xmax){
			dist_from_bottom_surf = (init[1] - local_ymax) - interior_radius;
		} else {
			dist_from_bottom_surf = sqrt((init[1]-local_ymax)*(init[1]-local_ymax)+(init[0]-local_xmax)*(init[0]-local_xmax))-interior_radius;
		}
		ref[0] = local_xmax + interior_radius + dist_from_bottom_surf;
		ref[1] = param.X + interior_radius;
		moved[0] = ref[0];
		double local_radius = interior_radius + dist_from_bottom_surf;

		if(init[0]<=local_xmax){
			double theta = PI/2.0;
			moved[1] = ref[1] + param.R*theta + abs(init[0]-ref[0])-local_radius;
		} else {
			double d = sqrt((init[1]-ref[1])*(init[1]-ref[1])+(init[0]-ref[0])*(init[0]-ref[0]));
			double theta = 2*asin(d/2.0/local_radius);

			moved[1] = ref[1] + theta*param.R;
		}
	}
	else if (init[1]<=local_ymin){
		if(init[0]<=local_xmax){
			dist_from_bottom_surf = - init[1] + local_ymin - interior_radius;
		} else {
			dist_from_bottom_surf = sqrt((init[1]-local_ymin)*(init[1]-local_ymin)+(init[0]-local_xmax)*(init[0]-local_xmax))-interior_radius;
		}
		ref[0] = local_xmax + interior_radius + dist_from_bottom_surf;
		ref[1] = interior_radius;
		moved[0] = ref[0];
		double local_radius = interior_radius + dist_from_bottom_surf;

		if(init[0]<=local_xmax){
			double theta = PI/2.0;
			moved[1] = ref[1] - (param.R*theta + abs(init[0]-ref[0])-local_radius);
		} else {
			double d = sqrt((init[1]-ref[1])*(init[1]-ref[1])+(init[0]-ref[0])*(init[0]-ref[0]));
			double theta = 2*asin(d/2.0/local_radius);

			moved[1] = ref[1] - theta*param.R;
		}
	} else {
		dist_from_bottom_surf = init[0] - local_xmax - interior_radius;
	}

	Flat_map[node] = moved;
	distFBS[node] = dist_from_bottom_surf;
	weight[node] = dist_from_bottom_surf/param.Y;

	if (dist_from_bottom_surf>(param.Y-0.001)){
		Vector2d tmp;
		tmp[0]=moved[1];
		tmp[1]=moved[2];
		N.push_back(tmp);
		indices_map_top_surf.push_back(node);
	}

}

void GridTransformation::prepare_GT_for_Flat_model(int& node, Vector3d& point, Parameters& param, std::tuple<double,double,double>& ramp_param){

	local_ymin = 0.0;
	double local_ymid = (local_ymin+local_ymax)/2.0;

	Vector3d init = point;
	
	double dist_from_bottom_surf = param.Y - abs(init[1]);
	if(abs(dist_from_bottom_surf)<0.0000001)
		dist_from_bottom_surf=0.0;
	// std::cout << dist_from_bottom_surf << std::endl;

	Flat_map[node] = init; // No transformation, same coordinate
	distFBS[node] = dist_from_bottom_surf;

	// weight[node] = dist_from_bottom_surf/param.Y; // Linear variation through thickness
	double barrier = param.Y*1.0/2.0;
	if(dist_from_bottom_surf<=barrier){
		// weight[node] = pow(dist_from_bottom_surf / barrier, 1/2);
		weight[node] = sqrt(dist_from_bottom_surf / barrier);
	} else if(dist_from_bottom_surf<param.Y) {
		// weight[node] = pow((param.Y-dist_from_bottom_surf) / barrier, 1/2);
		weight[node] = sqrt((param.Y-dist_from_bottom_surf) / barrier);
	} else {
		weight[node] = 0.0;
	}


	std::cout << "param.Y: " << param.Y << std::endl;
	std::cout << "dist_from_bottom_surf: " << dist_from_bottom_surf << std::endl;
	std::cout << "weight[node]: " << weight[node] << std::endl;

	if (dist_from_bottom_surf>(param.Y-0.001)){
		Vector2d tmp;
		tmp[0]=init[0];
		tmp[1]=init[2];
		N.push_back(tmp);
		indices_map_top_surf.push_back(node);
	}

}

void GridTransformation::AssociatedtoTopSurfaceNode(Mesh& m, Parameters& param){

	// omp_set_dynamic(0);
	// omp_set_num_threads(10);

	int nb_vertice_top_plan = m.Nb_vertices()/(m.Nb_plies()+1);
	Vector3d ext;
	if (param.Shape == 0) {
		ext = {1.0, 0.0, 0.0}; // For flat domain of the C-spar
	} else if (param.Shape == 1) {
		ext = {0.0, 1.0, 0.0}; // For flat beam
	}

	associated_top_surf_node.resize(m.Nb_vertices());

	// std::cout << " size : " << nb_vertice_top_plan << std::endl;
	// std::cout << " Flat_map size : " << Flat_map.size() << std::endl;
	// std::cout << " distFBS size : " << distFBS.size() << std::endl;

	// #pragma omp parallel for schedule(guided)
	for(int i=0; i<m.Nb_vertices(); i++){
		Vector3d proj2topSurf = Flat_map[i] + (param.Y-distFBS[i])*ext;
		// std::cout << "i : " << i << std::endl;
		for(int j=0; j<indices_map_top_surf.size(); j++){
			// if(j<20)
			// 	std::cout << " 	j : " << j << std::endl;
			double dist = distance<Vector3d>(proj2topSurf, Flat_map[indices_map_top_surf[j]]);
			if (dist<2){
				associated_top_surf_node[i] = indices_map_top_surf[j];
				// std::cout << "Distance ok : " << dist << std::endl;
			}
		}
	}

	for(int i=0; i<m.Nb_vertices(); i++){
		// std::cout << i << " : " <<  associated_top_surf_node[i] << std::endl;
		m.Vertices[i].top_surf_indice = associated_top_surf_node[i];
		m.Vertices[i].weight = distFBS[i];
	}

}

void GridTransformation::Gaussian_random_field_K_initialisation(Parameters& param){
	K.resize(N.size(), N.size());
	Z.resize(N.size());
	Y.resize(N.size());

	for(int i=0; i< N.size(); i++){
		for(int j=0; j< N.size(); j++){
			double dist2 = (N[i](1)-N[j](1))*(N[i](1)-N[j](1))+(N[i](0)-N[j](0))*(N[i](0)-N[j](0));
			K(i,j) = pow(param.sigma,2) * exp( - dist2 / (2.0*pow(param.length,2)) );
			if(i==j){
				K(i,j) += 1e-8;
			}
		}
	}

	LLT<MatrixXd> lltOfA(K); // compute the Cholesky decomposition of K
	L = lltOfA.matrixL(); // retrieve factor L in the decomposition

	std::random_device rd;
	std::normal_distribution<> d{0,1};
	std::map<int, int> hist{};

	for(int n=0; n<N.size(); ++n) {
		Z(n) = d(rd);
		// ++hist[std::round(d(rd))];
		// std::cout << Z[n] << std::endl;
	}

	/* Plot the histogram */
	// for(int n=0; n<10000; ++n) {
    //     ++hist[std::round(d(rd))];
    // }

	// for(auto p : hist) {
	// 	std::cout << std::setw(2)
	// 				<< p.first << ' ' << std::string(p.second/200, '*') << '\n';
	// }

	Y = L*Z;
}

void GridTransformation::Apply_Gaussian_random_field(Vector3d& point, Vector3d& normal, double& value){

	// std::cout << Y(index) << std::endl;
	// std::cout << normal[0] << std::endl;
	// std::cout << normal[1] << std::endl;
	// std::cout << normal[2] << std::endl;

	point[0] += normal[0]*value;
	point[1] += normal[1]*value;
	point[2] += normal[2]*value;
}

void GridTransformation::Corner_thickness(Vector3d& point, Vector3d& normal, Parameters& param, std::tuple<double,double,double>& ramp_param, double CTV){

	double local_ymid = (local_ymin+local_ymax)/2.0;

	Vector3d ramp;
	ramp[2] = 0.0;
	ramp[0] = -std::get<0>(ramp_param);
	if(point[1]>=local_ymid)
		ramp[1] = -std::get<1>(ramp_param);
	else
		ramp[1] = std::get<2>(ramp_param);

	Vector3d init = point-ramp;

	double angle = 0.0;
	double sign=0.0;
	if (init[1]>=local_ymax){ // Flange + right corner
		if(init[0]>=local_xmax){ // right corner
			angle = atan((init[1]-local_ymax)/(init[0]-local_xmax));// * PI / 180.0;
			sign = 1;
		}
	}
	else if (init[1]<=local_ymin){ // Flange + left corner
		if(init[0]>=local_xmax){ // left corner
			angle = atan((init[1]-local_ymin)/(init[0]-local_xmax));// * PI / 180.0;
			sign = 1;
		}
	}


	double delta_val = - CTV * (1 - cos(4*angle))/2;

	// std::cout << "[ " << moved[0] << " ; " << moved[1] << " ]" << std::endl;

	point[0] += normal[0]*delta_val;
	point[1] += normal[1]*delta_val;
	point[2] += normal[2]*delta_val;

	return;
}

std::tuple<double,double,double> GridTransformation::ramp(Vector3d& point){
	double ymid = (ymax+ymin)/2.0;
	double decrease_0=0.0; double decrease_1=0.0; double increase_1=0.0;

	// std::vector<double> X = {0.0, 0.1*xmax, 0.75*xmax, xmax}; // must be increased

	if (point(2)>z1 && point(2)<z4){ // section 1, 2, 3

		if (point(2)>=z2 && point(2)<=z3){ // section 2

			double initial_0 = point(0);
			if (initial_0 > (xmax - threshold_avoiding_deformation_of_corners))
				decrease_0 = delta_max; // cnst
			else
				decrease_0 = initial_0 * delta_max / (xmax-threshold_avoiding_deformation_of_corners); // Linear

			point(0) -= decrease_0;

			double initial_1 = point(1);
			if(initial_1 > ymid) {

				if (initial_1 > (ymax - threshold_avoiding_deformation_of_corners))
					decrease_1 = delta_max;
				else
					decrease_1 = (delta_max/(ymax-threshold_avoiding_deformation_of_corners-ymid)) * (initial_1 - ymid); // Linear

				point(1) -= decrease_1;
			} else {

				if (initial_1 < ymin + threshold_avoiding_deformation_of_corners)
					increase_1 = delta_max;
				else
					increase_1 = (delta_max/(ymid-(ymin+threshold_avoiding_deformation_of_corners))) * (ymid - initial_1); // Linear

				point(1) += increase_1;
			}

		} else if (point(2)<z2){ // section 1

			double local_delta = (delta_max / (z2-z1)) * point(2) + ((z1*delta_max)/(z1-z2));

			double initial_0 = point(0);

			if (initial_0 > (xmax - threshold_avoiding_deformation_of_corners))
				decrease_0 = local_delta; // cnst
			else
				decrease_0 = initial_0 * local_delta / (xmax-threshold_avoiding_deformation_of_corners); // Linear


			point(0) -= decrease_0;

			double initial_1 = point(1);
			if(initial_1 > ymid) {

				if (initial_1 > (ymax - threshold_avoiding_deformation_of_corners))
					decrease_1 = local_delta;
				else
					decrease_1 = (local_delta/(ymax-threshold_avoiding_deformation_of_corners-ymid)) * (initial_1 - ymid); // Linear

				point(1) -= decrease_1;
			} else {

				if (initial_1 < ymin + threshold_avoiding_deformation_of_corners)
					increase_1 = local_delta;
				else
					increase_1 = (local_delta/(ymid-(ymin+threshold_avoiding_deformation_of_corners))) * (ymid - initial_1); // Linear

				point(1) += increase_1;
			}


		} else if (point(2)>z3){ // section 3

			double local_delta = (delta_max / (z3-z4)) * point(2) + ((z4*delta_max)/(z4-z3));

			double initial_0 = point(0);

			if (initial_0 > (xmax - threshold_avoiding_deformation_of_corners))
				decrease_0 = local_delta; // cnst
			else
				decrease_0 = initial_0 * local_delta / (xmax-threshold_avoiding_deformation_of_corners); // Linear

			point(0) -= decrease_0;

			double initial_1 = point(1);
			if(initial_1 > ymid) {

				if (initial_1 > (ymax - threshold_avoiding_deformation_of_corners))
					decrease_1 = local_delta;
				else
					decrease_1 = (local_delta/(ymax-threshold_avoiding_deformation_of_corners-ymid)) * (initial_1 - ymid); // Linear

				point(1) -= decrease_1;
			} else {
				if (initial_1 < ymin + threshold_avoiding_deformation_of_corners)
					increase_1 = local_delta;
				else
					increase_1 = (local_delta/(ymid-(ymin+threshold_avoiding_deformation_of_corners))) * (ymid - initial_1); // Linear

				point(1) += increase_1;
			}
		}
	}
	std::tuple<double, double, double> ramp_param = std::make_tuple(decrease_0,decrease_1,increase_1);
	return ramp_param;
}

void GridTransformation::Csection_wrinkles(Vector3d& point, Vector3d& normal, int number, Parameters& param, std::tuple<double,double,double>& ramp_param){

	double local_ymid = (local_ymin+local_ymax)/2.0;

	Vector3d ramp;
	ramp[2] = 0.0;
	ramp[0] = -std::get<0>(ramp_param);
	if(point[1]>=local_ymid)
		ramp[1] = -std::get<1>(ramp_param);
	else
		ramp[1] = std::get<2>(ramp_param);

	Vector3d init = point-ramp;

	defect_location = param.wrinklePos[number];
	wrinkleOri = param.wrinkleOri[number];
	defect_size = param.wrinkleSize[number];
	if (param.Shape == 0) { // C-spar
		damping = {param.wrinkleDamp[number](0)/(local_xmax-xmin), param.wrinkleDamp[number](1)/(local_ymax-local_ymin), param.wrinkleDamp[number](2)/(zmax-zmin)}; // the more damping is big the small is the wrinkle in that direction
	} else if (param.Shape == 1) { // Flat laminate
		damping = {param.wrinkleDamp[number](0), param.wrinkleDamp[number](1), param.wrinkleDamp[number](2)}; // the more damping is big the small is the wrinkle in that direction
	}

	Matrix3d B;
	B(0,0)=1.0;
	B(0,1)=0.0;
	B(0,2)=0.0;
	B(1,0)=0.0;
	B(1,1)=cos(wrinkleOri*PI/180.0);
	B(1,2)=sin(wrinkleOri*PI/180.0);
	B(2,0)=0.0;
	B(2,1)=-sin(wrinkleOri*PI/180.0);
	B(2,2)=cos(wrinkleOri*PI/180.0);

	Vector3d ref, moved;
	double dist_from_bottom_surf = 0.0;
	moved = init;
	// Moving in the flat coordinate system
	if (param.Shape == 0) { // C-spar
		if (init[1]>=local_ymax){
			if(init[0]<=local_xmax){
				dist_from_bottom_surf = (init[1] - local_ymax) - interior_radius;
			} else {
				dist_from_bottom_surf = sqrt((init[1]-local_ymax)*(init[1]-local_ymax)+(init[0]-local_xmax)*(init[0]-local_xmax))-interior_radius;
			}
			ref[0] = local_xmax + interior_radius + dist_from_bottom_surf;
			ref[1] = param.X + interior_radius;
			moved[0] = ref[0];
			double local_radius = interior_radius + dist_from_bottom_surf;

			if(init[0]<=local_xmax){
				double theta = PI/2.0;
				moved[1] = ref[1] + param.R*theta + abs(init[0]-ref[0])-local_radius;
			} else {
				double d = sqrt((init[1]-ref[1])*(init[1]-ref[1])+(init[0]-ref[0])*(init[0]-ref[0]));
				double theta = 2*asin(d/2.0/local_radius);

				moved[1] = ref[1] + theta*param.R;
			}
		}
		else if (init[1]<=local_ymin){
			if(init[0]<=local_xmax){
				dist_from_bottom_surf = - init[1] + local_ymin - interior_radius;
			} else {
				dist_from_bottom_surf = sqrt((init[1]-local_ymin)*(init[1]-local_ymin)+(init[0]-local_xmax)*(init[0]-local_xmax))-interior_radius;
			}
			ref[0] = local_xmax + interior_radius + dist_from_bottom_surf;
			ref[1] = interior_radius;
			moved[0] = ref[0];
			double local_radius = interior_radius + dist_from_bottom_surf;

			if(init[0]<=local_xmax){
				double theta = PI/2.0;
				moved[1] = ref[1] - (param.R*theta + abs(init[0]-ref[0])-local_radius);
			} else {
				double d = sqrt((init[1]-ref[1])*(init[1]-ref[1])+(init[0]-ref[0])*(init[0]-ref[0]));
				double theta = 2*asin(d/2.0/local_radius);

				moved[1] = ref[1] - theta*param.R;
			}
		}
	} else {
		// Nothing to do for the flat laminate
	}

	moved -= defect_location;
	Vector3d rot = B*moved;

	// std::cout << "––––––––" << std::endl;
	// std::cout << defect_size << std::endl;
	// std::cout << normal[0] << std::endl;
	// std::cout << normal[1] << std::endl;
	// std::cout << normal[2] << std::endl;
	// std::cout << damping[0] << std::endl;
	// std::cout << damping[1] << std::endl;
	// std::cout << damping[2] << std::endl;
	// std::cout << rot[0] << std::endl;
	// std::cout << rot[1] << std::endl;
	// std::cout << rot[2] << std::endl;

	point[0] += normal[0]*defect_size \
		* 1./pow(cosh(damping[0] * PI * rot[0]),2.)\
		* 1./pow(cosh(damping[1] * PI * rot[1]),2.)\
		* 1./pow(cosh(damping[2] * PI * rot[2]),2.);
	point[1] += normal[1]*defect_size \
		* 1./pow(cosh(damping[0] * PI * rot[0]),2.)\
		* 1./pow(cosh(damping[1] * PI * rot[1]),2.)\
		* 1./pow(cosh(damping[2] * PI * rot[2]),2.);
	point[2] += normal[2]*defect_size \
		* 1./pow(cosh(damping[0] * PI * rot[0]),2.)\
		* 1./pow(cosh(damping[1] * PI * rot[1]),2.)\
		* 1./pow(cosh(damping[2] * PI * rot[2]),2.);

	// if(abs(point[1]-init[1])>0.00001){
	// 	// std::cout << "Node: " <<  << std::endl;
	// 	std::cout << moved[0] << "; " << moved[1] << ", " << moved[2] << std::endl;
	// 	// std::cout << point[0] << "; " << point[1] << ", " << point[2] << std::endl;
	// }
	return;
}

void GridTransformation::Csection_wrinkles_2(Vector3d& point, Vector3d& normal, int number, Parameters& param, std::tuple<double,double,double>& ramp_param){

	double local_ymid = (local_ymin+local_ymax)/2.0;

	Vector3d ramp;
	ramp[2] = 0.0;
	ramp[0] = -std::get<0>(ramp_param);
	if(point[1]>=local_ymid)
		ramp[1] = -std::get<1>(ramp_param);
	else
		ramp[1] = std::get<2>(ramp_param);
	Vector3d init = point-ramp;

	defect_location = param.wrinklePos[number];
	wrinkleOri = param.wrinkleOri[number];
	defect_size = param.wrinkleSize[number];
	// if (param.Shape == 0) { // C-spar
	// 	damping = {param.wrinkleDamp[number](0)/(local_xmax-xmin), param.wrinkleDamp[number](1)/(local_ymax-local_ymin), param.wrinkleDamp[number](2)/(zmax-zmin)}; // the more damping is big the small is the wrinkle in that direction
	// } else if (param.Shape == 1) { // Flat laminate
	damping = {param.wrinkleDamp[number](0), param.wrinkleDamp[number](1), param.wrinkleDamp[number](2)}; // the more damping is big the small is the wrinkle in that direction
	// }

	Matrix3d B;
	B(0,0)=1.0;
	B(0,1)=0.0;
	B(0,2)=0.0;
	B(1,0)=0.0;
	B(1,1)=cos(wrinkleOri*PI/180.0);
	B(1,2)=sin(wrinkleOri*PI/180.0);
	B(2,0)=0.0;
	B(2,1)=-sin(wrinkleOri*PI/180.0);
	B(2,2)=cos(wrinkleOri*PI/180.0);

	Vector3d ref, moved;
	double dist_from_bottom_surf = 0.0;
	moved = init;
	// Moving in the flat coordinate system
	if (param.Shape == 0) { // C-spar
		if (init[1]>=local_ymax){
			if(init[0]<=local_xmax){
				dist_from_bottom_surf = (init[1] - local_ymax) - interior_radius;
			} else {
				dist_from_bottom_surf = sqrt((init[1]-local_ymax)*(init[1]-local_ymax)+(init[0]-local_xmax)*(init[0]-local_xmax))-interior_radius;
			}
			ref[0] = local_xmax + interior_radius + dist_from_bottom_surf;
			ref[1] = param.X + interior_radius;
			moved[0] = ref[0];
			double local_radius = interior_radius + dist_from_bottom_surf;

			if(init[0]<=local_xmax){
				double theta = PI/2.0;
				moved[1] = ref[1] + param.R*theta + abs(init[0]-ref[0])-local_radius;
			} else {
				double d = sqrt((init[1]-ref[1])*(init[1]-ref[1])+(init[0]-ref[0])*(init[0]-ref[0]));
				double theta = 2*asin(d/2.0/local_radius);

				moved[1] = ref[1] + theta*param.R;
			}
		}
		else if (init[1]<=local_ymin){
			if(init[0]<=local_xmax){
				dist_from_bottom_surf = - init[1] + local_ymin - interior_radius;
			} else {
				dist_from_bottom_surf = sqrt((init[1]-local_ymin)*(init[1]-local_ymin)+(init[0]-local_xmax)*(init[0]-local_xmax))-interior_radius;
			}
			ref[0] = local_xmax + interior_radius + dist_from_bottom_surf;
			ref[1] = interior_radius;
			moved[0] = ref[0];
			double local_radius = interior_radius + dist_from_bottom_surf;

			if(init[0]<=local_xmax){
				double theta = PI/2.0;
				moved[1] = ref[1] - (param.R*theta + abs(init[0]-ref[0])-local_radius);
			} else {
				double d = sqrt((init[1]-ref[1])*(init[1]-ref[1])+(init[0]-ref[0])*(init[0]-ref[0]));
				double theta = 2*asin(d/2.0/local_radius);

				moved[1] = ref[1] - theta*param.R;
			}
		}
	} else {
		// Nothing to do for the flat laminate
	}

	moved -= defect_location;
	Vector3d rot = B*moved;

	double x = moved[0]*cos((wrinkleOri*PI/180.0))+moved[2]*sin((wrinkleOri*PI/180.0));
	/* Non symetric */
	// defect_size  *= (1.0/PI) * (PI/2.0 + atan(2*x));

	/* Aguinea -- Serpentine cubic */
	double d = 2.0;
	// defect_size  *= d*x / (x*x + d*d);

	point[0] += normal[0]*defect_size \
		* 1./pow(cosh(damping[0] * PI * rot[0]),2.)\
		* 1./pow(cosh(damping[1] * PI * rot[1]),2.)\
		* 1./pow(cosh(damping[2] * PI * rot[2]),2.);
	point[1] += normal[1]*defect_size \
		* 1./pow(cosh(damping[0] * PI * rot[0]),2.)\
		* 1./pow(cosh(damping[1] * PI * rot[1]),2.)\
		* 1./pow(cosh(damping[2] * PI * rot[2]),2.);
	point[2] += normal[2]*defect_size \
		* 1./pow(cosh(damping[0] * PI * rot[0]),2.)\
		* 1./pow(cosh(damping[1] * PI * rot[1]),2.)\
		* 1./pow(cosh(damping[2] * PI * rot[2]),2.);




	// if(abs(point[1]-init[1])>0.00001){
	// 	// std::cout << "Node: " <<  << std::endl;
	// 	std::cout << moved[0] << "; " << moved[1] << ", " << moved[2] << std::endl;
	// 	// std::cout << point[0] << "; " << point[1] << ", " << point[2] << std::endl;
	// }
	return;
}

void GridTransformation::Csection_wrinkles_splined(Vector3d& point, Vector3d& normal, int number, Parameters& param, std::tuple<double,double,double>& ramp_param){

	double local_ymid = (local_ymin+local_ymax)/2.0;

	Vector3d ramp;
	ramp[2] = 0.0;
	ramp[0] = -std::get<0>(ramp_param);
	if(point[1]>=local_ymid)
		ramp[1] = -std::get<1>(ramp_param);
	else
		ramp[1] = std::get<2>(ramp_param);
	Vector3d init = point-ramp;

	defect_location = param.wrinklePos[number];

	Vector3d a(0.0, 0.0, 0.0);
	Vector3d b(0.0, -5.0, -6.0);
	Vector3d c(0.0, 10.0, -5.0);
	Vector3d e(0.0, 20.0, 10.0);
	Vector3d spline_point_1 = b;
	Vector3d spline_point_2 = c;
	Vector3d spline_point_3 = e;

	
	std::vector<double> X = {spline_point_1(1), a(1), spline_point_2(1), spline_point_3(1)};
	std::vector<double> Y = {spline_point_1(2), a(2), spline_point_2(2), spline_point_3(2)};

	Kluge::tk::spline s(X,Y);

	// std::vector<double> Values = arange<double>(spline_point_1(1), spline_point_3(1), 0.2);
	// std::cout << Values.size() << std::endl;
	// for (int i=0; i<Values.size();i++){
	// 	std::cout << Values[i]  << " : " << s(Values[i]) << std::endl;
	// }


	wrinkleOri = param.wrinkleOri[number];
	defect_size = param.wrinkleSize[number];
	// if (param.Shape == 0) { // C-spar
	// 	damping = {param.wrinkleDamp[number](0)/(local_xmax-xmin), param.wrinkleDamp[number](1)/(local_ymax-local_ymin), param.wrinkleDamp[number](2)/(zmax-zmin)}; // the more damping is big the small is the wrinkle in that direction
	// } else if (param.Shape == 1) { // Flat laminate
	damping = {param.wrinkleDamp[number](0), param.wrinkleDamp[number](1), param.wrinkleDamp[number](2)}; // the more damping is big the small is the wrinkle in that direction
	// }

	Matrix3d B;
	B(0,0)=1.0;
	B(0,1)=0.0;
	B(0,2)=0.0;
	B(1,0)=0.0;
	B(1,1)=cos(wrinkleOri*PI/180.0);
	B(1,2)=sin(wrinkleOri*PI/180.0);
	B(2,0)=0.0;
	B(2,1)=-sin(wrinkleOri*PI/180.0);
	B(2,2)=cos(wrinkleOri*PI/180.0);

	Vector3d ref, moved;
	double dist_from_bottom_surf = 0.0;
	moved = init;
	// Moving in the flat coordinate system
	if (param.Shape == 0) { // C-spar
		if (init[1]>=local_ymax){
			if(init[0]<=local_xmax){
				dist_from_bottom_surf = (init[1] - local_ymax) - interior_radius;
			} else {
				dist_from_bottom_surf = sqrt((init[1]-local_ymax)*(init[1]-local_ymax)+(init[0]-local_xmax)*(init[0]-local_xmax))-interior_radius;
			}
			ref[0] = local_xmax + interior_radius + dist_from_bottom_surf;
			ref[1] = param.X + interior_radius;
			moved[0] = ref[0];
			double local_radius = interior_radius + dist_from_bottom_surf;

			if(init[0]<=local_xmax){
				double theta = PI/2.0;
				moved[1] = ref[1] + param.R*theta + abs(init[0]-ref[0])-local_radius;
			} else {
				double d = sqrt((init[1]-ref[1])*(init[1]-ref[1])+(init[0]-ref[0])*(init[0]-ref[0]));
				double theta = 2*asin(d/2.0/local_radius);

				moved[1] = ref[1] + theta*param.R;
			}
		}
		else if (init[1]<=local_ymin){
			if(init[0]<=local_xmax){
				dist_from_bottom_surf = - init[1] + local_ymin - interior_radius;
			} else {
				dist_from_bottom_surf = sqrt((init[1]-local_ymin)*(init[1]-local_ymin)+(init[0]-local_xmax)*(init[0]-local_xmax))-interior_radius;
			}
			ref[0] = local_xmax + interior_radius + dist_from_bottom_surf;
			ref[1] = interior_radius;
			moved[0] = ref[0];
			double local_radius = interior_radius + dist_from_bottom_surf;

			if(init[0]<=local_xmax){
				double theta = PI/2.0;
				moved[1] = ref[1] - (param.R*theta + abs(init[0]-ref[0])-local_radius);
			} else {
				double d = sqrt((init[1]-ref[1])*(init[1]-ref[1])+(init[0]-ref[0])*(init[0]-ref[0]));
				double theta = 2*asin(d/2.0/local_radius);

				moved[1] = ref[1] - theta*param.R;
			}
		}
	} else {
		// Nothing to do for the flat laminate
	}

	moved -= defect_location;
	moved(2) += s(moved(1));
	// std::cout << Values[i]  << " : " << s(Values[i]) << std::endl;

	Vector3d rot = B*moved;

	double x = moved[0]*cos((wrinkleOri*PI/180.0))+moved[2]*sin((wrinkleOri*PI/180.0));
	/* Non symetric */
	// defect_size  *= (1.0/PI) * (PI/2.0 + atan(2*x));

	/* Aguinea -- Serpentine cubic */
	double d = 2.0;
	// defect_size  *= d*x / (x*x + d*d);

	point[0] += normal[0]*defect_size \
		* 1./pow(cosh(damping[0] * PI * rot[0]),2.)\
		* 1./pow(cosh(damping[1] * PI * rot[1]),2.)\
		* 1./pow(cosh(damping[2] * PI * rot[2]),2.);
	point[1] += normal[1]*defect_size \
		* 1./pow(cosh(damping[0] * PI * rot[0]),2.)\
		* 1./pow(cosh(damping[1] * PI * rot[1]),2.)\
		* 1./pow(cosh(damping[2] * PI * rot[2]),2.);
	point[2] += normal[2]*defect_size \
		* 1./pow(cosh(damping[0] * PI * rot[0]),2.)\
		* 1./pow(cosh(damping[1] * PI * rot[1]),2.)\
		* 1./pow(cosh(damping[2] * PI * rot[2]),2.);


	// if(abs(point[1]-init[1])>0.00001){
	// 	// std::cout << "Node: " <<  << std::endl;
	// 	std::cout << moved[0] << "; " << moved[1] << ", " << moved[2] << std::endl;
	// 	// std::cout << point[0] << "; " << point[1] << ", " << point[2] << std::endl;
	// }
	return;
}

void GridTransformation::wrinkles(Vector3d& point){

	// Wrinkle in plan XY
	Vector3d init = point;

	double xref = defect_location[0];
	double yref = defect_location[1];
	double zref = defect_location[2] + (ymid - init[1])*tan(wrinkleOri*PI/180.0);

	point[0] += defect_size \
		* 1./pow(cosh(damping[0] * PI * (init[0] - xref)),2.)\
		* 1./pow(cosh(damping[1] * PI * (init[1] - yref)),2.)\
		* 1./pow(cosh(damping[2] * PI * (init[2] - zref)),2.);

	point[1] += 0.5 * defect_size \
		* 1./pow(cosh(damping[0] * PI * (init[0] - xref)),2.)\
		* 1./pow(cosh(damping[1] * PI * (init[1] - yref)),2.)\
		* 1./pow(cosh(damping[2] * PI * (init[2] - zref)),2.);


	return;
}