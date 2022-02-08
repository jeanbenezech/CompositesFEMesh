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

	// General parameters
	int verbosity = 0;

	void initialise(Mesh& m, Parameters& param);
	std::tuple<double,double,double> ramp(Vector3d& point);
	void Csection_wrinkles(Vector3d& point, Vector3d& normal, int number, Parameters& param, std::tuple<double,double,double>& ramp_param);

	void flat_transformation(int& node, Vector3d& point, Parameters& param, std::tuple<double,double,double>& ramp_param);
	void AssociatedtoTopSurfaceNode(Mesh& m, Parameters& param);
	void Gaussian_random_field_K_initialisation(Parameters& param);
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
		std::cout << "Error: Cannot open file" << filename << std::endl;
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

void GridTransformation::AssociatedtoTopSurfaceNode(Mesh& m, Parameters& param){

	int nb_vertice_top_plan = m.Nb_vertices()/(m.Nb_plies()+1);
	Vector3d ext (1.0, 0.0, 0.0);

	associated_top_surf_node.resize(m.Nb_vertices());

	// std::cout << " size : " << nb_vertice_top_plan << std::endl;
	// std::cout << " Flat_map size : " << Flat_map.size() << std::endl;
	// std::cout << " distFBS size : " << distFBS.size() << std::endl;

	for(int i=0; i<m.Nb_vertices(); i++){
		Vector3d proj2topSurf = Flat_map[i] + (param.Y-distFBS[i])*ext;
		// std::cout << "i : " << i << std::endl;
		for(int j=0; j<indices_map_top_surf.size(); j++){
			// std::cout << " 	j : " << j << std::endl;
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

				if (initial_1 < threshold_avoiding_deformation_of_corners)
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
	damping = {param.wrinkleDamp[number](0)/(local_xmax-xmin), param.wrinkleDamp[number](1)/(local_ymax-local_ymin), param.wrinkleDamp[number](2)/(zmax-zmin)}; // the more damping is big the small is the wrinkle in that direction

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

	moved -= defect_location;
	Vector3d rot = B*moved;

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

/////////////////////////////////////////////////////////////////////////////////////

void find_U(Parameters& param, Vector3d& u, Vector3d bary) {
	Vector3d Z (0.0, 0.0, 1.0);
	Vector3d Z_ (0.0, 0.0, -1.0);

	double xmin = param.Height;
	double xmax = param.Height;

	double ymin = param.R;//-param.Y;
	double ymax = param.X+ymin;

	if (bary(0) < xmin) {
		if (bary(1) < ymin) {
			u(0) = 1.0;
		} else if (bary(1) > ymax) {
			u(0) = -1.0;
		}
	} else {
		if (ymin <= bary(1) && bary(1) <= ymax) {
			u(1) = 1.0;
		} else if (bary(1) < ymin) {
			Vector3d tmp (0.0, 0.0, 0.0);
			// tmp(0) = std::abs(bary(0) - xmin);
			// tmp(1) = std::abs(bary(1) - ymin);
			tmp(0) = bary(0) - xmin;
			tmp(1) = bary(1) - ymin;
			tmp.normalize();
			u=tmp.cross(Z_);
		} else if (bary(1) > ymax) {
			Vector3d tmp (0.0, 0.0, 0.0);
			tmp(0) = bary(0) - xmax;
			tmp(1) = bary(1) - ymax;
			tmp.normalize();
			u=tmp.cross(Z_);
		}
	}
	return;
}

void localCoorSyst(Mesh& m, Parameters& param) {

	for(auto& elem : m.Elements){ // Loop over element types
		for(int l=0; l < elem.nb; ++l) { // Loop over elements for each type

			Vector3d u, v, urot, w;
			u=Vector3d::Zero();
			v=Vector3d::Zero();
			w=Vector3d::Zero();
			urot=Vector3d::Zero();

			if(elem.Markers(0,l) == param.resin_id){ // Resin layer got global coord sys
				elem.U(0,l) = 1.0;
				elem.U(1,l) = 0.0;
				elem.U(2,l) = 0.0;

				elem.V(0,l) = 0.0;
				elem.V(1,l) = 1.0;
				elem.V(2,l) = 0.0;

				elem.W(0,l) = 0.0;
				elem.W(1,l) = 0.0;
				elem.W(2,l) = 1.0;

			} else {

				find_U(param, u, elem.Initial_barycenter.col(l));

				Vector3d Z (0.0, 0.0, 1.0);

				// u.normalize();
				v = u.cross(Z);
				v.normalize();
				w=u.cross(v);
				w.normalize();

				for(int i=0; i<3;++i) {
					elem.U(i,l) = u(i);
					elem.V(i,l) = v(i);
					elem.W(i,l) = w(i);
				}

				// orientation for the wrinkles
				if(param.add_wrinkles || param.GaussianThickness)
					for(int nodeId=0; nodeId<elem.size; nodeId++)
						m.Vertices[elem.Nodes(nodeId,l)].normal += v;

			}
		}
	}

	// orientation for the wrinkles
	if(param.add_wrinkles || param.GaussianThickness){
		for(int nodeId=0; nodeId<m.Vertices.size(); nodeId++){
			m.Vertices[nodeId].normal.normalize();
		}
	}

}

void globalCoorSyst(Mesh& m, Parameters& param) {

	for(auto& elem : m.Elements){
		for(int l=0; l < elem.nb; ++l) {

			elem.U(0,l) = 1.0;
			elem.U(1,l) = 0.0;
			elem.U(2,l) = 0.0;

			elem.V(0,l) = 0.0;
			elem.V(1,l) = 1.0;
			elem.V(2,l) = 0.0;

			elem.W(0,l) = 0.0;
			elem.W(1,l) = 0.0;
			elem.W(2,l) = 1.0;

			// orientation for the wrinkles
			if(param.add_wrinkles){
				for(int nodeId=0; nodeId<elem.size; nodeId++)
					m.Vertices[elem.Nodes(nodeId,l)].normal += elem.V.col(l);
			}
		}
	}

	// orientation for the wrinkles
	if(param.add_wrinkles)
		for(int nodeId=0; nodeId<m.Vertices.size(); nodeId++)
			m.Vertices[nodeId].normal.normalize();

}

void attribute_weight(Mesh& m, Parameters& param) {
	double xmax = param.Height;
	double ymin = param.R-param.Y;
	double ymax = param.X+ymin;
	double ymid = (ymin+ymax)/2.0;

	for(auto& elem : m.Elements){

		double min_weight = 0;
		double max_weight = 0;
		std::vector<double> weight(elem.nb);

		for(int l=0; l < elem.nb; ++l) {
			Vector3d init=Vector3d::Zero();
			for(int i=0; i<elem.size; i++)
				init += m.Vertices[elem.Nodes(i,l)].coord;
			init *= 1.0 / (elem.size+0.0);


			Vector3d ref, moved;
			double dist_from_bottom_surf;
			double interior_radius = param.R - param.Y;
			moved = init;
			if (init[1]>=ymax){
				if(init[0]<=xmax){
					dist_from_bottom_surf = (init[1] - ymax) - interior_radius;
				} else {
					dist_from_bottom_surf = sqrt((init[1]-ymax)*(init[1]-ymax)+(init[0]-xmax)*(init[0]-xmax))-interior_radius;
				}
				ref[0] = param.Height + interior_radius + dist_from_bottom_surf;
				ref[1] = param.X + interior_radius;
				moved[0] = ref[0];
				double local_radius = interior_radius + dist_from_bottom_surf;

				if(init[0]<=xmax){
					double theta = PI/2.0;
					moved[1] = ref[1] + param.R*theta + abs(init[0]-ref[0])-local_radius;
				} else {
					double d = sqrt((init[1]-ref[1])*(init[1]-ref[1])+(init[0]-ref[0])*(init[0]-ref[0]));
					double theta = 2*asin(d/2.0/local_radius);

					moved[1] = ref[1] + theta*param.R;
				}
			}
			else if (init[1]<=ymin){
				if(init[0]<=xmax){
					dist_from_bottom_surf = - init[1] + ymin - interior_radius;
				} else {
					dist_from_bottom_surf = sqrt((init[1]-ymin)*(init[1]-ymin)+(init[0]-xmax)*(init[0]-xmax))-interior_radius;
				}
				ref[0] = param.Height + interior_radius + dist_from_bottom_surf;
				ref[1] = interior_radius;
				moved[0] = ref[0];
				double local_radius = interior_radius + dist_from_bottom_surf;

				if(init[0]<=xmax){
					double theta = PI/2.0;
					moved[1] = ref[1] - (param.R*theta + abs(init[0]-ref[0])-local_radius);
				} else {
					double d = sqrt((init[1]-ref[1])*(init[1]-ref[1])+(init[0]-ref[0])*(init[0]-ref[0]));
					double theta = 2*asin(d/2.0/local_radius);

					moved[1] = ref[1] - theta*param.R;
				}
			}

			weight[l] = moved[1]*init[2];
			if (weight[l]<min_weight)
				min_weight = weight[l];
			if (weight[l]>max_weight)
				max_weight = weight[l];

		}

		for(int l=0; l < elem.nb; ++l) {
			weight[l] -= min_weight;
			weight[l] /= (max_weight - min_weight);
		}

		std::vector<double> conversion(0);
		conversion.push_back(weight[0]);
		for(int l=1; l < elem.nb; ++l) {
			bool is_in = false;
			for(int iter = 0; iter<conversion.size();iter++){
				if(abs(weight[l]-conversion[iter]) < 0.00001 ){
					is_in=true;
					break;
				}
			}
			if(!is_in)
				conversion.push_back(weight[l]);
		}


		std::random_shuffle ( conversion.begin(), conversion.end() );


		for(int l=0; l < elem.nb; ++l) {
			for(int iter = 0; iter<conversion.size();iter++){
				if(abs(weight[l]-conversion[iter]) < 0.00001 ){
					elem.DD_weight(0,l) = 100*iter;
					break;
				}
			}
		}
	}
}

void GeometricTransformation(Mesh& m, Parameters& param) {

	if(param.add_ramp==false && param.add_wrinkles==0){
		std::cout << "GeoTransformation is used for nothing." << std::endl;
		return;
	}

	GridTransformation GT;
	GT.initialise(m, param);

	double no_ramp=0.0;
	std::tuple<double, double, double> ramp_param = std::make_tuple(no_ramp,no_ramp,no_ramp);

	if(param.GaussianThickness || param.CornerThickness){
		GT.Flat_map.resize(m.Nb_vertices());
		GT.distFBS.resize(m.Nb_vertices());
		GT.weight.resize(m.Nb_vertices());
		GT.N.resize(0);
		for(int node=0; node < m.Nb_vertices(); ++node) {
			Vector3d point;
			point(0) = m.Vertices[node].coord(0);
			point(1) = m.Vertices[node].coord(1);
			point(2) = m.Vertices[node].coord(2);

			GT.flat_transformation(node, point, param, ramp_param); // TODO: ramp_param is null here
		}
		GT.AssociatedtoTopSurfaceNode(m, param);
		GT.Gaussian_random_field_K_initialisation(param);
	}


	// Interpolate the random field, defined on node, on elements
	if(param.GaussianThickness || param.CornerThickness){
		m.exportRandomFieldOnElement=true;
		for(auto& elem : m.Elements){
			for(int l=0; l < elem.nb; ++l) {
				int cnt=0;
				for(int i=0; i<elem.size; i++){
					auto iter = std::find (GT.indices_map_top_surf.begin(), GT.indices_map_top_surf.end(),GT.associated_top_surf_node[elem.Nodes(i,l)]);
					if(iter != GT.indices_map_top_surf.end()){
						int index = iter - GT.indices_map_top_surf.begin();

						elem.Random_Field[l].values.push_back(GT.Y(index));
						elem.Random_Field[l].nodes.push_back(m.Vertices[elem.Nodes(i,l)].coord);

						elem.Random_Field[l].average += GT.Y(index);
						cnt+=1;
					}
				}
				if(cnt>0)
					elem.Random_Field[l].average/=(cnt+0.0);
			}
		}
	}

	// TRANSFORMATION OF NODES
	for(int node=0; node < m.Nb_vertices(); ++node) {

		Vector3d point;
		point(0) = m.Vertices[node].coord(0);
		point(1) = m.Vertices[node].coord(1);
		point(2) = m.Vertices[node].coord(2);

		Vector3d normal = m.Vertices[node].normal;

		if(param.add_ramp)
			ramp_param = GT.ramp(point);
		if(param.add_wrinkles)
			for(int i=0; i<param.add_wrinkles;i++)
				GT.Csection_wrinkles(point, normal, i, param, ramp_param);

		if(param.CornerThickness){
			// double tmp = pow(GT.weight[node],2)*param.ThicknessVar;
			double tmp = GT.weight[node]*param.ThicknessVar;
			GT.Corner_thickness(point, normal, param, ramp_param, tmp);
		}

		if(param.GaussianThickness){
			auto iter = std::find (GT.indices_map_top_surf.begin(), GT.indices_map_top_surf.end(), GT.associated_top_surf_node[node]);
			if(iter != GT.indices_map_top_surf.end()){
				int index = iter - GT.indices_map_top_surf.begin();	
				// double tmp =pow(GT.weight[node],2)*GT.Y(index);
				double tmp =GT.weight[node]*GT.Y(index);
				GT.Apply_Gaussian_random_field(point, normal, tmp);
			}
		}

		m.Vertices[node].coord(0) = point(0);
		m.Vertices[node].coord(1) = point(1);
		m.Vertices[node].coord(2) = point(2);
	}

	// TRANSFORMATION OF LOCAL ORIENTATIONS
	for(auto& elem : m.Elements){
		for(int l=0; l < elem.nb; ++l) {

			Vector3d u, v, w, urot, vrot, wrot;
			for(int i=0; i<3;++i) {
				u(i)=elem.U(i,l);
				v(i)=elem.V(i,l);
				w(i)=elem.W(i,l);
			}
			urot=Vector3d::Zero();
			vrot=Vector3d::Zero();
			wrot=Vector3d::Zero();

			double h = 0.1;

			Vector3d normal = elem.V.col(l);

			double local_weight = 0.0;
			for(int i=0; i<elem.size; i++){
				local_weight += GT.weight[elem.Nodes(i,l)];
			}
			local_weight/=(elem.size+0.0);

			//////////////////////////////////////////////////////////////////////////////////////////
			std::tuple<double, double, double> ramp_param_Pu = std::make_tuple(no_ramp,no_ramp,no_ramp);
			std::tuple<double, double, double> ramp_param_Mu = std::make_tuple(no_ramp,no_ramp,no_ramp);

			Vector3d baryPu = elem.Initial_barycenter.col(l) + h * u;
			Vector3d baryMu = elem.Initial_barycenter.col(l) - h * u;
			if(param.add_ramp){
				ramp_param_Pu = GT.ramp(baryPu);
				ramp_param_Mu = GT.ramp(baryMu);
			}
			if(param.add_wrinkles){
				for(int i=0; i<param.add_wrinkles;i++){
					GT.Csection_wrinkles(baryPu, normal, i, param, ramp_param_Pu);
					GT.Csection_wrinkles(baryMu, normal, i, param, ramp_param_Mu);
				}
			}

			if(param.CornerThickness){
				double tmp = local_weight*param.ThicknessVar;
				GT.Corner_thickness(baryPu, normal, param, ramp_param_Pu, tmp);
				GT.Corner_thickness(baryMu, normal, param, ramp_param_Mu, tmp);
			}
			if(param.GaussianThickness){
				double value=0;
				double sum=0;
				for (int indexRF=0; indexRF<elem.Random_Field[l].nodes.size(); indexRF++){
					double dist = 1.0/distance<Vector3d>(baryPu,elem.Random_Field[l].nodes[indexRF]);
					value+= dist * elem.Random_Field[l].values[indexRF];
					sum+=dist;
				}
				value/=sum;
				value*=local_weight;

				GT.Apply_Gaussian_random_field(baryPu, normal, value);

				value=0;
				sum=0;
				for (int indexRF=0; indexRF<elem.Random_Field[l].nodes.size(); indexRF++){
					double dist = 1.0/distance<Vector3d>(baryMu,elem.Random_Field[l].nodes[indexRF]);
					value+= dist * elem.Random_Field[l].values[indexRF];
					sum+=dist;
				}
				value/=sum;
				value*=local_weight;
				GT.Apply_Gaussian_random_field(baryMu, normal, value);
			}

			//////////////////////////////////////////////////////////////////////////////////////////
			std::tuple<double, double, double> ramp_param_Pv = std::make_tuple(no_ramp,no_ramp,no_ramp);
			std::tuple<double, double, double> ramp_param_Mv = std::make_tuple(no_ramp,no_ramp,no_ramp);

			Vector3d baryPv = elem.Initial_barycenter.col(l) + h * v;
			Vector3d baryMv = elem.Initial_barycenter.col(l) - h * v;
			if(param.add_ramp){
				ramp_param_Pv = GT.ramp(baryPv);
				ramp_param_Mv = GT.ramp(baryMv);
			}
			if(param.add_wrinkles){
				for(int i=0; i<param.add_wrinkles;i++){
					GT.Csection_wrinkles(baryPv, normal, i, param, ramp_param_Pv);
					GT.Csection_wrinkles(baryMv, normal, i, param, ramp_param_Mv);
				}
			}

			if(param.CornerThickness){
				double tmp = local_weight*param.ThicknessVar;
				GT.Corner_thickness(baryPv, normal, param, ramp_param_Pv, tmp);
				GT.Corner_thickness(baryMv, normal, param, ramp_param_Mv, tmp);
			}
			if(param.GaussianThickness){
				double value=0;
				double sum=0;
				for (int indexRF=0; indexRF<elem.Random_Field[l].nodes.size(); indexRF++){
					double dist = 1.0/distance<Vector3d>(baryPv,elem.Random_Field[l].nodes[indexRF]);
					value+= dist * elem.Random_Field[l].values[indexRF];
					sum+=dist;
				}
				value/=sum;
				value*=local_weight;

				GT.Apply_Gaussian_random_field(baryPv, normal, value);

				value=0;
				sum=0;
				for (int indexRF=0; indexRF<elem.Random_Field[l].nodes.size(); indexRF++){
					double dist = 1.0/distance<Vector3d>(baryMv,elem.Random_Field[l].nodes[indexRF]);
					value+= dist * elem.Random_Field[l].values[indexRF];
					sum+=dist;
				}
				value/=sum;
				value*=local_weight;
				GT.Apply_Gaussian_random_field(baryMv, normal, value);
			}


			//////////////////////////////////////////////////////////////////////////////////////////
			std::tuple<double, double, double> ramp_param_Pw = std::make_tuple(no_ramp,no_ramp,no_ramp);
			std::tuple<double, double, double> ramp_param_Mw = std::make_tuple(no_ramp,no_ramp,no_ramp);

			Vector3d baryPw = elem.Initial_barycenter.col(l) + h * w;
			Vector3d baryMw = elem.Initial_barycenter.col(l) - h * w;
			if(param.add_ramp){
				ramp_param_Pw = GT.ramp(baryPw);
				ramp_param_Mw = GT.ramp(baryMw);
			}
			if(param.add_wrinkles){
				for(int i=0; i<param.add_wrinkles;i++){
					GT.Csection_wrinkles(baryPw, normal, i, param, ramp_param_Pw);
					GT.Csection_wrinkles(baryMw, normal, i, param, ramp_param_Mw);
				}
			}

			
			if(param.CornerThickness){
				double tmp = local_weight*param.ThicknessVar;
				GT.Corner_thickness(baryPw, normal, param, ramp_param_Pw, tmp);
				GT.Corner_thickness(baryMw, normal, param, ramp_param_Mw, tmp);
			}
			if(param.GaussianThickness){
				double value=0;
				double sum=0;
				for (int indexRF=0; indexRF<elem.Random_Field[l].nodes.size(); indexRF++){
					double dist = 1.0/distance<Vector3d>(baryPw,elem.Random_Field[l].nodes[indexRF]);
					value+= dist * elem.Random_Field[l].values[indexRF];
					sum+=dist;
				}
				value/=sum;
				value*=local_weight;

				GT.Apply_Gaussian_random_field(baryPw, normal, value);

				value=0;
				sum=0;
				for (int indexRF=0; indexRF<elem.Random_Field[l].nodes.size(); indexRF++){
					double dist = 1.0/distance<Vector3d>(baryMw,elem.Random_Field[l].nodes[indexRF]);
					value+= dist * elem.Random_Field[l].values[indexRF];
					sum+=dist;
				}
				value/=sum;
				value*=local_weight;
				GT.Apply_Gaussian_random_field(baryMw, normal, value);
			}

			double angle0=0.0, angle1=0.0, angle2=0.0;

			// No shift in the V direction which is the normal to the ply
			// if(dim==1){
			Vector3d local_dir_use_U = baryPu-baryMu;
			angle1 = -atan2(local_dir_use_U.dot(u)/h,local_dir_use_U.dot(w)/h)*180.0/PI+90.0; // atan(dir0/dir2) --> adjacent dir is 2

			// Vector3d local_dir_use_V = baryPv-baryMv;
			angle2 = -atan2(local_dir_use_U.dot(v)/h,local_dir_use_U.dot(u)/h)*180.0/PI;

			Vector3d local_dir_use_W = baryPw-baryMw;
			angle0 = -atan2(local_dir_use_W.dot(w)/h,local_dir_use_W.dot(v)/h)*180.0/PI+90.0; // atan(dir2/dir1) --> adjacent dir is 1

			// ///
			// Vector3d local_dir_use_U = baryPu-baryMu;
			// angle1 = -atan2(local_dir_use_U.dot(u)/h,local_dir_use_U.dot(w)/h)*180.0/PI+90.0; // atan(dir0/dir2) --> adjacent dir is 2

			// ///
			// Vector3d local_dir_use_V = baryPv-baryMv;
			// angle2 = -atan2(local_dir_use_V.dot(v)/h,local_dir_use_V.dot(u)/h)*180.0/PI+90.0; // atan(dir1/dir0) --> adjacent dir is 0

			// ///
			// Vector3d local_dir_use_W = baryPw-baryMw;
			// angle0 = -atan2(local_dir_use_W.dot(w)/h,local_dir_use_W.dot(v)/h)*180.0/PI+90.0; // atan(dir2/dir1) --> adjacent dir is 1

			// if (angle1 > 0.2)
			// 	std::cout << "angle0 : " << angle0 << ", angle1 : " << angle1 << ", angle2 : " << angle2 << std::endl;

			urot = rotate_hors_axes(v, angle1, u);
			urot = rotate_hors_axes(w, angle2, urot);

			vrot = rotate_hors_axes(u, angle0, v);
			vrot = rotate_hors_axes(w, angle2, vrot);

			wrot = rotate_hors_axes(u, angle0, w);
			wrot = rotate_hors_axes(v, angle1, wrot);

			for(int i=0; i<3;++i) {
				elem.U(i,l) = urot(i);
				elem.V(i,l) = vrot(i);
				elem.W(i,l) = wrot(i);
			}
		}
	}

	return;
}

void StackingSequence(Mesh& m, Parameters& param) {

	for(auto& elem : m.Elements){
		for(int l=0; l < elem.nb; ++l) {

			if(elem.Markers(0,l) != param.resin_id && elem.Markers(0,l) != param.cz_id){ // Resin layer or CZ do not need StaqSeq
				Vector3d u, v, w, urot, wrot;
				for(int i=0; i<3;++i) {
					u(i)=elem.U(i,l);
					v(i)=elem.V(i,l);
					w(i)=elem.W(i,l);
				}
				urot=Vector3d::Zero();
				wrot=Vector3d::Zero();

				Vector3d Z (0.0, 0.0, 1.0);

				urot = rotate_hors_axes(v, m.stacking_sequence[elem.Markers(0,l)-1], u);

				wrot=urot.cross(v);
				wrot.normalize();

				for(int i=0; i<3;++i) {
					elem.U(i,l) = urot(i);
					elem.W(i,l) = wrot(i);
				}
			}
		}
	}
}