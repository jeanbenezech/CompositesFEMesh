#pragma once

#include "mesh.h"
#include "element.h"
#include "parameters.h"
#include "vector_tools.h"
#include "spline.h"
#include "gridTransformation.h"
#include "CTData.h"

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

				elem.tmp_U(0,l) = 1.0;
				elem.tmp_U(1,l) = 0.0;
				elem.tmp_U(2,l) = 0.0;
				elem.tmp_V(0,l) = 0.0;
				elem.tmp_V(1,l) = 1.0;
				elem.tmp_V(2,l) = 0.0;
				elem.tmp_W(0,l) = 0.0;
				elem.tmp_W(1,l) = 0.0;
				elem.tmp_W(2,l) = 1.0;

				for (int i=0; i<3; i++){
					elem.normal(i,l) = elem.V(i,l);
				}
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
					elem.tmp_U(i,l) = u(i);
					elem.tmp_V(i,l) = v(i);
					elem.tmp_W(i,l) = w(i);

					elem.normal(i,l) = v(i);
				}

				// orientation for the wrinkles
				if(param.add_wrinkles || param.GaussianThickness || param.CornerThickness)
					for(int nodeId=0; nodeId<elem.size; nodeId++)
						m.Vertices[elem.Nodes(nodeId,l)].normal += v;

			}
		}
	}

	// orientation for the wrinkles
	if(param.add_wrinkles || param.GaussianThickness || param.CornerThickness){
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

			elem.tmp_U(0,l) = 1.0;
			elem.tmp_U(1,l) = 0.0;
			elem.tmp_U(2,l) = 0.0;
			elem.tmp_V(0,l) = 0.0;
			elem.tmp_V(1,l) = 1.0;
			elem.tmp_V(2,l) = 0.0;
			elem.tmp_W(0,l) = 0.0;
			elem.tmp_W(1,l) = 0.0;
			elem.tmp_W(2,l) = 1.0;

			for (int i=0; i<3; i++){
				elem.normal(i,l) = elem.V(i,l);
			}

			// orientation for the wrinkles
			// if(param.add_wrinkles){
			for(int nodeId=0; nodeId<elem.size; nodeId++)
				m.Vertices[elem.Nodes(nodeId,l)].normal += elem.V.col(l);
			// }
		}
	}

	// orientation for the wrinkles
	// if(param.add_wrinkles)
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

	if(param.add_ramp==false && param.add_wrinkles==0 && param.rotateRVE==0){
		std::cout << "GeoTransformation is used for nothing." << std::endl;
		GridTransformation GT;
		GT.initialise(m, param);
		GT.finalise(m, param);
		return;
	}

	GridTransformation GT;
	GT.initialise(m, param);

	double no_ramp=0.0;
	std::tuple<double, double, double> ramp_param = std::make_tuple(no_ramp,no_ramp,no_ramp);

	GT.Flat_map.resize(m.Nb_vertices());
	GT.distFBS.resize(m.Nb_vertices());
	GT.weight.resize(m.Nb_vertices());
	GT.N.resize(0);
	if(param.GaussianThickness || param.CornerThickness){

		// #pragma omp parallel for schedule(guided)
		for(int node=0; node < m.Nb_vertices(); ++node) {
			Vector3d point;
			point(0) = m.Vertices[node].coord(0);
			point(1) = m.Vertices[node].coord(1);
			point(2) = m.Vertices[node].coord(2);

			GT.flat_transformation(node, point, param, ramp_param); // TODO: ramp_param is null here
		}
		GT.AssociatedtoTopSurfaceNode(m, param);
		if(param.GaussianThickness){
			GT.Gaussian_random_field_K_initialisation(param);
		}
	}


	// Interpolate the random field, defined on node, on elements
	// if(param.GaussianThickness || param.CornerThickness){
	if(param.GaussianThickness){
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
	// #pragma omp parallel for schedule(guided)
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
				// GT.Csection_wrinkles_2(point, normal, i, param, ramp_param);
				GT.Csection_wrinkles_splined(point, normal, i, param, ramp_param);
				// GT.Csection_wrinkles(point, normal, i, param, ramp_param);

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

		if(param.RotateFlanges){
			GT.RotateFlanges(point, param, ramp_param);
		}

		if(param.rotateRVE){
			GT.RotateRVE(point, param);
		}

		m.Vertices[node].coord(0) = point(0);
		m.Vertices[node].coord(1) = point(1);
		m.Vertices[node].coord(2) = point(2);
	}

	// // TRANSFORMATION OF LOCAL ORIENTATIONS ALL BUT ROTATE FLANGES
	for(auto& elem : m.Elements){
		for(int l=0; l < elem.nb; ++l) {

			Vector3d u, v, w, urot, vrot, wrot;
			for(int i=0; i<3;++i) {
				u(i)=elem.tmp_U(i,l);
				v(i)=elem.tmp_V(i,l);
				w(i)=elem.tmp_W(i,l);
			}
			urot=Vector3d::Zero();
			vrot=Vector3d::Zero();
			wrot=Vector3d::Zero();

			double h = 0.1;

			Vector3d normal = elem.normal.col(l);

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
					// GT.Csection_wrinkles_2(baryPu, normal, i, param, ramp_param_Pu);
					GT.Csection_wrinkles_splined(baryPu, normal, i, param, ramp_param_Pu);
					
					// GT.Csection_wrinkles_2(baryMu, normal, i, param, ramp_param_Mu);
					GT.Csection_wrinkles_splined(baryMu, normal, i, param, ramp_param_Mu);
					// GT.Csection_wrinkles(baryPu, normal, i, param, ramp_param_Pu);
					// GT.Csection_wrinkles(baryMu, normal, i, param, ramp_param_Mu);
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

			if(param.rotateRVE){
				GT.RotateRVE(baryPu, param);
				GT.RotateRVE(baryMu, param);
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
					// GT.Csection_wrinkles_2(baryPv, normal, i, param, ramp_param_Pv);
					GT.Csection_wrinkles_splined(baryMv, normal, i, param, ramp_param_Mv);
					// GT.Csection_wrinkles_2(baryMv, normal, i, param, ramp_param_Mv);
					GT.Csection_wrinkles_splined(baryMv, normal, i, param, ramp_param_Mv);
					// GT.Csection_wrinkles(baryPv, normal, i, param, ramp_param_Pv);
					// GT.Csection_wrinkles(baryMv, normal, i, param, ramp_param_Mv);
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

			if(param.rotateRVE){
				GT.RotateRVE(baryPv, param);
				GT.RotateRVE(baryMv, param);
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
					// GT.Csection_wrinkles_2(baryPw, normal, i, param, ramp_param_Pw);
					GT.Csection_wrinkles_splined(baryPw, normal, i, param, ramp_param_Pw);
					// GT.Csection_wrinkles_2(baryMw, normal, i, param, ramp_param_Mw);
					GT.Csection_wrinkles_splined(baryMw, normal, i, param, ramp_param_Mw);
					// GT.Csection_wrinkles(baryPw, normal, i, param, ramp_param_Pw);
					// GT.Csection_wrinkles(baryMw, normal, i, param, ramp_param_Mw);
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

			if(param.rotateRVE){
				GT.RotateRVE(baryPw, param);
				GT.RotateRVE(baryMw, param);
			}

			double angle0=0.0, angle1=0.0, angle2=0.0;

			Vector3d local_dir_use_U = baryPu-baryMu;
			angle1 = -atan2(local_dir_use_U.dot(u)/h,local_dir_use_U.dot(w)/h)*180.0/PI+90.0; // atan(dir0/dir2) --> adjacent dir is 2

			// Vector3d local_dir_use_V = baryPv-baryMv;
			angle2 = -atan2(local_dir_use_U.dot(v)/h,local_dir_use_U.dot(u)/h)*180.0/PI;

			Vector3d local_dir_use_W = baryPw-baryMw;
			angle0 = -atan2(local_dir_use_W.dot(w)/h,local_dir_use_W.dot(v)/h)*180.0/PI+90.0; // atan(dir2/dir1) --> adjacent dir is 1

			urot = rotate_hors_axes(v, angle1, u);
			urot = rotate_hors_axes(w, angle2, urot);

			vrot = rotate_hors_axes(u, angle0, v);
			vrot = rotate_hors_axes(w, angle2, vrot);

			wrot = rotate_hors_axes(u, angle0, w);
			wrot = rotate_hors_axes(v, angle1, wrot);

			for(int i=0; i<3;++i) {
				elem.tmp_U(i,l) = urot(i);
				elem.tmp_V(i,l) = vrot(i);
				elem.tmp_W(i,l) = wrot(i);
			}

		}
	}

	if(param.RotateFlanges){
		// TRANSFORMATION OF LOCAL ORIENTATIONS ROTATE FLANGES
		for(auto& elem : m.Elements){
			for(int l=0; l < elem.nb; ++l) {

				Vector3d u, v, w, urot, vrot, wrot;
				for(int i=0; i<3;++i) {
					u(i)=elem.tmp_U(i,l);
					v(i)=elem.tmp_V(i,l);
					w(i)=elem.tmp_W(i,l);
				}
				urot=Vector3d::Zero();
				vrot=Vector3d::Zero();
				wrot=Vector3d::Zero();

				double h = 0.1;

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

				GT.RotateFlanges(baryPu, param, ramp_param_Pu);
				GT.RotateFlanges(baryMu, param, ramp_param_Mu);

				//////////////////////////////////////////////////////////////////////////////////////////
				std::tuple<double, double, double> ramp_param_Pv = std::make_tuple(no_ramp,no_ramp,no_ramp);
				std::tuple<double, double, double> ramp_param_Mv = std::make_tuple(no_ramp,no_ramp,no_ramp);

				Vector3d baryPv = elem.Initial_barycenter.col(l) + h * v;
				Vector3d baryMv = elem.Initial_barycenter.col(l) - h * v;
				if(param.add_ramp){
					ramp_param_Pv = GT.ramp(baryPv);
					ramp_param_Mv = GT.ramp(baryMv);
				}

				GT.RotateFlanges(baryPv, param, ramp_param_Pv);
				GT.RotateFlanges(baryMv, param, ramp_param_Mv);

				//////////////////////////////////////////////////////////////////////////////////////////
				std::tuple<double, double, double> ramp_param_Pw = std::make_tuple(no_ramp,no_ramp,no_ramp);
				std::tuple<double, double, double> ramp_param_Mw = std::make_tuple(no_ramp,no_ramp,no_ramp);

				Vector3d baryPw = elem.Initial_barycenter.col(l) + h * w;
				Vector3d baryMw = elem.Initial_barycenter.col(l) - h * w;
				if(param.add_ramp){
					ramp_param_Pw = GT.ramp(baryPw);
					ramp_param_Mw = GT.ramp(baryMw);
				}


				GT.RotateFlanges(baryPw, param, ramp_param_Pw);
				GT.RotateFlanges(baryPw, param, ramp_param_Mw);

				double angle0=0.0, angle1=0.0, angle2=0.0;

				Vector3d local_dir_use_V = baryPv-baryMv;
				// Vector3d local_dir_use_W = baryPw-baryMw;
				angle0 = -atan2(local_dir_use_V.dot(w)/h,local_dir_use_V.dot(v)/h)*180.0/PI; // atan(dir2/dir1) --> adjacent dir is 1

				Vector3d local_dir_use_U = baryPu-baryMu;
				angle1 = -atan2(local_dir_use_U.dot(u)/h,local_dir_use_U.dot(w)/h)*180.0/PI + 90.0; // atan(dir0/dir2) --> adjacent dir is 2
				// Vector3d local_dir_use_V = baryPv-baryMv;
				angle2 = -atan2(local_dir_use_U.dot(v)/h,local_dir_use_U.dot(u)/h)*180.0/PI;


				urot = rotate_hors_axes(v, angle1, u);
				urot = rotate_hors_axes(w, angle2, urot);

				vrot = rotate_hors_axes(u, angle0, v);
				vrot = rotate_hors_axes(w, angle2, vrot);

				wrot = rotate_hors_axes(u, angle0, w);
				wrot = rotate_hors_axes(v, angle1, wrot);

				for(int i=0; i<3;++i) {
					// elem.U(i,l) = urot(i);
					// elem.V(i,l) = wrot(i);
					// elem.W(i,l) = -vrot(i);
					elem.tmp_U(i,l) = urot(i);
					elem.tmp_V(i,l) = vrot(i);
					elem.tmp_W(i,l) = wrot(i);
				}

			}
		}
	}

	for(auto& elem : m.Elements){
		elem.U= elem.tmp_U;
		elem.V= elem.tmp_V;
		elem.W= elem.tmp_W;
	}

	GT.finalise(m, param);
	return;
}

void CT_GeometricTransformation(Mesh& m, CT_DATA& ct, Parameters& param) {

	GridTransformation GT;
	GT.initialise(m, param);

	double no_ramp=0.0;
	std::tuple<double, double, double> ramp_param = std::make_tuple(no_ramp,no_ramp,no_ramp);

	// GT variable initiation
	GT.Flat_map.resize(m.Nb_vertices());
	GT.distFBS.resize(m.Nb_vertices());
	GT.weight.resize(m.Nb_vertices());
	GT.N.resize(0);
	// if(param.GaussianThickness || param.CornerThickness){

	for(int node=0; node < m.Nb_vertices(); ++node) {
		Vector3d point;
		point(0) = m.Vertices[node].coord(0);
		point(1) = m.Vertices[node].coord(1);
		point(2) = m.Vertices[node].coord(2);

		GT.prepare_GT_for_Flat_model(node, point, param, ramp_param); // TODO: ramp_param is null here
	}
	std::cout << "After the creation of Grid Transformation class." << std::endl;

	GT.AssociatedtoTopSurfaceNode(m, param);

	std::cout << "After pile connection creation." << std::endl;

	Mesh CTmesh;
	ct.init_mesh(CTmesh, GT, param);
	CTmesh.write_point_cloud_vtk("test", 1);
	// std::cout << CTmesh.Nb_vertices() << std::endl;

	// GT.Transfo_CT(CTmesh);
	// m.exportRandomFieldOnElement=true;
	for(auto& elem : m.Elements){
		for(int l=0; l < elem.nb; ++l) {
			int cnt=0;
			for(int i=0; i<elem.size; i++){
				auto iter = std::find (GT.indices_map_top_surf.begin(), GT.indices_map_top_surf.end(),GT.associated_top_surf_node[elem.Nodes(i,l)]);
				if(iter != GT.indices_map_top_surf.end()){
					int index = iter - GT.indices_map_top_surf.begin();

					elem.CT_Field[l].values.push_back(GT.Y(index));
					elem.CT_Field[l].nodes.push_back(m.Vertices[elem.Nodes(i,l)].coord);

					elem.CT_Field[l].average += GT.Y(index);
					cnt+=1;
				}
			}
			if(cnt>0)
				elem.CT_Field[l].average/=(cnt+0.0);
		}
	}


	// TRANSFORMATION OF NODES
	for(int node=0; node < m.Nb_vertices(); ++node) {

		Vector3d point;
		point(0) = m.Vertices[node].coord(0);
		point(1) = m.Vertices[node].coord(1);
		point(2) = m.Vertices[node].coord(2);

		Vector3d normal = m.Vertices[node].normal;

		// std::cout << normal(0) << ", " << normal(1) << ", " << normal(2) << std::endl;

		auto iter = std::find (GT.indices_map_top_surf.begin(), GT.indices_map_top_surf.end(), GT.associated_top_surf_node[node]);
		if(iter != GT.indices_map_top_surf.end()){
			int index = iter - GT.indices_map_top_surf.begin();
			// double tmp =pow(GT.weight[node],2)*GT.Y(index);
			// std::cout << GT.weight[node] << " |||| " << GT.Y(index) << std::endl;
			double tmp =GT.weight[node]*GT.Y(index);
			GT.Apply_Gaussian_random_field(point, normal, tmp);
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
			// if(param.add_ramp){
			// 	ramp_param_Pu = GT.ramp(baryPu);
			// 	ramp_param_Mu = GT.ramp(baryMu);
			// }
			// if(param.add_wrinkles){
			// 	for(int i=0; i<param.add_wrinkles;i++){
			// 		GT.Csection_wrinkles_2(baryPu, normal, i, param, ramp_param_Pu);
			// 		GT.Csection_wrinkles_2(baryMu, normal, i, param, ramp_param_Mu);
			// 		// GT.Csection_wrinkles(baryPu, normal, i, param, ramp_param_Pu);
			// 		// GT.Csection_wrinkles(baryMu, normal, i, param, ramp_param_Mu);
			// 	}
			// }

			// if(param.CornerThickness){
			// 	double tmp = local_weight*param.ThicknessVar;
			// 	GT.Corner_thickness(baryPu, normal, param, ramp_param_Pu, tmp);
			// 	GT.Corner_thickness(baryMu, normal, param, ramp_param_Mu, tmp);
			// }
		
			double value=0;
			double sum=0;
			for (int indexRF=0; indexRF<elem.CT_Field[l].nodes.size(); indexRF++){
				double dist = 1.0/distance<Vector3d>(baryPu,elem.CT_Field[l].nodes[indexRF]);
				value+= dist * elem.CT_Field[l].values[indexRF];
				sum+=dist;
			}
			value/=sum;
			value*=local_weight;

			GT.Apply_Gaussian_random_field(baryPu, normal, value);

			value=0;
			sum=0;
			for (int indexRF=0; indexRF<elem.CT_Field[l].nodes.size(); indexRF++){
				double dist = 1.0/distance<Vector3d>(baryMu,elem.CT_Field[l].nodes[indexRF]);
				value+= dist * elem.CT_Field[l].values[indexRF];
				sum+=dist;
			}
			value/=sum;
			value*=local_weight;
			GT.Apply_Gaussian_random_field(baryMu, normal, value);


			//////////////////////////////////////////////////////////////////////////////////////////
			std::tuple<double, double, double> ramp_param_Pv = std::make_tuple(no_ramp,no_ramp,no_ramp);
			std::tuple<double, double, double> ramp_param_Mv = std::make_tuple(no_ramp,no_ramp,no_ramp);

			Vector3d baryPv = elem.Initial_barycenter.col(l) + h * v;
			Vector3d baryMv = elem.Initial_barycenter.col(l) - h * v;

			value=0;
			sum=0;
			for (int indexRF=0; indexRF<elem.CT_Field[l].nodes.size(); indexRF++){
				double dist = 1.0/distance<Vector3d>(baryPv,elem.CT_Field[l].nodes[indexRF]);
				value+= dist * elem.CT_Field[l].values[indexRF];
				sum+=dist;
			}
			value/=sum;
			value*=local_weight;

			GT.Apply_Gaussian_random_field(baryPv, normal, value);

			value=0;
			sum=0;
			for (int indexRF=0; indexRF<elem.CT_Field[l].nodes.size(); indexRF++){
				double dist = 1.0/distance<Vector3d>(baryMv,elem.CT_Field[l].nodes[indexRF]);
				value+= dist * elem.CT_Field[l].values[indexRF];
				sum+=dist;
			}
			value/=sum;
			value*=local_weight;
			GT.Apply_Gaussian_random_field(baryMv, normal, value);



			//////////////////////////////////////////////////////////////////////////////////////////
			std::tuple<double, double, double> ramp_param_Pw = std::make_tuple(no_ramp,no_ramp,no_ramp);
			std::tuple<double, double, double> ramp_param_Mw = std::make_tuple(no_ramp,no_ramp,no_ramp);

			Vector3d baryPw = elem.Initial_barycenter.col(l) + h * w;
			Vector3d baryMw = elem.Initial_barycenter.col(l) - h * w;

			value=0;
			sum=0;
			for (int indexRF=0; indexRF<elem.CT_Field[l].nodes.size(); indexRF++){
				double dist = 1.0/distance<Vector3d>(baryPw,elem.CT_Field[l].nodes[indexRF]);
				value+= dist * elem.CT_Field[l].values[indexRF];
				sum+=dist;
			}
			value/=sum;
			value*=local_weight;

			GT.Apply_Gaussian_random_field(baryPw, normal, value);

			value=0;
			sum=0;
			for (int indexRF=0; indexRF<elem.CT_Field[l].nodes.size(); indexRF++){
				double dist = 1.0/distance<Vector3d>(baryMw,elem.CT_Field[l].nodes[indexRF]);
				value+= dist * elem.CT_Field[l].values[indexRF];
				sum+=dist;
			}
			value/=	sum;
			value*=local_weight;
			GT.Apply_Gaussian_random_field(baryMw, normal, value);

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

void Rigid_Boundary(Mesh& m, Parameters& param) {


	int beg, end;
	if (param.SizeIntervaldZ==2){
		beg = param.intervaldZ[0];
		end = param.intervaldZ[1];
	} else if (param.SizeIntervaldZ==4){
		beg = param.intervaldZ[0];
		end = param.intervaldZ[3];
	}

	// std::cout << "Steel behaviour from z=0 to z=" << beg << " and z=" << end << " to z=zmax" << std::endl;

	for(auto& elem : m.Elements){
		for(int l=0; l < elem.nb; ++l) {

			if (elem.Initial_barycenter.col(l)(2) < beg || 
				elem.Initial_barycenter.col(l)(2) > end){

				elem.Markers(1,l) = 4; // Marker for rigid boundary element
				elem.U(0,l) = 1.0;
				elem.U(1,l) = 0.0;
				elem.U(2,l) = 0.0;
				elem.V(0,l) = 0.0;
				elem.V(1,l) = 1.0;
				elem.V(2,l) = 0.0;
				elem.W(0,l) = 0.0;
				elem.W(1,l) = 0.0;
				elem.W(2,l) = 1.0;

			}

		}
	}
}


void TESTinitialiseDamage(Mesh& m, Parameters& param){ 
                // double locX = 100., double locZ = 370.,
                // double radius = 20., double radius_offset_small_axis = 0.,
                // double radius_offset_big_axis = 0.) {

    // Localisation initial delamination
    double radius = 20.;
    // radius_offset_small_axis = 0.;
    // radius_offset_big_axis = 0.;
    // locX = 100.;
    // double locZ = 370.;

    std::vector<double> locZ = param.centersZ;
    std::vector<double> locX = param.centersX;

    GridTransformation GT;
    GT.initialise(m,param);
    // double no_ramp=0.0;
    // std::tuple<double, double, double> ramp_param = std::make_tuple(no_ramp,no_ramp,no_ramp);


    // for(int node=0; node < m.Nb_vertices(); ++node) {
    //     Vector3d point;
    //     point(0) = m.Vertices[node].coord(0);
    //     point(1) = m.Vertices[node].coord(1);
    //     point(2) = m.Vertices[node].coord(2);
    //     GT.prepare_GT_for_Flat_model(node, point, param, ramp_param); // TODO: ramp_param is null here
    // }



    for(auto& elem : m.Elements){
        for(int l=0; l < elem.nb; ++l) {

            Vector3d c = elem.center(m.Vertices,l);

            // Vector3d point = GT.TESTtoFlatCoordinateSys(c, param,0,0); // unShrink along X and Y = true
            // Vector3d point = GT.TESTtoFlatCoordinateSys(c, param,1,1); // unShrink along X and Y = true
            Vector3d point = GT.toFlatCoordinateSys_VCarl(c, param,1,1); // unShrink along X and Y = true
            // std::cout << point[0] << std::endl;
            // std::cout << point[1] << std::endl;
            // std::cout << point[2] << std::endl;

            int plan_width_coord = 1;
            for (int iter_locZ=0;  iter_locZ<locZ.size();iter_locZ++ ){
                for (int iter_locX=0;  iter_locX<locX.size();iter_locX++ ){
                    /* Circular initial delam */
                    double dist=  pow(locX[iter_locX]-point[plan_width_coord],2) / pow(radius,2)
                                + pow(locZ[iter_locZ]-point[2],2)  / pow(radius,2);
                    // double dist=  pow(locX[iter_locX+6*iter_locZ]-point[plan_width_coord],2) / pow(radius,2)
                    //             + pow(locZ[iter_locZ]-point[2],2)  / pow(radius,2);
                    // std::cout << dist << std::endl;


                    if (dist <= 1 ){
                        elem.preDamage(0, l) = 1;
                        // std::cout << "ici" << std::endl;
                    } else {
                        // elem.preDamage(0, l) = 0;
                    }
                }
            }
        } // endfor
    } // endfor

	if (param.do_flatten){
		for(int node=0; node < m.Nb_vertices(); ++node) {
			Vector3d point;
			point(0) = m.Vertices[node].coord(0);
			point(1) = m.Vertices[node].coord(1);
			point(2) = m.Vertices[node].coord(2);
			// m.Vertices[node].coord = GT.TESTtoFlatCoordinateSys(point, param,1,1);
			m.Vertices[node].coord = GT.toFlatCoordinateSys_VCarl(point, param,1,1);
		}
	}
} //end initialiseDamage

