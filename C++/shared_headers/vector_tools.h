#pragma once

#include <cmath>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Geometry>

using namespace Eigen;

Vector3f rotate_hors_axes(const Vector3f& u, float theta, const Vector3f& v) { // ref vector, angle, vector to be rot
	Vector3f urot;
	if (std::abs(theta)>0.0) {
		Matrix3f Q,R,id,Q2;
		Q(0,0)=0.0;
		Q(1,1)=0.0;
		Q(2,2)=0.0;
		Q(1,0)=-u(2);
		Q(2,0)=u(1);
		Q(2,1)=-u(0);
		Q(0,1)=-Q(1,0);
		Q(0,2)=-Q(2,0);
		Q(1,2)=-Q(2,1);
		id.setIdentity();
		Q2 = Q*Q;
		theta=theta*M_PI/180;
		R = id + std::sin(theta) * Q + (1-std::cos(theta)) * Q2;
		urot = R*v;
		urot.normalize();
		for (int i = 0; i < 3; i++) {
			if (std::abs(urot(i))<0.000001) {
				urot(i) = 0.0;
			}
		}
	} else {
		urot=v;
	}
	return urot;
}
