/*
 * vector3.cpp
 *
 *  Created on: Feb 12, 2023
 *      Author: norms
 */

// LINUX
#include "../headers/vector3.h"

// WINDOWS
//#include "..\headers\vector3.h"

#include "math.h"

double dotProd(vector3 v1, vector3 v2)
{
	return (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);
}

// Vector cross product
vector3 crossProd(vector3 v1, vector3 v2)
{
	return (vector3){v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x};
}

double mag(vector3 v)
{
	return sqrt((v.x * v.x) + (v.y * v.y) + (v.z * v.z));
}

// Returns unit vector of a vector3
vector3 unitVector(vector3 v)
{
	double vmag = mag(v);
	vector3 uv = {v.x / vmag, v.y / vmag, v.z / vmag};
	return uv;
}

// Returns angle between two vector3s
double angle(vector3 v1, vector3 v2)
{
	return acos(dotProd(v1, v2) / (mag(v1) * (mag(v2))));
}

// Returns angle between unit vectors of v1 and v2, good if worried about overflow
double unitAngle(vector3 v1, vector3 v2)
{
	vector3 uv1 = unitVector(v1); vector3 uv2 = unitVector(v2);
	return acos(dotProd(uv1, uv2));
}



