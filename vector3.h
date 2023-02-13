/*
 * vector3.h
 *
 *  Created on: Feb 12, 2023
 *      Author: norms
 */

#ifndef HEADERS_VECTOR3_H_
#define HEADERS_VECTOR3_H_


typedef struct
{
	double x;
	double y;
	double z;
} vector3;

double dotProd(vector3 v1, vector3 v2);
double mag(vector3 v);
// Returns unit vector of a vector3
vector3 unitVector(vector3 v);
// Returns angle between two vector3s
double angle(vector3 v1, vector3 v2);
// Returns angle between unit vectors of v1 and v2, good if worried about overflow
double unitAngle(vector3 v1, vector3 v2);



#endif /* HEADERS_VECTOR3_H_ */
