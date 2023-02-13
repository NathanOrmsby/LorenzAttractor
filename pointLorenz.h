/*
 * pointLorenz.h
 *
 *  Created on: Feb 12, 2023
 *      Author: norms
 */

#ifndef HEADERS_POINTLORENZ_H_
#define HEADERS_POINTLORENZ_H_

#include "ProbeAndPoint.h"

// Main function pointLorenz: Outputs positions over time to file for a number of points through the Lorenz attractor
void pointLorenz(Point *points, int pointNum, double dt, int totalSteps, int ITERATIONS);
// RK4 for points in lorenz
void pointsRK4(Point *points, int pointNum, double dt);
// Stores lorenz change on x, y, z using lorenz equations: Uses Lorenz attract initial conditions
void extractLorenzK(Point *points, int pointNum, vector3 *k);
// Moves Points one step forward. Helper function to pointsRK4
void pointsRK4Step(Point *pointCopies, int pointNum, vector3 *k, double dt);
// Write paths of points to a different file per point
void pathsToFile(vector3 **data, int n, int pointNum);

#endif /* HEADERS_POINTLORENZ_H_ */
