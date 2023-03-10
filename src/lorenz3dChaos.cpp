//============================================================================
// Name        : lorenz3dChaos.cpp
// Author      : Nathan Ormsby
// Version     :
// Copyright   : DO NOT COPY MY CODE, it probably wont work
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <math.h>

#include "headers\ProbeAndPoint.h"
#include "headers\sphereGenerator.h"
#include "headers\vector3.h"
#include "headers\pointLorenz.h"
#include "headers\singleProbeChaosLorenz.h"

int main()
{
	// Timestep
	double dt = 0.01;
	int totalSteps = 310000;

	// Data storage
	int ITERATIONS = 280000;

	// Probe Lorenz Main function
	// Initialize Probe
	Probe probe;
	probe.pos = {0.01, 0, 0};
	int NUMPERTURBED = 20; 					// Number of perturbations
	double *lyapunovs = singleProbeChaosLorenz(probe, dt, totalSteps, ITERATIONS, NUMPERTURBED, 1.0);
	std::cout << "h1: " << lyapunovs[0] << " h2: " << lyapunovs[1] << " h3: " << lyapunovs[2] << std::endl;
	free(lyapunovs);
	return 0;

	// Point Lorenz Main function
	// Initialize points
//	int pointNum = 3;
//	Point *points = (Point *)malloc(pointNum * sizeof(Point));
//	points[0].pos = {0.01, 0, 0};
//	points[1].pos = {1.0, 0, 0};
//	points[2].pos = {-1.0, 0, 0};
//
//	// Run the stuff
//	pointLorenz(points, pointNum, dt, totalSteps, ITERATIONS);
//	std::cout << "done" << std::endl;
	return 0;

}
