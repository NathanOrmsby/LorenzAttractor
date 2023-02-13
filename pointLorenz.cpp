/*
 * planetLorenz.cpp
 *
 *  Created on: Feb 12, 2023
 *      Author: norms
 */

#include "..\headers\pointLorenz.h"
#include "..\headers\ProbeAndPoint.h"
#include "..\headers\sphereGenerator.h"
#include "..\headers\vector3.h"

#include <iostream>
#include <fstream>

// Main function pointLorenz: Outputs positions over time to file for a number of points through the Lorenz attractor
void pointLorenz(Point *points, int pointNum, double dt, int totalSteps, int ITERATIONS)
{
	// Loop
	int c = 0, dc = 0;
	// Store data: Store all x,y,z positions of every ITERATIONS step for every Point object
	int n = ((totalSteps + ITERATIONS - 1) / ITERATIONS) + 1;
	vector3 **data = (vector3 **)malloc(n * sizeof(vector3 *));
	for (int i = 0; i < n; ++i) { data[i] = (vector3 *)malloc(pointNum * sizeof(vector3)); }

	while (c < totalSteps)
	{
		// Store data of JW and its perturbations every ITERATIONS step
		if (c % ITERATIONS == 0)
		{
			for (int i = 0; i < pointNum; ++i)
			{
				data[dc][i] = {points[i].pos.x, points[i].pos.y, points[i].pos.z};
//				std::cout << "Point " << i << " Pos: x: " << data[dc][3 * i] << " y: " << data[dc][3 * i + 1] << " z: " << data[dc][3 * i + 2] << std::endl;
			}
			dc++;
		}
		// Numerical Integration: rk4
		pointsRK4(points, pointNum, dt);
		c++;
	}

	// Write data to file
	pathsToFile(data, n, pointNum);

	// Free stuff
	for (int i = 0; i < n; ++i) { free(data[i]); } free(data);
}

// RK4 for points in lorenz
void pointsRK4(Point *points, int pointNum, double dt)
{
	// Initialize Object arrays
	Point k2Points[pointNum], k3Points[pointNum], k4Points[pointNum];
	// Initialize k's
	vector3 *k1 = (vector3 *)malloc(pointNum * sizeof(vector3)); vector3 *k2 = (vector3 *)malloc(pointNum * sizeof(vector3)); vector3 *k3 = (vector3 *)malloc(pointNum * sizeof(vector3)); vector3 *k4 = (vector3 *)malloc(pointNum * sizeof(vector3)); vector3 *weightedK = (vector3 *)malloc(pointNum * sizeof(vector3));
	// Copy data into Point arrays for all k
	for (int i = 0; i < pointNum; ++i)
	{
		copyPoint(&k2Points[i], points[i]); copyPoint(&k3Points[i], points[i]);  copyPoint(&k4Points[i], points[i]);
	}

	// k1 step
	// Extract k1 and put in arrays
	extractLorenzK(points, pointNum, k1);

	// k2 Step
	// Move k2 object copies dt / 2 using k1
	pointsRK4Step(k2Points, pointNum, k1, dt / 2.0);
	extractLorenzK(k2Points, pointNum, k2);

	// k3 Step
	// Move k3 object copies dt / 2 using k2 Accelerations and k2 Velocities
	pointsRK4Step(k3Points, pointNum, k2, dt / 2.0);
	extractLorenzK(k3Points, pointNum, k3);

	// k4 Step
	// Move k4 object copies dt using k3 Accelerations and k3 Velocities
	pointsRK4Step(k4Points, pointNum, k3, dt);
	extractLorenzK(k4Points, pointNum, k4);

	// rk4 step
	for (int i = 0; i < pointNum; ++i)
	{
		weightedK[i] = {k1[i].x + 2 * k2[i].x + 2 * k3[i].x + k4[i].x, k1[i].y + 2 * k2[i].y + 2 * k3[i].y + k4[i].y, k1[i].z + 2 * k2[i].z + 2 * k3[i].z + k4[i].z};
	}
	// Move points
	pointsRK4Step(points, pointNum, weightedK, dt / 6.0);

	// Free stuff
	free(k1); free(k2); free(k3); free(k4); free(weightedK);
}

// Moves Points one step forward. Helper function to pointsRK4
void pointsRK4Step(Point *pointCopies, int pointNum, vector3 *k, double dt)
{
	for (int i = 0; i < pointNum; ++i)
	{
		pointCopies[i].pos.x += k[i].x * dt; pointCopies[i].pos.y += k[i].y * dt; pointCopies[i].pos.z += k[i].z * dt;
	}
}

// Stores lorenz change on x, y, z using lorenz equations: Uses Lorenz attract initial conditions
void extractLorenzK(Point *points, int pointNum, vector3 *k)
{
	// Defined to use Lorenz Attractor constants
	double ROW = 28.0;
	double SIGMA = 10.0;
	double BETA = 8.0 / 3.0;
	for (int i = 0; i < pointNum; ++i)
	{
		k[i] = {SIGMA * (points[i].pos.y - points[i].pos.x), points[i].pos.x * (ROW - points[i].pos.z) - points[i].pos.y, (points[i].pos.x * points[i].pos.y) - (BETA * points[i].pos.z)};
	}
}

// Write paths of points to a different file per point
void pathsToFile(vector3 **data, int n, int pointNum)
{
	std::string fpath = "D:\\lorenz3dChaos\\pointPaths\\";

	// Write a file for each point
	for (int i = 0; i < pointNum; ++i)
	{
		std::ofstream file;
		std::string point = "point" + std::to_string(i) + ".csv";
		file.open(fpath + point);
		file << "x,y,z" << std::endl;
		for (int j = 0; j < n; ++j)
		{
			file << std::to_string(data[j][i].x) << "," << std::to_string(data[j][i].y) << "," << std::to_string(data[j][i].z) << std::endl;
		}
	}
}

