/*
 * singleProbeChaosLorenz.cpp
 *
 *  Created on: Feb 12, 2023
 *      Author: norms
 */


// LINUX
#include "../headers/singleProbeChaosLorenz.h"
#include "../headers/ProbeAndPoint.h"
#include "../headers/sphereGenerator.h"
#include "../headers/vector3.h"

// WINDOWS
//#include "..\headers\singleProbeChaosLorenz.h"
//#include "..\headers\ProbeAndPoint.h"
//#include "..\headers\sphereGenerator.h"
//#include "..\headers\vector3.h"

#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>

// Main function, probe with sphere of perturbations
double *singleProbeChaosLorenz(Probe probe, double dt, int totalSteps, int ITERATIONS, int NUMPERTURBED, int RADIUS)
{
	// Write the sphere of points to file
	writeSphereToFile(NUMPERTURBED, RADIUS);

	// Initialize perturbed probes
	initializePerturbations(&probe, NUMPERTURBED);

	// Initialize h2 and h3 probes
	probe.h2 = (Probe *)malloc(2 * sizeof(Probe)); copyProbe(&probe.h2[0], probe); copyProbe(&probe.h2[1], probe);
	probe.h3 = (Probe *)malloc(3 * sizeof(Probe)); copyProbe(&probe.h3[0], probe); copyProbe(&probe.h3[1], probe); copyProbe(&probe.h3[2], probe);

	// Calculating h2
	// Initialize two orthonormal vectors:
	probe.h2[0].pos.x += 1.0; probe.h2[1].pos.y += 1.0;

	// Calculating h3
	// Initialize three orthonormal vectors: Memory is already allocated and probe position copied
	probe.h3[0].pos.x += 1.0; probe.h3[1].pos.y += 1.0; probe.h3[2].pos.z += 1.0;


//	std::cout << "Orthonormal vectors initial" << std::endl;
//
//	vector3 uv1 = {probe.perturbed[NUMPERTURBED].pos.x - probe.pos.x, probe.perturbed[NUMPERTURBED].pos.y - probe.pos.y, probe.perturbed[NUMPERTURBED].pos.z - probe.pos.z};
//	vector3 uv2 = {probe.perturbed[NUMPERTURBED + 1].pos.x - probe.pos.x, probe.perturbed[NUMPERTURBED + 1].pos.y - probe.pos.y, probe.perturbed[NUMPERTURBED + 1].pos.z - probe.pos.z};
//	std::cout << "uv1: " << uv1.x << " , " << uv1.y << " , " << uv1.z << std::endl;
//	std::cout << "uv2: " << uv2.x << " , " << uv2.y << " , " << uv2.z << std::endl;
//	std::cout << "Initial angle between vectors: " << angle(uv1, uv2) << std::endl;

	// TODO: Will work on three orthonormal vectors for calculating third lyapunov later:

	// Storage for min and max lyapunov exponents
	double *lyapunovs = (double *)malloc(3 * sizeof(double));

	// Sum of log of relative separation every timestep for each perturbation: ORBITAL SEPARATION METHOD FOR CALCULATING LYAPUNOVS
	// Max lyapunov will be extracted from this at the end
	double *LfSums = (double *)malloc(NUMPERTURBED * sizeof(double)); for (int i = 0; i < NUMPERTURBED; ++i) { LfSums[i] = 0; }
	// h2 lyapunov will be extracted from this at the end
	double areaSum = 0;
	double initialArea = RADIUS * RADIUS;
	// h3 lyapunov will be extracted from this
	double volumeSum = 0;
	double initialVolume = RADIUS * RADIUS * RADIUS;

	// Store data: Store all data of probe and perturbations every ITERATIONS step: Useful for plotting
	int n = ((totalSteps + ITERATIONS - 1) / ITERATIONS) + 1;
	vector3 **data = (vector3 **)malloc(n * sizeof(vector3 *));
	// 0: probe, 1 - NUMPERTURBED + 1: perturbations, NUMPERTURBED + 1 - NUMPERTURBED + 5: h2 and h3
	for (int i = 0; i < n; ++i) { data[i] = (vector3 *)malloc((1 + NUMPERTURBED + 5) * sizeof(vector3)); }

	// Loop
	int c = 0;
	int dc = 0;
	while (c < totalSteps)
	{
		// Store data of Probe and its perturbations every ITERATIONS step
		if (c % ITERATIONS == 0)
		{
			// Probe
			data[dc][0] = {probe.pos.x, probe.pos.y, probe.pos.z};
			// Perturbations
			for (int i = 1; i < NUMPERTURBED + 1; ++i)
			{
				data[dc][i] = {probe.perturbed[i - 1].pos.x, probe.perturbed[i - 1].pos.y, probe.perturbed[i - 1].pos.z};
			}
			// h2
			data[dc][NUMPERTURBED + 1] = {probe.h2[0].pos.x, probe.h2[0].pos.y, probe.h2[0].pos.z}; data[dc][NUMPERTURBED + 2] = {probe.h2[1].pos.x, probe.h2[1].pos.y, probe.h2[1].pos.z};
			// h3
			data[dc][NUMPERTURBED + 3] = {probe.h3[0].pos.x, probe.h3[0].pos.y, probe.h3[0].pos.z}; data[dc][NUMPERTURBED + 4] = {probe.h3[1].pos.x, probe.h3[1].pos.y, probe.h3[1].pos.z}; data[dc][NUMPERTURBED + 5] = {probe.h3[2].pos.x, probe.h3[2].pos.y, probe.h3[2].pos.z};
			dc++;
		}
		// Numerical Integrator: RK4
		probeRK4(&probe, NUMPERTURBED, dt);

		// Function for LfSums, Gram Schmidt2, Gram Schmidt3, and Renormalization of everything
		lyapunovChaosStuffSingle(probe, NUMPERTURBED, RADIUS, LfSums, &areaSum, initialArea, &volumeSum, initialVolume);
		c++;
	}

	// Store latest info
	data[n - 1][0] = {probe.pos.x, probe.pos.y, probe.pos.z};
	for (int i = 1; i < NUMPERTURBED + 3; ++i) { data[n - 1][i] = {probe.perturbed[i - 1].pos.x, probe.perturbed[i - 1].pos.y, probe.perturbed[i - 1].pos.z}; }

	// Calculate and store lyapunov exponents
	extractLyapunovSingle(LfSums, areaSum, volumeSum, NUMPERTURBED, lyapunovs, dt * (double)totalSteps);

	// Write to file
	// perturbationsToFile(data, n, NUMPERTURBED);

	// Free stuff
	for (int i = 0; i < n; ++i) { free(data[i]); } free(data); free(LfSums);

	return lyapunovs;


}

// Write probe and perturbation x,y,z values to file
void perturbationsToFile(vector3 **data, int dataLen, int NUMPERTURBED)
{
	// Desktop
	std::string fpath = "D:\\lorenz3dChaos\\singleProbe\\";
	std::string f1 = "probe.csv"; std::string f2 = "perturbations.csv"; std::string f3 = "gramSchmidtVectors.csv";
	std::ofstream file1; std::ofstream file2; std::ofstream file3;
	file1.open(fpath + f1); file2.open(fpath + f2); file3.open(fpath + f3);

	file1 << "x,y,z" << std::endl;
	file3 << "v1x, v1y, v1z, v2x, v2y, v2z" << std::endl;
	for (int i = 0; i < NUMPERTURBED; ++i)
	{
		file2 << "x" + std::to_string(i) << "," << "y" + std::to_string(i) << "," << "z" + std::to_string(i) << ",";
	}
	file2 << std::endl;
	for (int i = 0; i < dataLen; ++i)
	{
		// Probe
		file1 << std::to_string(data[i][0].x) << "," << std::to_string(data[i][0].y) << "," << std::to_string(data[i][0].z) << std::endl;
		for (int j = 1; j < NUMPERTURBED; ++j)
		{
			// x,y,z for each perturbation
			file2 << std::to_string(data[i][j].x) << "," << std::to_string(data[i][j].y) << "," << std::to_string(data[i][j].z) << ",";
		}
		// Last perturbation
		file2 << std::to_string(data[i][NUMPERTURBED].x) << "," << std::to_string(data[i][NUMPERTURBED].y) << "," << std::to_string(data[i][NUMPERTURBED].z) << std::endl;
		// Two Gram schmidt vectors: x1,y1,z1,x2,y2,z2. Offset by probe's current position for plotting direction
		file3 << std::to_string(data[i][NUMPERTURBED + 1].x - data[i][0].x) << "," << std::to_string(data[i][NUMPERTURBED + 1].y - data[i][0].y) << "," << std::to_string(data[i][NUMPERTURBED + 1].z - data[i][0].z) << "," << std::to_string(data[i][NUMPERTURBED + 2].x - data[i][0].x) << "," << std::to_string(data[i][NUMPERTURBED + 2].y - data[i][0].y) << "," << std::to_string(data[i][NUMPERTURBED + 2].z - data[i][0].z) << std::endl;
	}

	file1.close(); file2.close(); file3.close();
}
// Extracts minimum and maximum lyapunov for single probe with perturbations
void extractLyapunovSingle(double *LfSums, double areaSum, double volumeSum, int NUMPERTURBED, double *lyapunovs, double t)
{
	// Calculate and store max lyapunov
	double max = LfSums[0];
	for (int i = 1; i < NUMPERTURBED; ++i) { if (LfSums[i] > max) { max = LfSums[i]; } }

//	std::cout << std::endl << "Extracting lyapunov results" << std::endl;
//	std::cout << "Total time is: " << t << std::endl;
//	std::cout << "max LfSum is: " << max;
	lyapunovs[0] = max / t;
//	std::cout << "Max lyapunov exponent is: " << max / t << std::endl;
//	std::cout << "areaSum is: " << areaSum << std::endl;
//	std::cout << "RHS is: " << areaSum / t << std::endl;

	// Calculate and store h2
	lyapunovs[1] = (areaSum / t) - lyapunovs[0];
	std::cout << "h1 is: " << std::setprecision(5) << lyapunovs[0] << std::endl;
	std::cout << "areaSum / t is: " << std::setprecision(5) << areaSum / t << std::endl;

	// Calculate and store h3
	lyapunovs[2] = (volumeSum / t);
}

// Handles the timestep processes for calculating min and max lyapunov exponents: Might be broken
void lyapunovChaosStuffSingle(Probe probe, int NUMPERTURBED, double RADIUS, double *LfSums, double *areaSum, double initialArea, double *volumeSum, double initialVolume)
{
	// Max lyapunov
	double maxMag = mag((vector3){probe.perturbed[0].pos.x - probe.pos.x, probe.perturbed[0].pos.y - probe.pos.y, probe.perturbed[0].pos.z - probe.pos.z});
	for (int i = 0; i < NUMPERTURBED; ++i)
	{
		// Calculate Lf of each perturbation, and renormalize the position vector
		vector3 v = {probe.perturbed[i].pos.x - probe.pos.x, probe.perturbed[i].pos.y - probe.pos.y, probe.perturbed[i].pos.z - probe.pos.z};
		// RADIUS = 1.0
		double dist = mag(v);
		double logDist = log(dist);

		// DEBUGGING
		if (dist > maxMag) { maxMag = dist; }

		LfSums[i] += logDist;
		vector3 uv = unitVector(v);
		probe.perturbed[i].pos = {probe.pos.x + uv.x, probe.pos.y + uv.y, probe.pos.z + uv.z};
		// RADIUS > 1.0
//		 LfSums[i] += log(mag(v) / RADIUS);
//		 probe.perturbed[i].pos = {(uv.x + probe.pos.x) * RADIUS, (uv.y + probe.pos.y) * RADIUS, (uv.z + probe.pos.z) * RADIUS};
	}

	// DEBUGGING MAX LYAPUNOV
//	std::cout << std::endl << "For max lyapunov" << std::endl;
//	std::cout << "Distance is: " << maxMag << std::endl;
//	std::cout << "log distance is: " << log(maxMag) << std::endl;

	// h2 Lyapunov
	vector3 v1 = {probe.h2[0].pos.x - probe.pos.x, probe.h2[0].pos.y - probe.pos.y, probe.h2[0].pos.z - probe.pos.z};
	vector3 v2 = {probe.h2[1].pos.x - probe.pos.x, probe.h2[1].pos.y - probe.pos.y, probe.h2[1].pos.z - probe.pos.z};
	vector3 v1Crossv2 = crossProd(v1, v2);
	// RADIUS = 1
	double logArea = log(mag(v1Crossv2));
	// DEBUGGING

//	std::cout << "For min lyapunov" << std::endl;
//	std::cout << "Areasum is " << *areaSum << std::endl;

	(*areaSum) += logArea;
	// Debugging

//	std::cout << "Area of parallelogram is: " << magv1 * magv2 * sin(a) << std::endl;
//	std::cout << "log area of parallelogram is: " << logArea << std::endl << std::endl;
//	std::cout << "Areasum is now " << *areaSum << std::endl;
	// RADIUS > 1
	// areaSums += log(magv1 * magv2 * sin(a) / RADIUS);

	// Reorthonormalize h2 vectors
	gramSchmidtRenormalizationh2Single(&probe, NUMPERTURBED, v1, v2);

	//h3 lyapunov and renormalization: Assuming initial volume of 1
	gramSchmidtRenormalizationh3Single(&probe, NUMPERTURBED, volumeSum, logArea);
}

// Gram schmidt Renormalization for single probe: Renormalizes h3 vectors into three orthogonal unit vectors that preserve direction
void gramSchmidtRenormalizationh3Single(Probe *probe, int NUMPERTURBED, double *volumeSum, double logArea)
{
	// Three vectors
	vector3 v1 = {probe->h3[0].pos.x - probe->pos.x, probe->h3[0].pos.y - probe->pos.y, probe->h3[0].pos.z - probe->pos.z};
	vector3 v2 = {probe->h3[1].pos.x - probe->pos.x, probe->h3[1].pos.y - probe->pos.y, probe->h3[1].pos.z - probe->pos.z};
	vector3 v3 = {probe->h3[2].pos.x - probe->pos.x, probe->h3[2].pos.y - probe->pos.y, probe->h3[2].pos.z - probe->pos.z};

	// Add volume of parallelepiped to volume sum: V = |(v1 X v2) dot v3|
	vector3 v1Crossv2 = crossProd(v1, v2);
	double volume = abs(dotProd(v1Crossv2, v3));
	(*volumeSum) += (log(volume) - logArea);


	// Gram schmidt: Calculate orthonormal vectors that preserve direction
	// Step to get e2
	vector3 e1 = unitVector(v1);
	double v2Dote1 = dotProd(v2, e1);
	vector3 y2 = {v2.x - (v2Dote1) * e1.x, v2.y - (v2Dote1) * e1.y, v2.z - (v2Dote1) * e1.z};
	vector3 e2 = unitVector(y2);

	// Step to get e3
	double v3Dote1 = dotProd(v3, e1); vector3 j1 = {v3Dote1 * e1.x, v3Dote1 * e1.y, v3Dote1 * e1.z};
	double v3Dote2 = dotProd(v3, e2); vector3 j2 = {v3Dote2 * e2.x, v3Dote2 * e2.y, v3Dote2 * e2.z};

	vector3 j1Plusj2 = {j1.x + j2.x, j1.y + j2.y, j1.z + j2.z};
	vector3 y3 = {v3.x - j1Plusj2.x, v3.y - j1Plusj2.y, v3.z - j1Plusj2.z};
	vector3 e3 = unitVector(y3);

//	std::cout << std::endl << std::endl << "GRAM SCHMIDT H3" << std::endl;
//	std::cout << "Original vectors:" << std::endl;
//	std::cout << "v1: " << v1.x << " " << v1.y << " " << v1.z << std::endl;
//	std::cout << "v2: " << v2.x << " " << v2.y << " " << v2.z << std::endl;
//	std::cout << "v3: " << v3.x << " " << v3.y << " " << v3.z << std::endl;
//	std::cout << "Orthonormal vectors" << std::endl;
//	std::cout << "e1: " << e1.x << " " << e1.y << " " << e1.z << std::endl;
//	std::cout << "e2: " << e2.x << " " << e2.y << " " << e2.z << std::endl;
//	std::cout << "e3: " << e3.x << " " << e3.y << " " << e3.z << std::endl;

	probe->h3[0].pos = {probe->pos.x + e1.x, probe->pos.y + e1.y, probe->pos.z + e1.z};
	probe->h3[1].pos = {probe->pos.x + e2.x, probe->pos.y + e2.y, probe->pos.z + e2.z};
	probe->h3[2].pos = {probe->pos.x + e3.x, probe->pos.y + e3.y, probe->pos.z + e3.z};
}
// Gram schmidt Renormalization for single probe: Renormalizes h2 vectors into two orthogonal unit vectors that preserve direction
void gramSchmidtRenormalizationh2Single(Probe *probe, int NUMPERTURBED, vector3 v1, vector3 v2)
{
	// u1 = v1 / magv1: y2 = v2 - (v2 dot u1)*u1: u2 = y2 / magy2
	double magv1 = mag(v1); double magv2 = mag(v2);
	vector3 u1 = {v1.x / magv1, v1.y / magv1, v1.z / magv1};
	double v2Dotu1 = dotProd(v2, u1);
	vector3 y2 = {v2.x - (v2Dotu1) * u1.x, v2.y - (v2Dotu1) * u1.y, v2.z - (v2Dotu1) * u1.z};
	double magy2 = mag(y2);
	vector3 u2 = {y2.x / magy2, y2.y / magy2, y2.z / magy2};
	// u1 and u2 are new orthogonal unit vectors to assign to v1 and v2 of probe
//	std::cout << std::endl << "GRAM SCHMIDT ORTHOGONALIZATION" << std::endl;
//	std::cout << "Original vectors" << std::endl;
//	std::cout << "v1: " << v1.x << " , " << v1.y << " , " << v1.z << std::endl;
//	std::cout << "v2: " << v2.x << " , " << v2.y << " , " << v2.z << std::endl;
//	std::cout << "Angle between vectors: " << angle(v1, v2) * (180.0) / M_PI << std::endl;
//	std::cout << "Calculated Orthonormal vectors" << std::endl;
//	std::cout << "u1: " << u1.x << " , " << u1.y << " , " << u1.z << std::endl;
//	std::cout << "u2: " << u2.x << " , " << u2.y << " , " << u2.z << std::endl << std::endl;
//	std::cout << "Angle between orthonormal vectors: " << angle(u1, u2) * (180.0) / M_PI << std::endl;
	probe->h2[0].pos = {probe->pos.x + u1.x, probe->pos.y + u1.y, probe->pos.z + u1.z};
	probe->h2[1].pos = {probe->pos.x + u2.x, probe->pos.y + u2.y, probe->pos.z + u2.z};
}


// RK4 for probe and its perturbations in a lorenz system
void probeRK4(Probe *probe, int NUMPERTURBED, double dt)
{
	// Initialize Object arrays
	Probe k2Probe, k3Probe, k4Probe;

	// Initialize k's: Index 0 is for probe, last five are for gram schmidt vectors
	vector3 *k1 = (vector3 *)malloc((6 + NUMPERTURBED) * sizeof(vector3)); vector3 *k2 = (vector3 *)malloc((6 + NUMPERTURBED) * sizeof(vector3)); vector3 *k3 = (vector3 *)malloc((6 + NUMPERTURBED) * sizeof(vector3)); vector3 *k4 = (vector3 *)malloc((6 + NUMPERTURBED) * sizeof(vector3)); vector3 *weightedK = (vector3 *)malloc((6 + NUMPERTURBED) * sizeof(vector3));

	// Copy data into probeCopies, perturbations, h2, and h3 for all k
	copyProbe(&k2Probe, *probe); copyProbe(&k3Probe, *probe); copyProbe(&k4Probe, *probe);
	k2Probe.perturbed = (Probe *)malloc((NUMPERTURBED) * sizeof(Probe)); k3Probe.perturbed = (Probe *)malloc((NUMPERTURBED) * sizeof(Probe)); k4Probe.perturbed = (Probe *)malloc((NUMPERTURBED) * sizeof(Probe));
	k2Probe.h2 = (Probe *)malloc(2 * sizeof(Probe)); k3Probe.h2 = (Probe *)malloc(2 * sizeof(Probe)); k4Probe.h2 = (Probe *)malloc(2 * sizeof(Probe));
	k2Probe.h3 = (Probe *)malloc(3 * sizeof(Probe)); k3Probe.h3 = (Probe *)malloc(3 * sizeof(Probe)); k4Probe.h3 = (Probe *)malloc(3 * sizeof(Probe));
	for (int i = 0; i < NUMPERTURBED; ++i) { copyProbe(&k2Probe.perturbed[i], probe->perturbed[i]); copyProbe(&k3Probe.perturbed[i], probe->perturbed[i]); copyProbe(&k4Probe.perturbed[i], probe->perturbed[i]); }
	for (int i = 0; i < 2; ++i) { copyProbe(&k2Probe.h2[i], probe->h2[i]); copyProbe(&k3Probe.h2[i], probe->h2[i]); copyProbe(&k4Probe.h2[i], probe->h2[i]); }
	for (int i = 0; i < 3; ++i) { copyProbe(&k2Probe.h3[i], probe->h3[i]); copyProbe(&k3Probe.h3[i], probe->h3[i]); copyProbe(&k4Probe.h3[i], probe->h3[i]); }

	// k1 step
	// Extract k1 and put in arrays
	extractLorenzKProbe(*probe, NUMPERTURBED, k1);

	// k2 Step
	// Move k2 object copies dt / 2 using k1
	probeRK4Step(&k2Probe, NUMPERTURBED, k1, dt / 2.0);
	extractLorenzKProbe(k2Probe, NUMPERTURBED, k2);

	// k3 Step
	// Move k3 object copies dt / 2 using k2 Accelerations and k2 Velocities
	probeRK4Step(&k3Probe, NUMPERTURBED, k2, dt / 2.0);
	extractLorenzKProbe(k3Probe, NUMPERTURBED, k3);

	// k4 Step
	// Move k4 object copies dt using k3 Accelerations and k3 Velocities
	probeRK4Step(&k4Probe, NUMPERTURBED, k3, dt);
	extractLorenzKProbe(k4Probe, NUMPERTURBED, k4);

	//rk4 step
	// Fill up weighted k's
	for (int i = 0; i < NUMPERTURBED + 6; ++i)
	{
		weightedK[i] = {k1[i].x + 2 * k2[i].x + 2 * k3[i].x + k4[i].x, k1[i].y + 2 * k2[i].y + 2 * k3[i].y + k4[i].y, k1[i].z + 2 * k2[i].z + 2 * k3[i].z + k4[i].z};
	}
	probeRK4Step(probe, NUMPERTURBED, weightedK, dt / 6.0);

	// Free stuff
	free(k1); free(k2); free(k3); free(k4); free(weightedK);
}

// Single step for RK4: single probecopy with perturbations and two orthogonals: h2, and three orthogonals: h3
void probeRK4Step(Probe *probeCopy, int NUMPERTURBED, vector3 *k, double dt)
{
	// probe
	probeCopy->pos.x += k[0].x * dt; probeCopy->pos.y += k[0].y * dt; probeCopy->pos.z += k[0].z * dt;
	// Perturbations
	for (int i = 1; i < NUMPERTURBED + 1; ++i)
	{
		probeCopy->perturbed[i - 1].pos.x += k[i].x * dt; probeCopy->perturbed[i - 1].pos.y += k[i].y * dt; probeCopy->perturbed[i - 1].pos.z += k[i].z * dt;
	}
	// h2
	for (int i = 0; i < 2; ++i)
	{
		probeCopy->h2[i].pos.x += k[NUMPERTURBED + 1 + i].x * dt; probeCopy->h2[i].pos.y += k[NUMPERTURBED + 1 + i].y * dt; probeCopy->h2[i].pos.z += k[NUMPERTURBED + 1 + i].z * dt;
	}
	// h3
	for (int i = 0; i < 3; ++i)
	{
		probeCopy->h3[i].pos.x += k[NUMPERTURBED + 3 + i].x * dt; probeCopy->h3[i].pos.y += k[NUMPERTURBED + 3 + i].y * dt; probeCopy->h3[i].pos.z += k[NUMPERTURBED + 3 + i].z * dt;
	}
}
// Stores lorenz change on x, y, z using lorenz equations: Uses Lorenz attract initial conditions for single probe and its perturbations: h2 and h3 probes too
void extractLorenzKProbe(Probe probe, int NUMPERTURBED, vector3 *k)
{
	// Defined to use Lorenz Attractor constants
	double ROW = 28.0;
	double SIGMA = 10.0;
	double BETA = 8.0 / 3.0;

	// Probe
	k[0] = {SIGMA * (probe.pos.y - probe.pos.x), probe.pos.x * (ROW - probe.pos.z) - probe.pos.y, (probe.pos.x * probe.pos.y) - (BETA * probe.pos.z)};
	// Perturbations
	for (int i = 1; i < NUMPERTURBED + 1; ++i)
	{
		k[i] = {SIGMA * (probe.perturbed[i - 1].pos.y - probe.perturbed[i - 1].pos.x), probe.perturbed[i - 1].pos.x * (ROW - probe.perturbed[i - 1].pos.z) - probe.perturbed[i - 1].pos.y, (probe.perturbed[i - 1].pos.x * probe.perturbed[i - 1].pos.y) - (BETA * probe.perturbed[i - 1].pos.z)};
	}
	// h2
	for (int i = 0; i < 2; ++i)
	{
		k[NUMPERTURBED + 1 + i] = {SIGMA * (probe.h2[i].pos.y - probe.h2[i].pos.x), probe.h2[i].pos.x * (ROW - probe.h2[i].pos.z) - probe.h2[i].pos.y, (probe.h2[i].pos.x * probe.h2[i].pos.y) - (BETA * probe.h2[i].pos.z)};
	}
	// h3
	for (int i = 0; i < 3; ++i)
	{
		k[NUMPERTURBED + 3 + i] = {SIGMA * (probe.h3[i].pos.y - probe.h3[i].pos.x), probe.h3[i].pos.x * (ROW - probe.h3[i].pos.z) - probe.h3[i].pos.y, (probe.h3[i].pos.x * probe.h3[i].pos.y) - (BETA * probe.h3[i].pos.z)};
	}
}

// Initializes the sphere of probes around the probe object, which will then progress into an ellipsoid. Reads from file
void initializePerturbations(Probe *p, int NUMPOINTS)
{
	// Initialize all perturbation probes: Last five: 0,1 are for h2 lyapunov: 2,3,4 are for h3 lyapunov
	p->perturbed = (Probe *)malloc(NUMPOINTS * sizeof(Probe));
	for (int i = 0; i < NUMPOINTS; ++i) { copyProbe(&p->perturbed[i], *p); }

	// Read from file
	// Desktop path
	FILE *f = fopen("D:\\lorenz3dChaos\\unitSphere.csv", "r");
	fseek(f, 0L, SEEK_END);
	long int fsize = ftell(f);
	rewind(f);
	// Read data into buffer
	char *buffer = (char *)malloc(fsize * sizeof(char));
	fread(buffer, fsize, 1, f);
	fclose(f);
	// Read from buffer into perturbed probes
	int pointCount = 0;
	char arr[30];
	int arrCount = 0;
	int commaCount = 0;
	int i = 0;

	while (pointCount < NUMPOINTS)
	{
		if (buffer[i] == '\n')
		{
			arr[arrCount] = '\0';
			p->perturbed[pointCount].pos.z = std::atof(arr) + p->pos.z;
			pointCount++;
			arrCount = 0;
			commaCount = 0;

		}
		else if (buffer[i] == ',')
		{
			arr[arrCount] = '\0';
			if (commaCount == 0)
			{
				p->perturbed[pointCount].pos.x = std::atof(arr) + p->pos.x;
			}
			else
			{
				p->perturbed[pointCount].pos.y = std::atof(arr) + p->pos.y;
			}
			arrCount = 0;
			commaCount++;
		}
		else
		{
			arr[arrCount] = buffer[i];
			arrCount++;
		}
		++i;
	}
	free(buffer);
}
