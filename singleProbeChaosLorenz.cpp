/*
 * singleProbeChaosLorenz.cpp
 *
 *  Created on: Feb 12, 2023
 *      Author: norms
 */


#include "..\headers\singleProbeChaosLorenz.h"
#include "..\headers\ProbeAndPoint.h"
#include "..\headers\sphereGenerator.h"
#include "..\headers\vector3.h"

#include <math.h>
#include <iostream>
#include <fstream>

// Main function, probe with sphere of perturbations
double *singleProbeChaosLorenz(Probe probe, double dt, int totalSteps, int ITERATIONS, int NUMPERTURBED, int RADIUS)
{
	// Write the sphere of points to file
	writeSphereToFile(NUMPERTURBED, RADIUS);

	// Initialize perturbed probes
	initializePerturbations(&probe, NUMPERTURBED);


	// Initialize two orthonormal vectors: Memory is already allocated and probe position copied
	probe.perturbed[NUMPERTURBED].pos.x += 1.0; probe.perturbed[NUMPERTURBED + 1].pos.y += 1.0;

	// TODO: Will work on three orthonormal vectors for calculating third lyapunov later:

	// Storage for min and max lyapunov exponents
	double *lyapunovs = (double *)malloc(2 * sizeof(double));

	// Sum of log of relative separation every timestep for each perturbation: ORBITAL SEPARATION METHOD FOR CALCULATING LYAPUNOVS
	// Max lyapunov will be extracted from this at the end
	double *LfSums = (double *)malloc(NUMPERTURBED * sizeof(double)); for (int i = 0; i < NUMPERTURBED; ++i) { LfSums[i] = 0; }
	// Min lyapunov will be extracted from this at the end
	double areaSum = 0;
	double initialArea = RADIUS * RADIUS;

	// Store data: Store all data of probe and perturbations every ITERATIONS step: Useful for plotting
	int n = ((totalSteps + ITERATIONS - 1) / ITERATIONS) + 1;
	vector3 **data = (vector3 **)malloc(n * sizeof(vector3 *));
	for (int i = 0; i < n; ++i) { data[i] = (vector3 *)malloc((3 + NUMPERTURBED) * sizeof(vector3)); }

	// Loop
	int c = 0;
	int dc = 0;
	while (c < totalSteps)
	{
		// Store data of Probe and its perturbations every ITERATIONS step
		if (c % ITERATIONS == 0)
		{
			data[dc][0] = {probe.pos.x, probe.pos.y, probe.pos.z};
			for (int i = 1; i < NUMPERTURBED + 3; ++i)
			{
				data[dc][i] = {probe.perturbed[i - 1].pos.x, probe.perturbed[i - 1].pos.y, probe.perturbed[i - 1].pos.z};
			}
			dc++;
		}
		// Numerical Integrator: RK4
		probeRK4(&probe, NUMPERTURBED, dt);

		// Function for LfSums, Gram Schmidt, and Renormalization
		lyapunovChaosStuffSingle(probe, NUMPERTURBED, RADIUS, LfSums, areaSum, initialArea);
		c++;
	}

	// Store latest info
	data[n - 1][0] = {probe.pos.x, probe.pos.y, probe.pos.z};
	for (int i = 1; i < NUMPERTURBED + 3; ++i) { data[n - 1][i] = {probe.perturbed[i - 1].pos.x, probe.perturbed[i - 1].pos.y, probe.perturbed[i - 1].pos.z}; }

	// Calculate and store lyapunov exponents
	extractLyapunovSingle(LfSums, areaSum, NUMPERTURBED, lyapunovs, dt * (double)totalSteps);

	// Write to file
	perturbationsToFile(data, n, NUMPERTURBED);

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
void extractLyapunovSingle(double *LfSums, double areaSum, int NUMPERTURBED, double *lyapunovs, double t)
{
	// Calculate and store max lyapunov
	double max = LfSums[0];
	for (int i = 1; i < NUMPERTURBED; ++i) { if (LfSums[i] > max) { max = LfSums[i]; } }
	lyapunovs[1] = max / t;

	// Calculate and store min lyapunov
	lyapunovs[0] = (areaSum / t) - lyapunovs[1];
}

// Handles the timestep processes for calculating min and max lyapunov exponents: Might be broken
void lyapunovChaosStuffSingle(Probe probe, int NUMPERTURBED, double RADIUS, double *LfSums, double areaSum, double initialArea)
{
	// Max lyapunov
	for (int i = 0; i < NUMPERTURBED; ++i)
	{
		// Calculate Lf of each perturbation, and renormalize the position vector
		vector3 v = {probe.perturbed[i].pos.x - probe.pos.x, probe.perturbed[i].pos.y - probe.pos.y, probe.perturbed[i].pos.z - probe.pos.z};
		// RADIUS = 1.0
		LfSums[i] += log(mag(v));
		vector3 uv = unitVector(v);
		probe.perturbed[i].pos = {uv.x + probe.pos.x, uv.y + probe.pos.y, uv.z + probe.pos.z};
		// RADIUS > 1.0
		// LfSums[i] += log(mag(v) / RADIUS);
		// probe.perturbed[i].pos = {(uv.x + probe.pos.x) * RADIUS, (uv.y + probe.pos.y) * RADIUS, (uv.z + probe.pos.z) * RADIUS};
	}
	// Min lyapunov
	vector3 v1 = {probe.perturbed[NUMPERTURBED].pos.x - probe.pos.x, probe.perturbed[NUMPERTURBED].pos.y - probe.pos.y, probe.perturbed[NUMPERTURBED].pos.z - probe.pos.z};
	vector3 v2 = {probe.perturbed[NUMPERTURBED + 1].pos.x - probe.pos.x, probe.perturbed[NUMPERTURBED + 1].pos.y - probe.pos.y, probe.perturbed[NUMPERTURBED + 1].pos.z - probe.pos.z};
	double magv1 = mag(v1); double magv2 = mag(v2); double a = angle(v1, v2);
	// RADIUS = 1
	areaSum += log(magv1 * magv2 * sin(a));
	// RADIUS > 1
	// areaSums += log(magv1 * magv2 * sin(a) / RADIUS);

	// Reorthonormalize vectors
//	gramSchmidtRenormalizationSingle(&probe, NUMPERTURBED, v1, v2, magv1, magv2);
}

// Gram schmidt Renormalization for single probe: Renormalizes vectors into two orthogonal unit vectors that preserve direction
void gramSchmidtRenormalizationSingle(Probe *probe, int NUMPERTURBED, vector3 v1, vector3 v2, double magv1, double magv2)
{
	// u1 = v1 / magv1: y2 = v2 - (v2 dot u1)*u1: u2 = y2 / magy2
	vector3 u1 = {v1.x / magv1, v1.y / magv1, v1.z / magv1};
	double v2Dotu1 = dotProd(v2, u1);
	vector3 y2 = {v2.x - (v2Dotu1) * u1.x, v2.y - (v2Dotu1) * u1.y, v2.z - (v2Dotu1) * u1.z};
	double magy2 = mag(y2);
	// u1 and u2 are new orthogonal unit vectors to assign to v1 and v2 of probe
	probe->perturbed[NUMPERTURBED].pos = {probe->pos.x + u1.x, probe->pos.y + u1.y, probe->pos.z + u1.z};
	probe->perturbed[NUMPERTURBED + 1].pos = {probe->pos.x + (y2.x / magy2), probe->pos.y + (y2.y / magy2), probe->pos.z + (y2.z / magy2)};
}


// RK4 for probe and its perturbations in a lorenz system
void probeRK4(Probe *probe, int NUMPERTURBED, double dt)
{
	// Initialize Object arrays
	Probe k2Probe, k3Probe, k4Probe;
	// Initialize k's: Index 0 is for probe, last two are for gram schmidt vectors
	vector3 *k1 = (vector3 *)malloc((3 + NUMPERTURBED) * sizeof(vector3)); vector3 *k2 = (vector3 *)malloc((3 + NUMPERTURBED) * sizeof(vector3)); vector3 *k3 = (vector3 *)malloc((3 + NUMPERTURBED) * sizeof(vector3)); vector3 *k4 = (vector3 *)malloc((3 + NUMPERTURBED) * sizeof(vector3)); vector3 *weightedK = (vector3 *)malloc((3 + NUMPERTURBED) * sizeof(vector3));
	// Copy data into probeCopies and their perturbations for all k
	copyProbe(&k2Probe, *probe); copyProbe(&k3Probe, *probe); copyProbe(&k4Probe, *probe); k2Probe.perturbed = (Probe *)malloc((NUMPERTURBED + 2) * sizeof(Probe)); k3Probe.perturbed = (Probe *)malloc((NUMPERTURBED + 2) * sizeof(Probe)); k4Probe.perturbed = (Probe *)malloc((NUMPERTURBED + 2) * sizeof(Probe));
	for (int i = 0; i < NUMPERTURBED + 2; ++i) { copyProbe(&k2Probe.perturbed[i], probe->perturbed[i]); copyProbe(&k3Probe.perturbed[i], probe->perturbed[i]); copyProbe(&k4Probe.perturbed[i], probe->perturbed[i]); }

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
	for (int i = 0; i < NUMPERTURBED + 3; ++i)
	{
		weightedK[i] = {k1[i].x + 2 * k2[i].x + 2 * k3[i].x + k4[i].x, k1[i].y + 2 * k2[i].y + 2 * k3[i].y + k4[i].y, k1[i].z + 2 * k2[i].z + 2 * k3[i].z + k4[i].z};
	}
	probeRK4Step(probe, NUMPERTURBED, weightedK, dt / 6.0);

	// Free stuff
	free(k1); free(k2); free(k3); free(k4); free(weightedK);
}

// Single step for RK4: single probecopy with perturbations and two orthogonals
void probeRK4Step(Probe *probeCopy, int NUMPERTURBED, vector3 *k, double dt)
{
	// probe
	probeCopy->pos.x += k[0].x * dt; probeCopy->pos.y += k[0].y * dt; probeCopy->pos.z += k[0].z * dt;
	// Perturbations and gram schmidts
	for (int i = 1; i < NUMPERTURBED + 3; ++i)
	{
		probeCopy->perturbed[i - 1].pos.x += k[i].x * dt; probeCopy->perturbed[i - 1].pos.y += k[i].y * dt; probeCopy->perturbed[i - 1].pos.z += k[i].z * dt;
	}
}
// Stores lorenz change on x, y, z using lorenz equations: Uses Lorenz attract initial conditions for single probe and its perturbations
void extractLorenzKProbe(Probe probe, int NUMPERTURBED, vector3 *k)
{
	// Defined to use Lorenz Attractor constants
	double ROW = 28.0;
	double SIGMA = 10.0;
	double BETA = 8.0 / 3.0;

	// Probe
	k[0] = {SIGMA * (probe.pos.y - probe.pos.x), probe.pos.x * (ROW - probe.pos.z) - probe.pos.y, (probe.pos.x * probe.pos.y) - (BETA * probe.pos.z)};
	// Perturbations and orthogonals
	for (int i = 1; i < NUMPERTURBED + 3; ++i)
	{
		k[i] = {SIGMA * (probe.perturbed[i - 1].pos.y - probe.perturbed[i - 1].pos.x), probe.perturbed[i - 1].pos.x * (ROW - probe.perturbed[i - 1].pos.z) - probe.perturbed[i - 1].pos.y, (probe.perturbed[i - 1].pos.x * probe.perturbed[i - 1].pos.y) - (BETA * probe.perturbed[i - 1].pos.z)};
	}
}

// Initializes the sphere of probes around the probe object, which will then progress into an ellipsoid. Reads from file
void initializePerturbations(Probe *p, int NUMPOINTS)
{
	// Initialize all perturbation probes: Last two are for min lyapunov
	p->perturbed = (Probe *)malloc((NUMPOINTS + 2) * sizeof(Probe));
	for (int i = 0; i < NUMPOINTS + 2; ++i) { copyProbe(&p->perturbed[i], *p); }

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
