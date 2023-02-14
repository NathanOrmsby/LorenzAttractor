/*
 * singleProbeChaosLorenz.h
 *
 *  Created on: Feb 12, 2023
 *      Author: norms
 */

#ifndef HEADERS_SINGLEPROBECHAOSLORENZ_H_
#define HEADERS_SINGLEPROBECHAOSLORENZ_H_

#include "ProbeAndPoint.h"

// Main function, probe with sphere of perturbations
double *singleProbeChaosLorenz(Probe probe, double dt, int totalSteps, int ITERATIONS, int NUMPERTURBED, int RADIUS);
// Initializes the sphere of probes around the probe object, which will then progress into an ellipsoid. Reads from file
void initializePerturbations(Probe *p, int NUMPOINTS);
// RK4 for probe and its perturbations in a lorenz system
void probeRK4(Probe *probe, int NUMPERTURBED, double dt);
// Stores lorenz change on x, y, z using lorenz equations: Uses Lorenz attract initial conditions for single probe and its perturbations
void extractLorenzKProbe(Probe probe, int NUMPERTURBED, vector3 *k);
// Single step for RK4: single probe with perturbations and two orthogonals
void probeRK4Step(Probe *probe, int NUMPERTURBED, vector3 *k, double dt);
// Handles the timestep processes for calculating max and min lyapunov exponents
void lyapunovChaosStuffSingle(Probe probe, int NUMPERTURBED, double RADIUS, double *LfSums, double *areaSum, double initialArea, double *volumeSum, double initialVolume);
// Gram schmidt Renormalization for single probe: Renormalizes h2 vectors into two orthogonal unit vectors that preserve direction
void gramSchmidtRenormalizationh2Single(Probe *probe, int NUMPERTURBED, vector3 v1, vector3 v2);
// Gram schmidt Renormalization for single probe: Adds to the volume sum and renormalizes h3 vectors into three orthogonal unit vectors that preserve direction
void gramSchmidtRenormalizationh3Single(Probe *probe, int NUMPERTURBED, double *volumeSum, double logArea);
// Extracts minimum and maximum lyapunov for single probe with perturbations
void extractLyapunovSingle(double *LfSums, double areaSum, double volumeSum, int NUMPERTURBED, double *lyapunovs, double t);
// Write probe and perturbation x,y,z values to file
void perturbationsToFile(vector3 **data, int dataLen, int NUMPERTURBED);


#endif /* HEADERS_SINGLEPROBECHAOSLORENZ_H_ */
