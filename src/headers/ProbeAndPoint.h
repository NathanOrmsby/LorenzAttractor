/*
 * ProbeAndPoint.h
 *
 *  Created on: Feb 12, 2023
 *      Author: norms
 */

#ifndef HEADERS_PROBEANDPOINT_H_
#define HEADERS_PROBEANDPOINT_H_

#include "vector3.h"

// Contains perturbations
typedef struct Probe
{
	vector3 pos;
	Probe *perturbed;
	Probe *h2;
	Probe *h3;
} Probe;

// No perturbations
typedef struct
{
	vector3 pos;
} Point;


// Utility functions for probes and points

// Points
// Copies data from i to o
void copyPoint(Point *o, Point i);
// Probes
// Copies velocity data from a probe to its perturbed probe. Perturbed probe does not have further perturbations
void copyProbeToPerturbed(Probe *o, Probe i);
// Copies all data from one probe to another, doesn't worry about perturbed pointer
void copyProbe(Probe *o, Probe i);



#endif /* HEADERS_PROBEANDPOINT_H_ */
