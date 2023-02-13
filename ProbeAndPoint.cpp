/*
 * ProbeAndPoint.cpp
 *
 *  Created on: Feb 12, 2023
 *      Author: norms
 */


#include "..\headers\ProbeAndPoint.h"
#include "..\headers\vector3.h"

// Utility functions for probes and points

// Points
// Copies data from i to o
void copyPoint(Point *o, Point i)
{
	o->pos.x = i.pos.x; o->pos.y = i.pos.y; o->pos.z = i.pos.z;
}

// Probes
// Copies all data from one probe to another, doesn't worry about perturbed pointer
void copyProbe(Probe *o, Probe i)
{
	o->pos.x = i.pos.x; o->pos.y = i.pos.y; o->pos.z = i.pos.z;
}

