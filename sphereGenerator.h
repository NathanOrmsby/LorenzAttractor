/*
 * sphereGenerator.h
 *
 *  Created on: Feb 8, 2023
 *      Author: norms
 */

#ifndef SPHEREGENERATOR_H_
#define SPHEREGENERATOR_H_

#include <string>

void fibonacciSphere(int numPoints, double arr[][3], double r);
void toFile(std::string fileName, double data[][3], int dataLen);
void writeSphereToFile(int numPoints, double radius);
void writeCircleToFile(int numPoints, double radius);
void uniformCircle(int numPoints, double arr[][3], double r);


#endif /* SPHEREGENERATOR_H_ */
