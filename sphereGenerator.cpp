/*
 * sphereGenerator.cpp
 *
 *  Created on: Feb 8, 2023
 *      Author: norms
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>

#include "..\headers\sphereGenerator.h"


void writeSphereToFile(int numPoints, double radius)
{
	double data[numPoints][3];
	fibonacciSphere(numPoints, data, radius);

	std::string fname = "unitSphere.csv";
	toFile(fname, data, numPoints);
}

void writeCircleToFile(int numPoints, double radius)
{
	double data[numPoints][3];
	uniformCircle(numPoints, data, radius);

	std::string fname = "unitCircle.csv";
	toFile(fname, data, numPoints);
}

void uniformCircle(int numPoints, double arr[][3], double r)
{
	double d = (2.0 * M_PI) / (double)numPoints;
	for (int i = 0; i < numPoints; ++i)
	{
		double theta = i * d;
		arr[i][0] = r * cos(theta); arr[i][1] = 0; arr[i][2] = r * sin(theta);
	}
}

void fibonacciSphere(int numPoints, double arr[][3], double r)
{
	double phi = M_PI * (3.0 - sqrt(5.0));

	for (int i = 0; i < numPoints; i++)
	{
		double y = 1 - (i / (double)(numPoints - 1)) * 2;
		double radius = sqrt(1 - y * y);

		double theta = phi * i;
		double x = cos(theta) * radius;
		double z = sin(theta) * radius;

		arr[i][0] = r * x; arr[i][1] = r * y; arr[i][2] = r * z;
	}
}

void toFile(std::string fileName, double data[][3], int dataLen)
{
	std::ofstream file;
	file.open(fileName);

//	file << "x" << "," << "y" << "," << "z" << std::endl;
	for (int i = 0; i < dataLen; ++i)
	{
		file << std::to_string(data[i][0]) << "," << std::to_string(data[i][1]) << "," << std::to_string(data[i][2]) << std::endl;
	}
	file.close();
}


