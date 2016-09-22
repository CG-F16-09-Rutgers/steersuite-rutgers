//
// Copyright (c) 2009-2015 Glen Berseth, Mubbasir Kapadia, Shawn Singh, Petros Faloutsos, Glenn Reinman
// See license.txt for complete license.
// Copyright (c) 2015 Mahyar Khayatkhoei
//

#include <algorithm>
#include <vector>
#include <util/Geometry.h>
#include <util/Curve.h>
#include <util/Color.h>
#include <util/DrawLib.h>
#include "Globals.h"

using namespace Util;

Curve::Curve(const CurvePoint& startPoint, int curveType) : type(curveType)
{
	controlPoints.push_back(startPoint);
}

Curve::Curve(const std::vector<CurvePoint>& inputPoints, int curveType) : type(curveType)
{
	controlPoints = inputPoints;
	sortControlPoints();
}

// Add one control point to the vector controlPoints
void Curve::addControlPoint(const CurvePoint& inputPoint)
{
	controlPoints.push_back(inputPoint);
	sortControlPoints();
}

// Add a vector of control points to the vector controlPoints
void Curve::addControlPoints(const std::vector<CurvePoint>& inputPoints)
{
	for (int i = 0; i < inputPoints.size(); i++)
		controlPoints.push_back(inputPoints[i]);
	sortControlPoints();
}

// Draw the curve shape on screen, usign window as step size (bigger window: less accurate shape)
void Curve::drawCurve(Color curveColor, float curveThickness, int window)
{
#ifdef ENABLE_GUI
	if (!checkRobust()) {
		return;
	}

	Point start, end;
	float endTime = controlPoints[controlPoints.size() - 1].time;
	for (int t = 0; t < endTime; t += window) {
		calculatePoint(start, t);
		if (t + window >= endTime) {
			calculatePoint(end, endTime - 0.25);
		} else {
			calculatePoint(end, t + window);
		}
		DrawLib::drawLine( start, end, curveColor, curveThickness);
	}

	// Robustness: make sure there is at least two control point: start and end points
	// Move on the curve from t=0 to t=finalPoint, using window as step size, and linearly interpolate the curve points
	// Note that you must draw the whole curve at each frame, that means connecting line segments between each two points on the curve
	
	return;
#endif
}

// Sort controlPoints vector in ascending order: min-first
void Curve::sortControlPoints()
{
	for (int i = 1; i < controlPoints.size(); i++) {
		int j = i;
		while (j > 0 && controlPoints[j - 1].time >= controlPoints[j].time) {
			if (controlPoints[j - 1].time == controlPoints[j].time) {
				controlPoints.erase(controlPoints.begin() + j);
				i--;
				break;
			}
			Util::CurvePoint temp = controlPoints[j - 1];
			controlPoints[j - 1] = controlPoints[j];
			controlPoints[j] = temp;
			j--;
		}
	}

	return;
}

// Calculate the position on curve corresponding to the given time, outputPoint is the resulting position
// Note that this function should return false if the end of the curve is reached, or no next point can be found
bool Curve::calculatePoint(Point& outputPoint, float time)
{
	// Robustness: make sure there is at least two control point: start and end points
	if (!checkRobust())
		return false;

	// Define temporary parameters for calculation
	unsigned int nextPoint;

	// Find the current interval in time, supposing that controlPoints is sorted (sorting is done whenever control points are added)
	// Note that nextPoint is an integer containing the index of the next control point
	if (!findTimeInterval(nextPoint, time))
		return false;

	// Calculate position at t = time on curve given the next control point (nextPoint)
	if (type == hermiteCurve)
	{
		outputPoint = useHermiteCurve(nextPoint, time);
	}
	else if (type == catmullCurve)
	{
		outputPoint = useCatmullCurve(nextPoint, time);
	}

	// Return
	return true;
}

// Check Roboustness
bool Curve::checkRobust()
{
	if (controlPoints.size() < 2) {
		return false;
	}

	return true;
}

// Find the current time interval (i.e. index of the next control point to follow according to current time)
bool Curve::findTimeInterval(unsigned int& nextPoint, float time)
{
	for (int i = 0; i < controlPoints.size(); i++) {
		if (controlPoints[ i].time > time) {
			nextPoint = i;
			return true;
		}
	}

	return false;
}

// Implement Hermite curve
Point Curve::useHermiteCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;
	float normalTime, intervalTime;

	float blending[4];
	float G[4][3];
	float point[3];
	
	intervalTime = controlPoints[nextPoint].time - controlPoints[nextPoint - 1].time;
	normalTime = (time - controlPoints[nextPoint - 1].time) / intervalTime;
	
	float t2 = normalTime * normalTime;
	float t3 = t2 * normalTime;
	blending[0] = 2*t3 - 3*t2 + 1;
	blending[1] = -2*t3 + 3*t2;
	blending[2] = t3 - 2*t2 + normalTime;
	blending[3] = t3 - t2;

	Point tangent1 = Point(controlPoints[nextPoint - 1].tangent.x, controlPoints[nextPoint - 1].tangent.y, controlPoints[nextPoint - 1].tangent.z);
	Point tangent2 = Point(controlPoints[nextPoint].tangent.x, controlPoints[nextPoint].tangent.y, controlPoints[nextPoint].tangent.z);
	
	newPosition = blending[0] * controlPoints[nextPoint - 1].position;
	newPosition = newPosition + blending[1] * controlPoints[nextPoint].position;
	newPosition = newPosition + blending[2] * tangent1 * intervalTime;
	newPosition = newPosition + blending[3] * tangent2 * intervalTime;
	
	// Return result
	return newPosition;
}

// Implement Catmull-Rom curve
Point Curve::useCatmullCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;

	unsigned int nextPointq = nextPoint - 2;  //I made different assumptions about nextPoint than my partner, so I gotta use this

	float t[4];
	Point cp[4];
	int i = 0; //edge cases
	int max = 4;//used to handle edge cases of Catmull

	if (nextPointq == controlPoints.size() - 1) //Shouldn't ever run, I think, since it'll be the end of drawing, but just in case
		return controlPoints[nextPointq].position;

	if (nextPointq == 0) // edge case; first control point. Approach: add extra control point
	{
		cp[0] = -1 * (controlPoints[0].position + controlPoints[1].position) / 2;
		t[0] = -1 * (controlPoints[0].time + controlPoints[1].time) / 2;
		i++;
	}
	if (nextPointq == controlPoints.size() - 2)
	{
		cp[3] = -1 * (controlPoints[nextPointq].position + controlPoints[nextPointq - 1].position) / 2;
		t[3] = -1 * (controlPoints[nextPointq].time + controlPoints[nextPointq - 1].time) / 2;
		max--;
	}
	while (i < max) {
		t[i] = controlPoints[nextPointq - 1 + i].time;
		cp[i] = controlPoints[nextPointq - 1 + i].position;
		i++;
	}
	if ((int)t[0] + (int)t[1] == 0)
		t[0] = t[0] + 20;
	//std::cout << t[0] << " " << t[1] << " " << t[2] << " " << t[3] << " " << cp[0] << " " << cp[1] << " " << cp[2] << " " << cp[3] << " " << std::endl;
	Point A1 = (t[1] - time) / (t[1] - t[0])*cp[0] + (time - t[0]) / (t[1] - t[0])*cp[1];
	Point A2 = (t[2] - time) / (t[2] - t[1])*cp[1] + (time - t[1]) / (t[2] - t[1])*cp[2];
	Point A3 = (t[3] - time) / (t[3] - t[2])*cp[2] + (time - t[2]) / (t[3] - t[2])*cp[3];

	Point B1 = (t[2] - time) / (t[2] - t[0])*A1 + (time - t[0]) / (t[2] - t[0])*A2;
	Point B2 = (t[3] - time) / (t[3] - t[1])*A2 + (time - t[1]) / (t[3] - t[1])*A3;

	newPosition = (t[2] - time) / (t[2] - t[1])*B1 + (time - t[1]) / (t[2] - t[1])*B2;
	//std::cout << "Current time: " << time << "Current x pos: " << newPosition.x << std::endl;
	//std::cout << nextPointq << std::endl;

	// Calculate position at t = time on Catmull-Rom curve

	// Return result
	/*
	std::cout << newPosition.x << " " << newPosition.y << " " << newPosition.z << " " << std::endl;
	std::cout << cp[0].x << " " << cp[1].x << " " << cp[2].x << " " << cp[3].x << " " << std::endl;
	std::cout << t[0] << " " << t[1] << " " << t[2] << " " << t[3] << " " << std::endl;

	for (int i = 0; i < controlPoints.size(); i++) {
	std::cout << controlPoints.size() << " Info controlpoints: x pos: " << controlPoints.at(i).position.x << std::endl;
	std::cout << controlPoints.size() << " Info controlpoints: time: " << controlPoints.at(i).time << std::endl;
	//std::cout << controlPoints[nextPointq].time << std::endl;

	}
	*/
	return newPosition;
}