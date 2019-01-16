//
//  PointsProcessor.cpp
//
//  Created by Jeremi Niedziela on 16/01/2019.
//

#include "PointsProcessor.hpp"

PointsProcessor::PointsProcessor()
{
  
}

PointsProcessor::~PointsProcessor()
{
  
}

vector<vector<Point>> PointsProcessor::SplitPointsIntoLines(vector<Point> points, double tolerance) const
{
  vector<vector<Point>> pointsByLines;
  bool addedToExisting;
  double toleranceSquared = tolerance*tolerance;
  
  for(Point p : points){
    addedToExisting = false;
    
    // loop over existing lines and check if this point belongs to one of them
    for(vector<Point> &line : pointsByLines){
      // if distance to this line is small enough, just add the point to this line and go to next point
      if(distanceXYsquared(p, Point(line)) < toleranceSquared){
        line.push_back(p);
        addedToExisting = true;
        break;
      }
    }
    if(addedToExisting) continue;
    
    // If the point was not added to any line, create a new line for it
    vector<Point> line;
    line.push_back(p);
    sort(line.begin(), line.end(),[](Point p1, Point p2){return p1.GetZ() < p2.GetZ();});
    pointsByLines.push_back(line);
  }
  
  return pointsByLines;
}

vector<Point> PointsProcessor::GetRandomPoints(int nPoints) const
{
  vector<Point> points;
  double phi, R;
  int layerIndex;
  
  for(int i=0;i<nPoints;i++){
    phi = RandDouble(0, 2*TMath::Pi());
    layerIndex = RandDouble(0, 4);
    R = layerR[layerIndex];
    Point p(R*cos(phi), R*sin(phi), RandDouble(-pixelBarrelZsize, pixelBarrelZsize));
    points.push_back(p);
  }
  return points;
}

double PointsProcessor::distance(Point p1, Point p2) const
{
  return sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));
}

double PointsProcessor::distanceXY(Point p1, Point p2) const
{
  return sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2));
}

double PointsProcessor::distanceXYsquared(Point p1, Point p2) const
{
  return (p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y);
}
