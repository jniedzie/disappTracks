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

vector<vector<Point>> PointsProcessor::SplitPointsIntoLines(const vector<shared_ptr<Point>> &points, double tolerance) const
{
  vector<vector<Point>> pointsByLines;
  bool addedToExisting;
  double toleranceSquared = tolerance*tolerance;
  
  for(auto &p : points){
    addedToExisting = false;
    
    // loop over existing lines and check if this point belongs to one of them
    for(vector<Point> &line : pointsByLines){
      // if distance to this line is small enough, just add the point to this line and go to next point
      if(distanceXYsquared(*p, Point(line)) < toleranceSquared){
        line.push_back(*p);
        addedToExisting = true;
        break;
      }
    }
    if(addedToExisting) continue;
    
    // If the point was not added to any line, create a new line for it
    vector<Point> line;
    line.push_back(*p);
    sort(line.begin(), line.end(),[](Point p1, Point p2){return p1.GetZ() < p2.GetZ();});
    pointsByLines.push_back(line);
  }
  
  return pointsByLines;
}

vector<shared_ptr<Point>> PointsProcessor::GetRandomPoints(int nPoints) const
{
  vector<shared_ptr<Point>> points;
  double phi, R;
  
  for(int i=0;i<nPoints;i++){
    phi = RandDouble(0, 2*TMath::Pi());
    R = layerR[RandInt(0, config->nTrackerLayers)];
    auto p = make_shared<Point>(R*cos(phi), R*sin(phi), RandDouble(-pixelBarrelZsize, pixelBarrelZsize));
    points.push_back(p);
  }
  return points;
}

double PointsProcessor::distance(Point p1, Point p2) const
{
  return sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));
}

double PointsProcessor::distance(shared_ptr<Point> p1, shared_ptr<Point> p2) const
{
  return sqrt(pow(p1->x-p2->x,2)+pow(p1->y-p2->y,2)+pow(p1->z-p2->z,2));
}

double PointsProcessor::distanceWithUncertainZ(shared_ptr<Point> p1, shared_ptr<Point> p2,  double zTolerance) const
{
  double zDifference = fabs(p1->z - p2->z);
  
  if(zDifference < zTolerance) zDifference = 0;
  else                         zDifference -= zTolerance;
  
  return sqrt(pow(p1->x-p2->x, 2) + pow(p1->y-p2->y, 2) + pow(zDifference, 2));
}

double PointsProcessor::distanceXY(Point p1, Point p2) const
{
  return sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2));
}

double PointsProcessor::distanceXYsquared(Point p1, Point p2) const
{
  return (p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y);
}

vector<shared_ptr<Point>> PointsProcessor::FilterNearbyPoints(const vector<shared_ptr<Point>> &points, double minPointsSeparation)
{
  vector<shared_ptr<Point>> filteredPoints;
  
  // remove points that are too close to each other
  for(auto &point : points){
    bool tooClose = false;
    
    for(auto otherPoint : filteredPoints){
      if(distance(*point, *otherPoint) < minPointsSeparation) tooClose = true;
    }
    if(tooClose) continue;
    
    filteredPoints.push_back(point);
  }
  return filteredPoints;
}

TripletsVector PointsProcessor::BuildPointTriplets(const vector<shared_ptr<Point>> &points)
{
  int nPoints = (int)points.size();
  TripletsVector pointTriplets;
  
  for(int i=0;i<nPoints;i++){
    for(int j=i+1;j<nPoints;j++){
      PointsTriplet triplet = {make_shared<Point>(inf,inf,inf), points[i], points[j] };
      pointTriplets.push_back(triplet);
    }
  }
  return pointTriplets;
}

TripletPairsVector PointsProcessor::BuildPointTripletPairs(const vector<shared_ptr<Point>> &points,
                                                       const shared_ptr<Point> &originMin,
                                                       const shared_ptr<Point> &originMax)
{
  int nPoints = (int)points.size();
  TripletPairsVector pointTripletPairs;
  
  for(int i=0;i<nPoints;i++){
    for(int j=i+1;j<nPoints;j++){
      PointsTriplet tripletMin = {originMin, points[i], points[j] };
      PointsTriplet tripletMax = {originMax, points[i], points[j] };
      pointTripletPairs.push_back(make_pair(tripletMin, tripletMax));
    }
  }
  return pointTripletPairs;
}
