//  PointsProcessor.cpp
//
//  Created by Jeremi Niedziela on 16/01/2019.

#include "PointsProcessor.hpp"

PointsProcessor pointsProcessor = PointsProcessor();

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
    R = layerR[RandInt(0, config.nTrackerLayers)];
    auto p = make_shared<Point>(R*cos(phi), R*sin(phi), RandDouble(-pixelBarrelZsize, pixelBarrelZsize));
    points.push_back(p);
  }
  return points;
}

double PointsProcessor::distance(const Point &p1,const Point &p2) const
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

vector<PointsPair> PointsProcessor::BuildPointPairs(const vector<shared_ptr<Point>> &points)
{
  int nPoints = (int)points.size();
  vector<PointsPair> pointPairs;
  
  for(int i=0;i<nPoints;i++){
    for(int j=i+1;j<nPoints;j++){
      auto pointPair =  make_pair( points[i], points[j] );
      pointPairs.push_back(pointPair);
    }
  }
  return pointPairs;
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


vector<vector<shared_ptr<Point>>> PointsProcessor::SortByLayer(const vector<shared_ptr<Point>> &points)
{
  vector<vector<shared_ptr<Point>>> pointsByLayer;
  for(int iLayer=0; iLayer<nLayers; iLayer++) pointsByLayer.push_back(vector<shared_ptr<Point>>());
  
  for(auto &p : points){
    // find in which layer is the hit located
    int layer = -1;
    double minDist = inf;
    double pointR = sqrt(pow(p->GetX(), 2) + pow(p->GetY(), 2));
    
    for(int iLayer=0; iLayer<nLayers; iLayer++){
      double pointLayerDist = fabs(layerR[iLayer] - pointR);
      
      if(pointLayerDist < minDist){
        minDist = pointLayerDist;
        layer = iLayer;
      }
    }
    p->SetLayer(layer);
    pointsByLayer[layer].push_back(p);
  }
  
  return pointsByLayer;
}

double PointsProcessor::GetPointingAngle(const Point &p0, const Point &p1, const Point &p2)
{
  double x_v = p0.GetX();
  double y_v = p0.GetY();
  double z_v = p0.GetZ();
  
  double x_1 = p1.GetX();
  double y_1 = p1.GetY();
  double z_1 = p1.GetZ();
  
  double x_2 = p2.GetX();
  double y_2 = p2.GetY();
  double z_2 = p2.GetZ();
  
  double v1_x = 2*x_1-x_v;
  double v1_y = 2*y_1-y_v;
  double v1_z = 2*z_1-z_v;
  
  double v2_x = x_2-x_1;
  double v2_y = y_2-y_1;
  double v2_z = z_2-z_1;
  
  double num = v1_x*v2_x + v1_y*v2_y + v1_z*v2_z;
  double den = sqrt(v1_x*v1_x + v1_y*v1_y + v1_z*v1_z) * sqrt(v2_x*v2_x + v2_y*v2_y + v2_z*v2_z);
  return acos(num/den);
}

Point PointsProcessor::GetPointOnTrack(double L, const Track &track, const Point &eventVertex)
{
  Point p(L * cos(track.GetPhi())    + 10*eventVertex.GetX(),
          L * sin(track.GetPhi())    + 10*eventVertex.GetY(),
          L / tan(track.GetTheta())  + 10*eventVertex.GetZ());
  return p;
}

