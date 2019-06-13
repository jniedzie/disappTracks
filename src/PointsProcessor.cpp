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
  for(int iLayer=0; iLayer<layerRanges.size(); iLayer++) pointsByLayer.push_back(vector<shared_ptr<Point>>());
  
  for(auto &p : points){
    double pointR = sqrt(pow(p->GetX(), 2) + pow(p->GetY(), 2));
    
    for(int iLayer=0; iLayer<layerRanges.size(); iLayer++){
      if(layerRanges[iLayer].IsInside(pointR)){
        p->SetLayer(iLayer);
        pointsByLayer[iLayer].push_back(p);
        break;
      }
    }
  }
  
  return pointsByLayer;
}

void PointsProcessor::SetPointsLayers(vector<shared_ptr<Point>> &points)
{
  for(auto &p : points){
    double pointR = sqrt(pow(p->GetX(), 2) + pow(p->GetY(), 2));
    
    for(int iLayer=0; iLayer<layerRanges.size(); iLayer++){
      if(layerRanges[iLayer].IsInside(pointR)){
        p->SetLayer(iLayer);
        break;
      }
    }
  }
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

double PointsProcessor::GetPointingAngleXY(const Point &p0, const Point &p1, const Point &p2)
{
  double x_0 = p0.GetX();
  double y_0 = p0.GetY();
  
  double x_1 = p1.GetX();
  double y_1 = p1.GetY();
  
  double x_2 = p2.GetX();
  double y_2 = p2.GetY();
  
  double v1_x = x_1-x_0;
  double v1_y = y_1-y_0;
  
  double v2_x = x_2-x_1;
  double v2_y = y_2-y_1;
  
  // This would give abs(alpha) in fact:
//  double num = v1_x*v2_x + v1_y*v2_y;
//  double den = sqrt(v1_x*v1_x + v1_y*v1_y) * sqrt(v2_x*v2_x + v2_y*v2_y);
//  double alpha = acos(num/den);
  
  // This should give a signed angle
  double alpha2 = atan2(v2_y, v2_x) - atan2(v1_y, v1_x);
  while(alpha2 >=  TMath::Pi()) alpha2 -= 2*TMath::Pi();
  while(alpha2 <= -TMath::Pi()) alpha2 += 2*TMath::Pi();
  
  return alpha2;
}

double PointsProcessor::GetPointingAngleTZ(const Point &p0, const Point &p1, const Point &p2)
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
  
  double v1_x = x_1-x_v;
  double v1_y = y_1-y_v;
  double v1_z = z_1-z_v;
  
  double v2_x = x_2-x_1;
  double v2_y = y_2-y_1;
  double v2_z = z_2-z_1;
  
  // rotate both vectors
  double v1_y_rot = sqrt(v1_x*v1_x + v1_y*v1_y);
  double beta = atan2(v1_x, v1_y);
  double v2_y_rot = sin(beta)*v2_x + cos(beta)*v2_y;
  
  double num = v1_y_rot*v2_y_rot + v1_z*v2_z;
  double den = sqrt(v1_y_rot*v1_y_rot + v1_z*v1_z) * sqrt(v2_y_rot*v2_y_rot + v2_z*v2_z);
  return acos(num/den);
}

Point PointsProcessor::GetPointOnTrack(double L, const Track &track, const Point &eventVertex)
{
  Point p(L * cos(track.GetPhi())    + 10*eventVertex.GetX(),
          L * sin(track.GetPhi())    + 10*eventVertex.GetY(),
          L / tan(track.GetTheta())  + 10*eventVertex.GetZ());
  return p;
}

double PointsProcessor::GetTforPoint(Point &point, const Point &origin, int charge)
{
  double t;
  if(charge < 0) t =  atan2(point.GetY() - origin.GetY(), point.GetX() - origin.GetX());
  else           t = -atan2(point.GetX() - origin.GetX(), point.GetY() - origin.GetY());
  
  t =  atan2(point.GetY() - origin.GetY(), point.GetX() - origin.GetX());
  
  return t;
}


vector<vector<shared_ptr<Point>>> PointsProcessor::RegroupNerbyPoints(const vector<shared_ptr<Point>> &points,
                                                                      double threshold)
{
  vector<vector<shared_ptr<Point>>> regroupedPoints;
  vector<int> alreadyRegroupedIndices;
  
  for(int iPoint1=0; iPoint1<points.size(); iPoint1++){
    if(find(alreadyRegroupedIndices.begin(), alreadyRegroupedIndices.end(), iPoint1) != alreadyRegroupedIndices.end()){
      continue;
    }
    
    vector<shared_ptr<Point>> thisPointNeighbours;
    thisPointNeighbours.push_back(points[iPoint1]);
    
    for(int iPoint2=iPoint1+1; iPoint2<points.size(); iPoint2++){
      if(pointsProcessor.distance(*points[iPoint1], *points[iPoint2]) < threshold){
        thisPointNeighbours.push_back(points[iPoint2]);
        alreadyRegroupedIndices.push_back(iPoint2);
      }
    }
    
    regroupedPoints.push_back(thisPointNeighbours);
  }
  
  return regroupedPoints;
}

void PointsProcessor::SetPointsT(vector<shared_ptr<Point>> &points, const Point &origin, int charge)
{
  int previousPointLayer = -1;
  
  for(int iPoint=0; iPoint<points.size(); iPoint++){
    auto point = points[iPoint];
    double thisPointT = GetTforPoint(*point, origin, charge);
    int thisPointLayer = point->GetLayer();
    
    if(thisPointLayer != previousPointLayer){ // implement corner case!
      auto previousPoint = points[iPoint-1];
      double previousT = previousPoint->GetT();
      if(charge < 0)  while(thisPointT < previousT) thisPointT += 2*TMath::Pi();
      else            while(thisPointT > previousT) thisPointT -= 2*TMath::Pi();
      previousPointLayer = thisPointLayer;
    }
    else if(iPoint!=0){
      auto previousPoint = points[iPoint-1];
      double previousT = previousPoint->GetT();

      if(thisPointT < previousT){
        if(fabs(thisPointT+2*TMath::Pi()-previousT) < fabs(thisPointT-previousT)) thisPointT += 2*TMath::Pi();
      }
      else{
        if(fabs(thisPointT-2*TMath::Pi()-previousT) < fabs(thisPointT-previousT)) thisPointT -= 2*TMath::Pi();
      }
    }
    
    point->SetT(thisPointT);
  }
}

// this is for after turning
void PointsProcessor::SetPointsT(vector<shared_ptr<Point>> &points, const Point &origin, int charge,
                                 shared_ptr<Point> &lastPointBeforeTurning)
{
  int previousPointLayer = lastPointBeforeTurning->GetLayer();
  
  // skip first point as this is the reference
  for(int iPoint=0; iPoint<points.size(); iPoint++){
    auto point = points[iPoint];
    double thisPointT = GetTforPoint(*point, origin, charge);
    int thisPointLayer = point->GetLayer();
    
    if(thisPointLayer != previousPointLayer){ // implement corner case!
      auto previousPoint = points[iPoint-1];
      double previousT = previousPoint->GetT();
      if(charge < 0)  while(thisPointT < previousT) thisPointT += 2*TMath::Pi();
      else            while(thisPointT > previousT) thisPointT -= 2*TMath::Pi();
      previousPointLayer = thisPointLayer;
    }
    else if(iPoint!=0){
      auto previousPoint = points[iPoint-1];
      double previousT = previousPoint->GetT();
      if(thisPointT < previousT)  while(fabs(thisPointT-previousT) > TMath::Pi()/4.) thisPointT += 2*TMath::Pi();
      else                        while(fabs(thisPointT-previousT) > TMath::Pi()/4.) thisPointT -= 2*TMath::Pi();
    }
    else{ // iPoint==0
      auto previousPoint = lastPointBeforeTurning;
      double previousT = previousPoint->GetT();
      if(charge < 0)  while(thisPointT < previousT) thisPointT += 2*TMath::Pi();
      else            while(thisPointT > previousT) thisPointT -= 2*TMath::Pi();
    }
    
    point->SetT(thisPointT);
  }
}

bool PointsProcessor::IsPhiGood(const vector<shared_ptr<Point>> &lastPoints,
                                const vector<shared_ptr<Point>> &secondToLastPoints,
                                const shared_ptr<Point> &point)
{
  bool goodPhi=false;
  
  for(auto &lastPoint : lastPoints){
    for(auto &secondToLastPoint : secondToLastPoints){
      double deltaPhi = pointsProcessor.GetPointingAngleXY(*secondToLastPoint, *lastPoint, *point);
      if(config.nextPointDeltaPhi.IsOutside(fabs(deltaPhi))) continue;
      goodPhi = true;
      break;
    }
    if(goodPhi) break;
  }
  return goodPhi;
}

bool PointsProcessor::IsZgood(const vector<shared_ptr<Point>> &lastPoints,
                              const shared_ptr<Point> &point)
{
  bool goodZ = false;
  
  for(auto &lastPoint : lastPoints){
    double deltaZ = fabs(lastPoint->GetZ() - point->GetZ());
    if(deltaZ > config.nextPointMaxDeltaZ) continue;
    goodZ = true;
    break;
  }
  return goodZ;
}
