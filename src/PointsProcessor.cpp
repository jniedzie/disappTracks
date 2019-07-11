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

double PointsProcessor::distance(const Point &p1,const Point &p2) const
{
  return sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2)+pow(p1.z-p2.z,2));
}

double PointsProcessor::distance(shared_ptr<Point> p1, shared_ptr<Point> p2) const
{
  return sqrt(pow(p1->x-p2->x,2)+pow(p1->y-p2->y,2)+pow(p1->z-p2->z,2));
}

double PointsProcessor::distanceXY(Point p1, Point p2) const
{
  return sqrt(pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2));
}

vector<vector<shared_ptr<Point>>> PointsProcessor::SortByLayer(const vector<shared_ptr<Point>> &points)
{
  vector<vector<shared_ptr<Point>>> pointsByLayer;
  for(int iLayer=0; iLayer<layerRanges.size(); iLayer++) pointsByLayer.push_back(vector<shared_ptr<Point>>());
  
  for(auto &p : points){
    if(p->GetSubDetName() == "P1PXEC" ||
       p->GetSubDetName() == "TID" ||
       p->GetSubDetName() == "TEC"){
      continue;
    }
    
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

vector<vector<shared_ptr<Point>>> PointsProcessor::SortByDisk(const vector<shared_ptr<Point>> &points)
{
  vector<vector<shared_ptr<Point>>> pointsByDisk;
  for(size_t iDisk= 0; iDisk<2*diskRanges.size()+1; iDisk++){
    pointsByDisk.push_back(vector<shared_ptr<Point>>());
  }
  
  for(auto &p : points){
    if(p->GetSubDetName() == "P1PXB" ||
       p->GetSubDetName() == "TIB" ||
       p->GetSubDetName() == "TOB"){
      continue;
    }
    
    double pointZ = p->GetZ();
    int signZ = sgn(p->GetZ());
    
    for(size_t iDisk=0; iDisk<diskRanges.size(); iDisk++){
      range<double> diskRange(diskRanges[iDisk].front().GetMin(),
                              diskRanges[iDisk].back().GetMax());
      
      if(diskRange.IsInside(fabs(pointZ))){
        p->SetDisk((int)(signZ*(iDisk+1)));
        pointsByDisk[GetDisksArrayIndex((int)iDisk, signZ)].push_back(p);
        break;
      }
    }
  }
  
  return pointsByDisk;
}

void PointsProcessor::SetPointsLayers(vector<shared_ptr<Point>> &points)
{
  for(auto &p : points){
    if(p->GetSubDetName() == "P1PXEC" ||
       p->GetSubDetName() == "TID" ||
       p->GetSubDetName() == "TEC"){
      continue;
    }
    
    double pointR = sqrt(pow(p->GetX(), 2) + pow(p->GetY(), 2));
    
    for(int iLayer=0; iLayer<layerRanges.size(); iLayer++){
      if(layerRanges[iLayer].IsInside(pointR)){
        p->SetLayer(iLayer);
        break;
      }
    }
  }
}

void PointsProcessor::SetPointsDisks(vector<shared_ptr<Point>> &points)
{
  for(auto &p : points){
    if(p->GetSubDetName() == "P1PXB" ||
       p->GetSubDetName() == "TIB" ||
       p->GetSubDetName() == "TOB"){
      continue;
    }
    
    double pointZ = p->GetZ();
    int signZ = sgn(p->GetZ());
    
    for(int iDisk=0; iDisk<diskRanges.size(); iDisk++){
      range<double> diskRange(diskRanges[iDisk].front().GetMin(),
                              diskRanges[iDisk].back().GetMax());
      
      if(diskRange.IsInside(fabs(pointZ))){
        p->SetDisk((int)(signZ*(iDisk+1)));
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
  double v1_x = p1.GetX()-p0.GetX();
  double v1_y = p1.GetY()-p0.GetY();
  
  double v2_x = p2.GetX()-p1.GetX();
  double v2_y = p2.GetY()-p1.GetY();
  
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

double PointsProcessor::GetTforPoint(const Point &point, const Point &origin, int charge)
{
  double t = atan2(point.GetY() - origin.GetY(), point.GetX() - origin.GetX());

  // not sure if this should depend on charge or not...
  //  if(charge < 0) t =  atan2(point.GetY() - origin.GetY(), point.GetX() - origin.GetX());
  //  else           t = -atan2(point.GetX() - origin.GetX(), point.GetY() - origin.GetY());
  
  return t;
}


vector<vector<shared_ptr<Point>>> PointsProcessor::RegroupNerbyPoints(const vector<shared_ptr<Point>> &points)
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
      if(pointsProcessor.distance(*points[iPoint1], *points[iPoint2]) < config.params["double_hit_max_distance"]){
        thisPointNeighbours.push_back(points[iPoint2]);
        alreadyRegroupedIndices.push_back(iPoint2);
      }
    }
    
    regroupedPoints.push_back(thisPointNeighbours);
  }
  
  return regroupedPoints;
}

bool PointsProcessor::IsPhiGood(const vector<shared_ptr<Point>> &lastPoints,
                                const vector<shared_ptr<Point>> &secondToLastPoints,
                                const shared_ptr<Point> &point,
                                int charge)
{
  bool goodPhi=false;
  
  for(auto &lastPoint : lastPoints){
    for(auto &secondToLastPoint : secondToLastPoints){
      double deltaPhi = pointsProcessor.GetPointingAngleXY(*secondToLastPoint, *lastPoint, *point);
      
      range<double> nextPointDeltaPhi(-inf, inf);
      
      if(config.params["do_asymmetric_constraints"]){
        deltaPhi = charge*deltaPhi;
        nextPointDeltaPhi = range<double>(config.params["next_point_min_delta_phi"],
                                          config.params["next_point_max_delta_phi"]);
      }
      else{
        deltaPhi = fabs(deltaPhi);
        nextPointDeltaPhi = range<double>(0, config.params["next_point_max_delta_phi"]);
      }
      
      if(nextPointDeltaPhi.IsOutside(deltaPhi)) continue;
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
    if(deltaZ > config.params["next_point_max_delta_z"]) continue;
    goodZ = true;
    break;
  }
  return goodZ;
}

bool PointsProcessor::IsTgood(const vector<double> &lastPointsT, double pointT)
{
  bool goodT = false;
  
  for(auto &lastPointT : lastPointsT){
    double deltaT = fabs(lastPointT - pointT);
    if(deltaT > config.params["next_point_max_delta_t"]) continue;
    goodT = true;
    break;
  }
  return goodT;
}

bool PointsProcessor::IsGoodMiddleSeedHit(const Point &point,
                                          const Point &trackMidPoint,
                                          const Point &eventVertex,
                                          int charge)
{
  double middleHitDeltaPhi = GetPointingAngleXY(eventVertex, trackMidPoint, point);
  
  range<double> seedMiddleHitDeltaPhi(-inf, inf);
  
  if(config.params["do_asymmetric_constraints"]){
    middleHitDeltaPhi = charge*middleHitDeltaPhi;
    seedMiddleHitDeltaPhi = range<double>(config.params["seed_middle_hit_min_delta_phi"],
                                          config.params["seed_middle_hit_max_delta_phi"]);
  }
  else{
    middleHitDeltaPhi = fabs(middleHitDeltaPhi);
    seedMiddleHitDeltaPhi = range<double>(0, config.params["seed_middle_hit_max_delta_phi"]);
  }
  
  if(seedMiddleHitDeltaPhi.IsOutside(middleHitDeltaPhi)){
    if(config.params["verbosity_level"]>1) cout<<"Seed middle hit Δφ too large"<<endl;
    return false;
  }
  
  double middleHitDeltaZ = fabs(point.GetZ() - trackMidPoint.GetZ());
  if(middleHitDeltaZ > config.params["seed_middle_hit_max_delta_z"]){
    if(config.params["verbosity_level"]>1) cout<<"Seed middle hit Δz too large"<<endl;
    return false;
  }
  return true;
}


bool PointsProcessor::IsGoodLastSeedHit(const Point &point,
                                        const Point &middlePoint,
                                        const Point &trackMidPoint,
                                        int charge)
{
  double lastHitDeltaPhi = pointsProcessor.GetPointingAngleXY(trackMidPoint, middlePoint, point);
  
  range<double> seedLastHitDeltaPhi(-inf, inf);
  
  if(config.params["do_asymmetric_constraints"]){
    lastHitDeltaPhi = charge*lastHitDeltaPhi;
    seedLastHitDeltaPhi = range<double>(config.params["seed_last_hit_min_delta_phi"],
                                        config.params["seed_last_hit_max_delta_phi"]);
  }
  else{
    lastHitDeltaPhi = fabs(lastHitDeltaPhi);
    seedLastHitDeltaPhi = range<double>(0, config.params["seed_last_hit_max_delta_phi"]);
  }
  
  if(seedLastHitDeltaPhi.IsOutside(lastHitDeltaPhi) && !point.IsEndcapHit()){
    if(config.params["verbosity_level"]>1) cout<<"Seed last hit Δφ too large"<<endl;
    return false;
  }
  
  double lastPointDeltaZ = fabs(middlePoint.GetZ() - point.GetZ());
  if(lastPointDeltaZ > config.params["seed_last_hit_max_delta_z"]){
    if(config.params["verbosity_level"]>1) cout<<"Seed last hit Δz too large"<<endl;
    return false;
  }
  return true;
}

vector<shared_ptr<Point>> PointsProcessor::GetGoodMiddleSeedHits(const vector<shared_ptr<Point>> &middlePoints,
                                                                 const Point &trackMidPoint,
                                                                 const Point &eventVertex,
                                                                 int charge)
{
  vector<shared_ptr<Point>> goodMiddlePoints;
  for(auto &point : middlePoints){
    if(!pointsProcessor.IsGoodMiddleSeedHit(*point, trackMidPoint, eventVertex, charge)) continue;
    goodMiddlePoints.push_back(point);
  }
  return goodMiddlePoints;
}

vector<shared_ptr<Point>> PointsProcessor::GetGoodLastSeedHits(const vector<shared_ptr<Point>> &lastPoints,
                                                               const vector<shared_ptr<Point>> &goodMiddlePoints,
                                                               const Point &trackMidPoint,
                                                               int charge)
{
  vector<shared_ptr<Point>> goodLastPoints;
  for(auto &point : lastPoints){
    for(auto &middlePoint : goodMiddlePoints){
      if(!pointsProcessor.IsGoodLastSeedHit(*point, *middlePoint, trackMidPoint, charge)) continue;
      goodLastPoints.push_back(point);
      break;
    }
  }
  return goodLastPoints;
}
