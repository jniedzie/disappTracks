//
//  AcrSetProcessor.cpp
//
//  Created by Jeremi Niedziela on 26/03/2019.
//

#include "ArcSetProcessor.hpp"

ArcSetProcessor::ArcSetProcessor() :
pointsProcessor(make_unique<PointsProcessor>()),
circleProcessor(make_unique<CircleProcessor>())
{
  
}

ArcSetProcessor::~ArcSetProcessor()
{
  
}

vector<unique_ptr<ArcSet2D>> ArcSetProcessor::BuildArcSetsFromCircles(const vector<unique_ptr<Circle>> &circles)
{
  vector<unique_ptr<ArcSet2D>> arcs;

  for(auto &circle : circles){

    // re-order circle points if needed:
    double phi0 = circle->GetPointAngle(0);
    double phi1 = circle->GetPointAngle(1);
    double phi2 = circle->GetPointAngle(2);
    
    if((phi1 < phi2 && phi2 < phi0) ||
       (phi1 > phi2 && phi2 > phi0)
       ){
      auto tmp = circle->GetPoint(1);
      circle->points[1] = circle->points[2];
      circle->points[2] = tmp;
    }
    
    if(!IsValidSeed(circle)) continue;
    
    auto arcSet2D = make_unique<ArcSet2D>();
    
    // Add two first points, as they are not added automatically when adding a circle
    arcSet2D->points.push_back(circle->GetPoint(0));
    arcSet2D->points.push_back(circle->GetPoint(1));
    arcSet2D->AddCircle(circle); // this will also add third point
    
    arcs.push_back(move(arcSet2D));
  }
  
  return arcs;
}

bool ArcSetProcessor::IsValidSeed(const unique_ptr<Circle> &circle)
{
  double phiVertex = circle->GetPointAngle(0);
  double phi1      = circle->GetPointAngle(1);
  double phi2      = circle->GetPointAngle(2);
  
  // Reject cases where hits are on the both sides of the vertex point instead of forming a tracklet
  if((phi1 < phiVertex && phi2 > phiVertex) ||
     (phi2 < phiVertex && phi1 > phiVertex)){
    return false; // FILTER
  }
  
  double z0 = circle->GetPoint(0)->GetZ();
  double z1 = circle->GetPoint(1)->GetZ();
  double z2 = circle->GetPoint(2)->GetZ();
  
  // this is the range of where point 0 and point 1 can be located along Z axis
  double z0min = z0 - stripModuleZlength/2.;
  double z0max = z0 + stripModuleZlength/2.;
  double z1min = z1 - stripModuleZlength/2.;
  double z1max = z1 + stripModuleZlength/2.;
  
  // which determines range in the slope:
  double slopeMin = (z1min - z0max)/phi1;
  double slopeMax = (z1max - z0min)/phi1;
  
  // point 2 can be located somewhere between:
  double z2min = z0max + slopeMin * phi2;
  double z2max = z0min + slopeMax * phi2;
  
  if((z2 < z2min) || (z2 > z2max)){
    // this condition kills correct seed... has to be re-thought
    return false; // FILTER
  }
  
  return true;
}

TripletsVector ArcSetProcessor::BuildTripletsCompatibleWithArcSet(const unique_ptr<ArcSet2D> &arcSet,
                                                                  const vector<shared_ptr<Point>> &points)
{
  TripletsVector pointTriplets;
  
  double stripSensorHalfLength = stripModuleZlength/2.;
  
  unique_ptr<Circle> circle            = arcSet->GetLastCircle();
  shared_ptr<Point> point1             = arcSet->GetSecondToLastPoint();
  shared_ptr<Point> point2             = arcSet->GetLastPoint();
  vector<shared_ptr<Point>> pionPoints = arcSet->GetPoints();
  
  for(auto point : points){
    bool isValidPoint = true;
    
    // make sure that it's not the same point as already in the pion track candidate
    if(find(pionPoints.begin(), pionPoints.end(), point) != pionPoints.end()) continue; // FILTER
    
    // new point must be on the correct side of the previous arc in Z direction
    if(fabs(point1->GetZ() - point2->GetZ()) > stripSensorHalfLength){
      // we can check it only if two prevous hits are in different Z locations
      
      if(point1->GetZ() + stripSensorHalfLength < point2->GetZ() - stripSensorHalfLength  &&
         point->GetZ() < (point2->GetZ() - stripSensorHalfLength)){
        isValidPoint = false; // FILTER
      }
      
      if(point1->GetZ() - stripSensorHalfLength > point2->GetZ() + stripSensorHalfLength  &&
         point->GetZ() > (point2->GetZ() + stripSensorHalfLength)){
        isValidPoint = false; // FILTER
      }
    }
    
    // it also has to be within the radius of the helix
    double pointR = pointsProcessor->distanceXY(point, circle->GetCenter());
    if(pointR > 1.1*circle->GetRadius()) isValidPoint = false; // FILTER
    
    if(isValidPoint){
      PointsTriplet points = {point1, point2, point};
      pointTriplets.push_back(points);
    }
  }
  return pointTriplets;
}

vector<shared_ptr<Point>> ArcSetProcessor::FindPossibleNextPoints(const unique_ptr<ArcSet2D> &arcSet,
                                                                  const vector<shared_ptr<Point>> &points)
{
  vector<shared_ptr<Point>> possiblePoints;
  
  double stripSensorHalfLength = stripModuleZlength/2.;
  
  unique_ptr<Circle> circle            = arcSet->GetLastCircle();
  shared_ptr<Point> point1             = arcSet->GetSecondToLastPoint();
  shared_ptr<Point> point2             = arcSet->GetLastPoint();
  vector<shared_ptr<Point>> pionPoints = arcSet->GetPoints();
  
  for(auto point : points){
    bool isValidPoint = true;
    
    // make sure that it's not the same point as already in the pion track candidate
    if(find(pionPoints.begin(), pionPoints.end(), point) != pionPoints.end()) continue; // FILTER
    
    // remove point that are very close to the already existing track (and were not associated with it in prevous iterations)
    bool tooClose = false;
    
    for(auto &c : arcSet->GetCircles()){
      if( c->GetPointAngle(point) > c->GetRange().GetMin() &&
          c->GetPointAngle(point) < c->GetRange().GetMax()){
         
        if(c->GetDistanceToPoint(point) < config->circleThickness){
          tooClose = true;
          break;
        }
      }
    }
    if(tooClose) continue;
    
    // new point must be on the correct side of the previous arc in Z direction
    if(fabs(point1->GetZ() - point2->GetZ()) > stripSensorHalfLength){
      // we can check it only if two prevous hits are in different Z locations
      
      if(point1->GetZ() + stripSensorHalfLength < point2->GetZ() - stripSensorHalfLength  &&
         point->GetZ() < (point2->GetZ() - stripSensorHalfLength)){
        isValidPoint = false; // FILTER
      }
      
      if(point1->GetZ() - stripSensorHalfLength > point2->GetZ() + stripSensorHalfLength  &&
         point->GetZ() > (point2->GetZ() + stripSensorHalfLength)){
        isValidPoint = false; // FILTER
      }
    }
    
    // it also has to be within the radius of the helix
    double pointR = pointsProcessor->distanceXY(point, circle->GetCenter());
    if(pointR > circle->GetRadius() + config->circleThickness) isValidPoint = false; // FILTER
    
    if(isValidPoint) possiblePoints.push_back(point);
  }
  return possiblePoints;
}

unique_ptr<ArcSet2D> ArcSetProcessor::GetBestArcSet(const vector<unique_ptr<ArcSet2D>> &arcSets)
{
  unique_ptr<ArcSet2D> bestArcSet = nullptr;
//  int maxNarcs = 0;
  double bestChi2 = inf;
  
  for(auto &arcSet : arcSets){
    if(arcSet->GetNarcs() < 3) continue;
    
//    if(arcSet->GetNarcs() > maxNarcs){
//      maxNarcs = arcSet->GetNarcs();
//      bestArcSet = make_unique<ArcSet2D>(arcSet);
//    }
    double chi2 = arcSet->GetRadiiSlopeChi2();
    if(chi2 < bestChi2){
      bestChi2 = chi2;
      bestArcSet = make_unique<ArcSet2D>(arcSet);
    }
  }
  return bestArcSet;
}
