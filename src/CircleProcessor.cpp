//
//  CircleProcessor.cpp
//
//  Created by Jeremi Niedziela on 26/03/2019.
//

#include "CircleProcessor.hpp"

CircleProcessor::CircleProcessor() :
pointsProcessor(make_unique<PointsProcessor>())
{
  
}

CircleProcessor::~CircleProcessor()
{
  
}

void CircleProcessor::RemoveSimilarCircles(vector<unique_ptr<Circle>> &circles)
{
  if(circles.size() == 0) return;
  
  sort(circles.begin(), circles.end(),
       [](const auto &c1, const auto &c2) -> bool {return c1->GetRadius() < c2->GetRadius();});
  
  for(uint i=0; i<circles.size()-1; i++){
    if(fabs(circles[i]->GetRadius() - circles[i+1]->GetRadius()) < config->circleThickness){
      if((circles[i]->GetRange().GetMax()   - circles[i]->GetRange().GetMin()) >
         (circles[i+1]->GetRange().GetMax() - circles[i+1]->GetRange().GetMin())){
        circles.erase(circles.begin()+i);
      }
      else{
        circles.erase(circles.begin()+i+1);
      }
      i--;
    }
  }
}

vector<unique_ptr<Circle>> CircleProcessor::BuildCirclesFromPoints(const TripletsVector &points)
{
  vector<unique_ptr<Circle>> circles;
  
  for(auto &p : points){
    auto circle = GetCircleFromTriplet(p);
    circles.push_back(move(circle));
  }
  
  return circles;
}

vector<unique_ptr<Circle>> CircleProcessor::BuildCirclesFromTripletPairs(const TripletPairsVector &triplets,
                                                                         range<double> radiusRange)
{
  vector<unique_ptr<Circle>> circles;
  
  for(auto &p : triplets){
    auto circleMin = GetCircleFromTriplet(p.first);
    auto circleMax = GetCircleFromTriplet(p.second);
    
    // remove circles that are outside of the allowed range
    if((circleMin->GetRadius() < radiusRange.GetMin()) ||
       (circleMin->GetRadius() > radiusRange.GetMax())){
      circleMin = nullptr;
    }
    if((circleMax->GetRadius() < radiusRange.GetMin()) ||
       (circleMax->GetRadius() > radiusRange.GetMax())){
      circleMax = nullptr;
    }
    
    // if no circles survived, just skip this triplets pair
    if(!circleMin && !circleMax) continue;
    
    // if only one survived, it's already the better one
    if(!circleMin &&  circleMax){
      circles.push_back(move(circleMax));
      continue;
    }
    if( circleMin && !circleMax){
      circles.push_back(move(circleMin));
      continue;
    }
  
    // then, if both survived, keep the one with bigger radius
    if( circleMin->GetRadius() > circleMax->GetRadius() ) circles.push_back(move(circleMin));
    else                                                  circles.push_back(move(circleMax));
  }
  
  return circles;
}

unique_ptr<Circle> CircleProcessor::BuildCircleFromParams(const double *par,
                                                          const Point &vertex,
                                                          const Track &track)
{
  double R  = par[0];
  double px = par[1];
  double py = par[2];
  
  double phi   = track.GetPhi();
  double theta = track.GetTheta();
  
  double x0 = R*cos(phi)              + 10*vertex.GetX();
  double y0 = R*sin(phi)              + 10*vertex.GetY();
  double z0 = R/sin(theta)*cos(theta) + 10*vertex.GetZ();
  
  auto decayPoint  = Point(x0,y0,z0);
  auto momentum    = make_unique<Point>(px,py,0);
  
  return make_unique<Circle>(decayPoint, momentum);
}


unique_ptr<Circle> CircleProcessor::GetMostCompatibleCircle(const vector<unique_ptr<Circle>> &circles,
                                                            const unique_ptr<Circle> &theCircle,
                                                            vector<double> &alphaVector)
{
  double bestRadiiDifference = inf;
  unique_ptr<Circle> bestCircle = nullptr;
  
  // Calculate parameters of the circle
  double p_x  = theCircle->GetPoints()[2]->GetX();
  double p_y  = theCircle->GetPoints()[2]->GetY();
  double c1_x = theCircle->GetCenter().GetX();
  double c1_y = theCircle->GetCenter().GetY();
  double r1_x = c1_x - p_x;
  double r1_y = c1_y - p_y;
  double r1_mod = sqrt(r1_x*r1_x + r1_y*r1_y);
  
  for(auto &testingCircle : circles){
    // New track segment cannot have greater radius (within some tolerance)
    if(testingCircle->GetRadius() > 1.01*theCircle->GetRadius()) continue; // FILTER
  
    // make sure that the new circle doesn't go in an opposite direction
    if(testingCircle->GetRange().GetMin() > theCircle->GetRange().GetMin() ||
       testingCircle->GetRange().GetMax() > theCircle->GetRange().GetMin()){
//      continue;
    }
    
    
    // The center of the new circle should be withing the previous circle
    double centerDifference = pointsProcessor->distanceXY(theCircle->GetCenter(), testingCircle->GetCenter());
    if(centerDifference > theCircle->GetRadius()) continue; // FILTER
    
    // Here we calculate an angle between radius of the previous circle and radius of the testing circle
    // looking from the last point of previous circle. For perfectly matching circles that would be zero
    double r2_x   = testingCircle->GetCenter().GetX() - p_x;
    double r2_y   = testingCircle->GetCenter().GetY() - p_y;
    double r2_mod = sqrt(r2_x*r2_x + r2_y*r2_y);
    double alpha  = acos((r1_x*r2_x + r1_y*r2_y) / (r1_mod*r2_mod));
    alphaVector.push_back(alpha);
    
    if(alpha > 0.05) continue; // FILTER
    
    // This circle is better than previous if it's radius is closer to the desired one
    // TODO: This condition may not be the best one...
    double radiiDifference = theCircle->GetRadius() - testingCircle->GetRadius();
    
    if(radiiDifference < bestRadiiDifference){
      bestRadiiDifference = radiiDifference;
      bestCircle = make_unique<Circle>(testingCircle);
    }
  }
  
  return bestCircle;
}

unique_ptr<Circle> CircleProcessor::CopyCircleAddingRange(const unique_ptr<Circle> &circle,
                                                          const range<double> &phiRange)
{
  auto newCircle = make_unique<Circle>(circle);
  newCircle->phiRange = phiRange;
  return newCircle;
}

unique_ptr<Circle> CircleProcessor::GetParallelCircle(const unique_ptr<Circle> &circle,
                                                      const shared_ptr<Point> &point)
{
  // Calculate circle center and radius
  double p2_x = circle->GetLastPoint()->GetX();
  double p2_y = circle->GetLastPoint()->GetY();
  
  double r1_x = circle->GetCenter().GetX() - p2_x;
  double r1_y = circle->GetCenter().GetY() - p2_y;
  
  double delta_px = point->GetX() - p2_x;
  double delta_py = point->GetY() - p2_y;
  
  double A = (delta_px*delta_px + delta_py*delta_py)/(2*delta_py);
  double B = delta_px/delta_py;

  double r2_x = A/(r1_y/r1_x + B);
  double r2_y = A/(1 + r1_x/r1_y * B);
  
  double c2_x = r2_x + p2_x;
  double c2_y = r2_y + p2_y;
  
  
  auto newCircle = make_unique<Circle>(*point,
                                       make_unique<Point>(c2_x, c2_y, point->GetZ()),
                                       sqrt(r2_x*r2_x + r2_y*r2_y));
  
  vector<shared_ptr<Point>> circlePoints = {circle->GetLastPoint(), point };
  newCircle->SetPoints(circlePoints);
  
  
  // Calculate circle's range in phi
  
  double phiStart  = newCircle->GetPointAngle(0);
  double phiEnd    = newCircle->GetPointAngle(1);
  
  if(phiEnd > phiStart) phiEnd -= 2*TMath::Pi();
  
  double phiMin = min(phiEnd, phiStart)/TMath::Pi() * 180;
  double phiMax = max(phiEnd, phiStart)/TMath::Pi() * 180;
  
  auto r = range<double>(phiMin, phiMax);
  
  // add circle to the vector of segments
  newCircle = CopyCircleAddingRange(newCircle, r);
  
  return newCircle;
}


// Check the given point are perpendicular to x or y axis
bool CircleProcessor::IsPerpendicular(const Point &p1, const Point &p2,const Point &p3)
{
  double yDelta_a = p2.GetY() - p1.GetY();
  double xDelta_a = p2.GetX() - p1.GetX();
  double yDelta_b = p3.GetY() - p2.GetY();
  double xDelta_b = p3.GetX() - p2.GetX();
  
  if (fabs(xDelta_a) <= 0.000000001 && fabs(yDelta_b) <= 0.000000001) return false;
  
  if (fabs(yDelta_a) <= 0.0000001)        return true;
  else if (fabs(yDelta_b) <= 0.0000001)   return true;
  else if (fabs(xDelta_a)<= 0.000000001)  return true;
  else if (fabs(xDelta_b)<= 0.000000001)  return true;
  
  return false ;
}

unique_ptr<Circle> CircleProcessor::CalcCircle(const Point &p1, const Point &p2, const Point &p3)
{
  double yDelta_a= p2.GetY() - p1.GetY();
  double xDelta_a= p2.GetX() - p1.GetX();
  double yDelta_b= p3.GetY() - p2.GetY();
  double xDelta_b= p3.GetX() - p2.GetX();
  
  unique_ptr<Circle> circle = nullptr;
  
  double center_x;
  double center_y;
  double center_z;
  double radius;
  
  if (fabs(xDelta_a) <= 0.000000001 && fabs(yDelta_b) <= 0.000000001){
    center_x = 0.5*(p2.GetX() + p3.GetX());
    center_y = 0.5*(p1.GetY() + p2.GetY());
    center_z = p1.GetZ();
    
    Point center(center_x, center_y, center_z);
    radius = pointsProcessor->distanceXY(center, p1);
    
    circle = make_unique<Circle>(p1, center, radius);
    
    return circle;
  }
  
  double aSlope = yDelta_a/xDelta_a;
  double bSlope = yDelta_b/xDelta_b;
  if (fabs(aSlope-bSlope) <= 0.000000001) return circle;
  
  
  center_x = (aSlope*bSlope*(p1.GetY() - p3.GetY()) + bSlope*(p1.GetX() + p2 .GetX()) -
              aSlope*(p2.GetX()+p3.GetX()) )/(2* (bSlope-aSlope) );
  
  center_y  = -1*(center_x - (p1.GetX()+p2.GetX())/2)/aSlope +  (p1.GetY()+p2.GetY())/2;
  center_z = p1.GetZ();
  
  Point center(center_x, center_y, center_z);
  radius = pointsProcessor->distanceXY(center, p1);
  circle = make_unique<Circle>(p1, center, radius);
  
  return circle;
}

unique_ptr<Circle> CircleProcessor::GetCircleFromTriplet(const PointsTriplet &triplet)
{
  Point p1 = *triplet[0];
  Point p2 = *triplet[1];
  Point p3 = *triplet[2];
  
  unique_ptr<Circle> circle;
  
  if (!IsPerpendicular(p1, p2, p3) )         circle = CalcCircle(p1, p2, p3);
  else if (!IsPerpendicular(p1, p3, p2) )    circle = CalcCircle(p1, p3, p2);
  else if (!IsPerpendicular(p2, p1, p3) )    circle = CalcCircle(p2, p1, p3);
  else if (!IsPerpendicular(p2, p3, p1) )    circle = CalcCircle(p2, p3, p1);
  else if (!IsPerpendicular(p3, p2, p1) )    circle = CalcCircle(p3, p2, p1);
  else if (!IsPerpendicular(p3, p1, p2) )    circle = CalcCircle(p3, p1, p2);
  
  return circle;
  /*
  double x1 = triplet[0]->GetX();
  double y1 = triplet[0]->GetY();
  double x2 = triplet[1]->GetX();
  double y2 = triplet[1]->GetY();
  double x3 = triplet[2]->GetX();
  double y3 = triplet[2]->GetY();
  
  double x1_2 = x1*x1;
  double y1_2 = y1*y1;
  double x2_2 = x2*x2;
  double y2_2 = y2*y2;
  double x3_2 = x3*x3;
  double y3_2 = y3*y3;
  
  // Calculate center and radius of circle passing through 3 specified points
  double x0 = ((x1_2+y1_2)*(y2-y3)+(x2_2+y2_2)*(y3-y1)+(x3_2+y3_2)*(y1-y2))/(2*(x1*(y2-y3)-y1*(x2-x3)+x2*y3-x3*y2));
  double y0 = ((x1_2+y1_2)*(x3-x2)+(x2_2+y2_2)*(x1-x3)+(x3_2+y3_2)*(x2-x1))/(2*(x1*(y2-y3)-y1*(x2-x3)+x2*y3-x3*y2));
  double R  = sqrt(pow(x0 - x1, 2) + pow(y0 - y1, 2));
  
  auto center     = make_unique<Point>(x0, y0, 0.0);
  auto decayPoint = Point(*triplet[0]);
  
  auto circle = make_unique<Circle>(decayPoint, center, R);
  circle->SetPoints(triplet);
  circle->SetPz(triplet[1]->GetZ() - triplet[0]->GetZ());
  
  return circle;
   */
}
