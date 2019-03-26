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
      circles.erase(circles.begin()+i);
      i--;
    }
  }
}

vector<unique_ptr<Circle>> CircleProcessor::BuildCirclesFromPoints(const TripletsVector &points)
{
  vector<unique_ptr<Circle>> circles;
  
  for(auto &p : points){
    double x1 = p[0]->GetX();
    double y1 = p[0]->GetY();
    double x2 = p[1]->GetX();
    double y2 = p[1]->GetY();
    double x3 = p[2]->GetX();
    double y3 = p[2]->GetY();
    
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
    auto decayPoint = make_unique<Point>(p[0]);
    
    auto circle = make_unique<Circle>(decayPoint, center, R);
    circle->SetPoints(p);
    
    circles.push_back(move(circle));
  }
  
  return circles;
}

unique_ptr<Circle> CircleProcessor::BuildCircleFromParams(const double *par,
                                                          const unique_ptr<Point> &vertex,
                                                          const shared_ptr<Track> &track)
{
  double L  = par[0];
  double px = par[1];
  double py = par[2];
  
  double x0 = L*cos(track->GetPhi())                          + 10*vertex->GetX();
  double y0 = L*sin(track->GetPhi())                          + 10*vertex->GetY();
  double z0 = L/sin(track->GetTheta())*cos(track->GetTheta()) + 10*vertex->GetZ();
  
  auto decayPoint  = make_unique<Point>(x0,y0,z0);
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
  double c1_x = theCircle->GetCenter()->GetX();
  double c1_y = theCircle->GetCenter()->GetY();
  double r1_x = c1_x - p_x;
  double r1_y = c1_y - p_y;
  double r1_mod = sqrt(r1_x*r1_x + r1_y*r1_y);
  
  for(auto &testingCircle : circles){
    // New track segment cannot have greater radius (within some tolerance)
    if(testingCircle->GetRadius() > 1.1*theCircle->GetRadius()) continue; // FILTER
    
    // The center of the new circle should be withing the previous circle
    double centerDifference = pointsProcessor->distanceXY(theCircle->GetCenter(), testingCircle->GetCenter());
    if(centerDifference > 1.1*theCircle->GetRadius()) continue; // FILTER
    
    // Here we calculate an angle between radius of the previous circle and radius of the testing circle
    // looking from the last point of previous circle. For perfectly matching circles that would be zero
    double r2_x   = testingCircle->GetCenter()->GetX() - p_x;
    double r2_y   = testingCircle->GetCenter()->GetY() - p_y;
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
