//
//  CircleProcessor.cpp
//
//  Created by Jeremi Niedziela on 26/03/2019.
//

#include "CircleProcessor.hpp"

CircleProcessor::CircleProcessor()
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

vector<unique_ptr<Circle>> CircleProcessor::BuildCirclesFromPoints(const vector<vector<shared_ptr<Point>>> &points);
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
