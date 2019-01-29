//
//  Circle.cpp
//
//  Created by Jeremi Niedziela on 18/12/2018.
//

#include "Circle.hpp"

Circle::Circle(const unique_ptr<Point> &_decayPoint,
               const unique_ptr<Point> &_momentum,
               shared_ptr<FitterConfig> _config) :
decayPoint(make_unique<Point>(*_decayPoint)),
momentum(make_unique<Point>(*_momentum)),
config(_config)
{
  radius = GetRadiusInMagField(momentum->GetX(), momentum->GetY(), solenoidField);
  
  // take a vector perpendicular to the pion's momentum vector
  Point v = Point(-momentum->GetY(), momentum->GetX(), 0.0);
  double scale = radius/sqrt(pow(v.GetX(),2)+pow(v.GetY(),2));
  
  center = make_unique<Point>(*decayPoint);
  center->SetX(center->GetX() + scale*v.GetX());
  center->SetY(center->GetY() + scale*v.GetY());
  
  tShift = -atan2(v.GetY(), -v.GetX());
}

void Circle::Print()
{
  cout<<"Circle center ("<<center->GetX()<<","<<center->GetY()<<")\tR:"<<radius<<endl;
}

double Circle::GetDistanceToPoint(Point p)
{
  double t = atan2(p.GetY()-center->GetY(), p.GetX()-center->GetX());
  double x = radius*cos(t) + center->GetX();
  double y = radius*sin(t) + center->GetY();
  return sqrt(pow(p.GetX()-x,2)+pow(p.GetY()-y,2));
}

void Circle::SetPoints(const shared_ptr<vector<Point>> _points)
{
  points.clear();
  for(auto p : *_points){
    if(GetDistanceToPoint(p) < config->circleThickness) points.push_back(p);
  }
}

TArc* Circle::GetArc()
{
  return new TArc(center->GetX(),center->GetY(),radius);
}

void Circle::RemoveSimilarCircles(vector<unique_ptr<Circle>> &circles)
{
  if(circles.size() == 0) return;
  
  shared_ptr<FitterConfig> config = circles[0]->GetConfig();
  
  sort(circles.begin(), circles.end(),
       [](const auto &c1, const auto &c2) -> bool {return c1->GetRadius() < c2->GetRadius();});
  
  for(uint i=0; i<circles.size()-1; i++){
    if(fabs(circles[i]->GetRadius() - circles[i+1]->GetRadius()) < config->circleThickness){
      circles.erase(circles.begin()+i);
      i--;
    }
  }
}
