//
//  Circle.cpp
//
//  Created by Jeremi Niedziela on 18/12/2018.
//

#include "Circle.hpp"

Circle::Circle(const unique_ptr<Point> &_decayPoint,
               const unique_ptr<Point> &_momentum) :
decayPoint(make_unique<Point>(*_decayPoint)),
momentum(make_unique<Point>(*_momentum))
{
  radius = GetRadiusInMagField(momentum->GetX(), momentum->GetY(), solenoidField);
  
  // take a vector perpendicular to the pion's momentum vector
  Point v = Point(momentum->GetY(), -momentum->GetX(), 0.0);
  double scale = radius/sqrt(pow(v.GetX(),2)+pow(v.GetY(),2));
  
  center = make_unique<Point>(*decayPoint);
  center->SetX(center->GetX() + scale*v.GetX());
  center->SetY(center->GetY() + scale*v.GetY());
  
  tShift = -atan2(v.GetY(), -v.GetX());
}

Circle::Circle(const unique_ptr<Circle> &c)
{
  decayPoint = make_unique<Point>(c->decayPoint);
  center = make_unique<Point>(c->center);
  momentum = make_unique<Point>(c->momentum);
  
  for(auto p : c->points){ points.push_back(p);}
  radius = c->radius;
  tShift = c->tShift;
}

void Circle::Print()
{
  cout<<"Circle center ("<<center->GetX()<<","<<center->GetY()<<")\tR:"<<radius<<endl;
  cout<<"Momentum:";momentum->Print();cout<<endl;
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
  
  sort(circles.begin(), circles.end(),
       [](const auto &c1, const auto &c2) -> bool {return c1->GetRadius() < c2->GetRadius();});
  
  for(uint i=0; i<circles.size()-1; i++){
    if(fabs(circles[i]->GetRadius() - circles[i+1]->GetRadius()) < config->circleThickness){
      circles.erase(circles.begin()+i);
      i--;
    }
  }
}

double Circle::GetPointAngle(double x, double y)
{
  return TMath::Pi()/2. + atan2(-(x-center->GetX()),
                                 (y- center->GetY()));
}
