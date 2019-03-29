//
//  Circle.cpp
//
//  Created by Jeremi Niedziela on 18/12/2018.
//

#include "Circle.hpp"

Circle::Circle(const unique_ptr<Point> &_decayPoint,
               const unique_ptr<Point> &_momentum,
               const range<double> &_phiRange) :
decayPoint(make_unique<Point>(*_decayPoint)),
momentum(make_unique<Point>(*_momentum)),
phiRange(_phiRange)
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

Circle::Circle(const unique_ptr<Point> &_decayPoint,
               const unique_ptr<Point> &_center,
               double _radius,
               const range<double> &_phiRange) :
decayPoint(make_unique<Point>(*_decayPoint)),
center(make_unique<Point>(*_center)),
radius(_radius),
phiRange(_phiRange)
{
  double pt = GetPtInMagField(radius, solenoidField);
  
  tShift = TMath::Pi()/2. + GetPointAngle(make_shared<Point>(*decayPoint));
  
  momentum = make_unique<Point>(pt * sin(tShift), pt * cos(tShift), 0.0);
}

Circle::Circle(const unique_ptr<Circle> &c)
{
  decayPoint = make_unique<Point>(*c->decayPoint);
  center     = make_unique<Point>(*c->center);
  momentum   = make_unique<Point>(*c->momentum);
  
  for(auto p : c->points){ points.push_back(p);}
  radius   = c->radius;
  tShift   = c->tShift;
  phiRange = c->phiRange;
}

void Circle::Print()
{
  cout<<"Circle center ("<<center->GetX()<<","<<center->GetY()<<")\tR:"<<radius<<endl;
  cout<<"\tMomentum:";momentum->Print();cout<<endl;
  cout<<"\tRange:"<<phiRange.GetMin()<<" -- "<<phiRange.GetMax()<<endl;
}

double Circle::GetDistanceToPoint(Point p)
{
  double t = atan2(p.GetY()-center->GetY(), p.GetX()-center->GetX());
  double x = radius*cos(t) + center->GetX();
  double y = radius*sin(t) + center->GetY();
  return sqrt(pow(p.GetX()-x,2)+pow(p.GetY()-y,2));
}

void Circle::SetPoints(const vector<shared_ptr<Point>> &_points)
{
  points.clear();
  for(auto &p : _points){
    if(GetDistanceToPoint(*p) < config->circleThickness) points.push_back(p);
  }
}

TArc* Circle::GetArc()
{
  return new TArc(center->GetX(),center->GetY(),radius);
}

double Circle::GetPointAngle(uint i)
{
  return TMath::Pi()/2. - atan2( (points[i]->GetX() - center->GetX()),
                                 (points[i]->GetY() - center->GetY()));
}

double Circle::GetPointAngle(const shared_ptr<Point> &point)
{
  return TMath::Pi()/2. - atan2(  (point->GetX() - center->GetX()),
                                  (point->GetY() - center->GetY()));
}
