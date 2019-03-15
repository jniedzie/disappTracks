//
//  Point.cpp
//
//  Created by Jeremi Niedziela on 18/12/2018.
//

#include "Point.hpp"

Point::Point(double _x, double _y, double _z, double _value, string _subDetName) :
x(_x),
y(_y),
z(_z),
value(_value),
subDetName(_subDetName),
isPionHit(false)
{
  
}

Point::Point(const unique_ptr<Point> &p)
{
  x = p->x;
  y = p->y;
  z = p->z;
  value = p->value;
  isPionHit = p->isPionHit;
}

Point::Point(vector<Point> points)
{
  x=0; y=0; z=0;
  for(Point p : points){
    x += p.GetX();
    y += p.GetY();
    z += p.GetZ();
  }
  x /= points.size();
  y /= points.size();
  z /= points.size();
}

void Point::Print() const
{
  cout<<"("<<x<<","<<y<<","<<z<<")"<<endl;
}

double Point::GetVectorSlopeC() const
{
  return tan(TMath::Pi()/2.-acos(z/sqrt(x*x+y*y+z*z)));
}


