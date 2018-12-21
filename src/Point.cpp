//
//  Point.cpp
//
//  Created by Jeremi Niedziela on 18/12/2018.
//

#include "Point.hpp"

Point::Point(double _x, double _y, double _z, double _value) :
x(_x),
y(_y),
z(_z),
value(_value),
isPionHit(false)
{
  
}

void Point::Print()
{
  cout<<"("<<x<<","<<y<<","<<z<<")"<<endl;
}

double Point::distance(Point p)
{
  return sqrt(pow(x-p.x,2)+pow(y-p.y,2)+pow(z-p.z,2));
}

double Point::distanceXY(Point p)
{
  return sqrt(pow(x-p.x,2)+pow(y-p.y,2));
}

double Point::GetVectorSlopeC()
{
  return tan(TMath::Pi()/2.-acos(z/sqrt(x*x+y*y+z*z)));
}
