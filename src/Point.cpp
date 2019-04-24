//  Point.cpp
//
//  Created by Jeremi Niedziela on 18/12/2018.

#include "Point.hpp"

Point::Point() :
x(0),
y(0),
z(0),
errX(0),
errY(0),
errZ(0),
value(0),
t(0),
subDetName("")
{
  
}

Point::Point(double _x, double _y, double _z, double _value, string _subDetName,
             double _errX, double _errY, double _errZ, double _t) :
x(_x),
y(_y),
z(_z),
value(_value),
subDetName(_subDetName),
isPionHit(false),
errX(_errX),
errY(_errY),
errZ(_errZ),
t(_t)
{
  
}

Point::Point(const Point &p)
{
  x = p.x;
  y = p.y;
  z = p.z;
  value = p.value;
  isPionHit = p.isPionHit;
  errX = p.errX;
  errY = p.errY;
  errZ = p.errZ;
  subDetName = p.subDetName;
  t = p.t;
}

void Point::operator=(const Point &p)
{
  x = p.x;
  y = p.y;
  z = p.z;
  value = p.value;
  isPionHit = p.isPionHit;
  errX = p.errX;
  errY = p.errY;
  errZ = p.errZ;
  subDetName = p.subDetName;
  t = p.t;
}

Point::Point(vector<Point> points)
{
  x=0; y=0; z=0;
  for(Point p : points){
    x += p.GetX();
    y += p.GetY();
    z += p.GetZ();
    
    errX += p.GetXerr();
    errY += p.GetYerr();
    errZ += p.GetZerr();
    t += p.t;
  }
  x /= points.size();
  y /= points.size();
  z /= points.size();
  errX /= points.size();
  errY /= points.size();
  errZ /= points.size();
  t /= points.size();
}

void Point::Print() const
{
  cout<<"("<<x<<","<<y<<","<<z<<")";
  cout<<"\tVector length: "<<sqrt(x*x+y*y+z*z);
}

double Point::GetVectorSlopeC() const
{
  return tan(TMath::Pi()/2.-acos(z/sqrt(x*x+y*y+z*z)));
}


bool Point::operator==(const Point &p) const
{
  bool theSame = true;
  
  if(fabs(x - p.x) > 0.000001 ||
     fabs(y - p.y) > 0.000001 ||
     fabs(z - p.z) > 0.000001 ||
     fabs(value - p.value) > 0.000001) theSame = false;
  
  if(isPionHit != p.isPionHit ||
     subDetName != p.subDetName) theSame = false;
  
  return theSame;
}
