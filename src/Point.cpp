//
//  Point.cpp
//
//  Created by Jeremi Niedziela on 18/12/2018.
//

#include "Point.hpp"

Point::Point(double _x, double _y, double _z, double _val)
: x(_x), y(_y), z(_z), val(_val)
{
  
}

void Point::Print()
{
  cout<<"("<<x<<","<<y<<","<<z<<")"<<endl;
}

double Point::PerpendicularShift(double R,double c, int charge)
{
  int xSign=1, ySign=1;
  if(x> 0 && y> 0){xSign= 1; ySign=-1;}
  if(x<=0 && y> 0){xSign= 1; ySign= 1;}
  if(x<=0 && y<=0){xSign=-1; ySign= 1;}
  if(x> 0 && y<=0){xSign=-1; ySign=-1;}
  double tShift = acos(1/sqrt(pow(x/y,2)+1));
  double dx = x, dy=y;
  
  // the charge may be inverted here... to be checked later
  x += charge * xSign * R/sqrt(pow(dx/dy,2)+1);
  y += charge * ySign * R/sqrt(pow(dy/dx,2)+1);
  z += tShift*c;
  return tShift;
}

double Point::distance(Point p)
{
  return sqrt(pow(x-p.x,2)+pow(y-p.y,2)+pow(z-p.z,2));
}
