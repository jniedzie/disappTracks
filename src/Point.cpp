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

vector<vector<Point>> Point::SplitPointsIntoLines(vector<Point> points, double tolerance)
{
  vector<vector<Point>> pointsByLines;
  bool addedToExisting;
  
  for(Point p : points){
    addedToExisting = false;
    
    // loop over existing lines and check if this point belongs to one of them
    for(vector<Point> &line : pointsByLines){
      // if distance to this line is small enough, just add the point to this line and go to next point
      if(Point(line).distanceXY(p) < tolerance){
        line.push_back(p);
        addedToExisting = true;
        break;
      }
    }
    if(addedToExisting) continue;
    
    // If the point was not added to any line, create a new line for it
    vector<Point> line;
    line.push_back(p);
    pointsByLines.push_back(line);
  }
  
  return pointsByLines;
}

vector<Point> Point::GetRandomPoints(int nPoints)
{
  vector<Point> points;
  double phi, R;
  int layerIndex;
  
  for(int i=0;i<nPoints;i++){
    phi = RandDouble(0, 2*TMath::Pi());
    layerIndex = RandDouble(0, 4);
    R = layerR[layerIndex];
    Point p(R*cos(phi), R*sin(phi), RandDouble(-pixelBarrelZsize, pixelBarrelZsize));
    points.push_back(p);
  }
  return points;
}
