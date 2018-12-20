//
//  Helix.cpp
//
//  Created by Jeremi Niedziela on 18/12/2018.
//

#include "Helix.hpp"

Helix::Helix()
{
  
}

Helix::Helix(double _R, double _c, double _x0, double _y0, double _z0,int _nCycles, double _thickness)
: R(_R), c(_c), x0(_x0), y0(_y0), z0(_z0), nCycles(_nCycles), thickness(_thickness)
{
//  tShift = acos(1/sqrt(pow(x0/y0,2)+1));
  tShift = acos(-x0/sqrt(x0*x0+y0*y0));
  
  tMin = -tShift;
  tMax = nCycles*2*TMath::Pi();
  tStep = 0.01;
}

Helix::Helix(double _R, double _c, Point p,int _nCycles, double _thickness)
: R(_R), c(_c), x0(p.x), y0(p.y), z0(p.z), nCycles(_nCycles), thickness(_thickness)
{
//  tShift = acos(1/sqrt(pow(x0/y0,2)+1));
  tShift = acos(-x0/sqrt(x0*x0+y0*y0));
  tMin = -tShift;
  tMax = nCycles*2*TMath::Pi();
  tStep = 0.01;
}

Helix::Helix(double _c, Circle circle, int _nCycles, double _thickness)
: c(_c), nCycles(_nCycles), thickness(_thickness)
{
  x0 = circle.x;
  y0 = circle.y;
  R = circle.R;

  // not sure why tShifts like that - to be understood
//  tShift = acos(1/sqrt(pow(x0/y0,2)+1));
  tShift = circle.tShift;
  z0 = circle.z + tShift*c;
//  tShift = circle.tShift;
}

void Helix::Print()
{
  cout<<"R:"<<R<<"\tc:"<<c<<"\toffset:("<<x0<<","<<y0<<","<<z0<<")\tnPoints:"<<nPoints<<"\tnPionPoints:"<<nPionPoints;
  cout<<"\tpz:"<<pz<<endl;
  
}

void Helix::Shift(int charge)
{
  int xSign=1, ySign=1;
  if(x0> 0 && y0> 0){xSign= 1; ySign=-1;}
  if(x0<=0 && y0> 0){xSign= 1; ySign= 1;}
  if(x0<=0 && y0<=0){xSign=-1; ySign= 1;}
  if(x0> 0 && y0<=0){xSign=-1; ySign=-1;}
  double xInit = x0, yInit=y0;
  
  // the charge may be inverted here... to be checked later
  x0 += charge * xSign * R/sqrt(pow(xInit/yInit,2)+1);
  y0 += charge * ySign * R/sqrt(pow(yInit/xInit,2)+1);
  z0 += tShift*c;
}

void Helix::ShiftByVector(Point v, int charge)
{
  v = Point(charge * -v.y,charge * v.x, v.z); // take a vector perpendicular to the pion's momentum vector
  pz = v.z;
  
  double vTransverseLength = sqrt(v.x*v.x+v.y*v.y);
  tShift = acos(-v.x/vTransverseLength);
  
  double scale = R/vTransverseLength;
  
  v.x *= scale;
  v.y *= scale;
  
  x0 += v.x;
  y0 += v.y;
  z0 += tShift*c;
}

vector<Point> Helix::GetPointsHittingSilicon()
{
  vector<Point> points;
  
  double dh = sqrt(x0*x0+y0*y0);
  double Rl, C, delta;
  double x1,y1,x2,y2,z1,z2,t1,t2;
  
  for(int iLayer=0;iLayer<4/*nLayers*/;iLayer++){
    Rl = layerR[iLayer];
    C = (Rl*Rl+dh*dh-R*R)/2.;
    
    delta = 4*y0*y0*C*C-4*dh*dh*(C*C-Rl*Rl*x0*x0);
    if(delta < 0) continue;
    
    y1 = (2*y0*C+sqrt(delta))/(2*dh*dh);
    y2 = (2*y0*C-sqrt(delta))/(2*dh*dh);
    
    x1 = (C-y1*y0)/x0;
    x2 = (C-y2*y0)/x0;
    
    t1 = atan2(y1-y0,x1-x0);
    t2 = atan2(y2-y0,x2-x0);
    
    for(int n=0;n<nCycles;n++){
      if(n>0 || t1 > -tShift){
        z1 = z0 + c*(t1+n*2*TMath::Pi());
        points.push_back(Point(x1, y1, z1));
      }
      if(n>0 || t2 > -tShift){
        z2 = z0 + c*(t2+n*2*TMath::Pi());
        points.push_back(Point(x2, y2, z2));
      }
    }
    
  }
  return points;
}
/*
vector<Point> Helix::GetPointsHittingSilicon(double threshold)
{
  vector<Point> points;
  
  double x,y,z;
  for(double t=tMin;t<tMax;t+=tStep){
    x = R*cos(t) + x0;
    y = R*sin(t) + y0;
    z = c*t      + z0;
    
    for(int iLayer=0;iLayer<4/*nLayers;iLayer++){
      if(fabs(sqrt(x*x+y*y)-layerR[iLayer]) < threshold){
        points.push_back(Point(x,y,z));
      }
    }
  }
  return points;
}
*/
Point Helix::GetClosestPoint(Point p)
{
  double t = atan2(p.y-y0, p.x-x0);
  
  double x = R*cos(t) + x0;
  double y = R*sin(t) + y0;
  double z = c*t      + z0;
  
  double absC = fabs(c);
  
  while(fabs(z-p.z) >= absC*2*TMath::Pi()){
    if(z < p.z) z += absC*2*TMath::Pi();
    else        z -= absC*2*TMath::Pi();
  }
  
  double currentDistanceZ = fabs(p.z-z);
  
  if(fabs(p.z-(z-absC*2*TMath::Pi())) < currentDistanceZ){
    z -= absC*2*TMath::Pi();
    currentDistanceZ = fabs(p.z-z);
  }
  else if(fabs(p.z-(z+absC*2*TMath::Pi())) < currentDistanceZ){
    z += absC*2*TMath::Pi();
  }
  
  return Point(x,y,z);
}

vector<Point>* Helix::GetMatchingPoints(vector<Point> &points)
{
  vector<Point> *result = new vector<Point>();
  
  for(Point p : points){
    Point q = GetClosestPoint(p);
    if(p.distance(q) < thickness) result->push_back(p);
  }
  return result;
}

void Helix::CountMatchingPoints(const vector<Point> &points)
{
  nPoints = 0;
  nPionPoints = 0;
  this->points.clear();
  
  for(Point p : points){
    Point q = GetClosestPoint(p);
    if(p.distance(q) < thickness){
      nPoints++;
      if(p.isPionHit) nPionPoints++;
      this->points.push_back(p);
    }
  }
}

vector<vector<Point>> Helix::SplitPointsIntoLines()
{
  vector<vector<Point>> pointsByLines;
  bool addedToExisting;
  
  for(Point p : points){
    addedToExisting = false;
    
    // loop over existing lines and check if this point belongs to one of them
    for(vector<Point> &line : pointsByLines){
      // if distance to this line is small enough, just add the point to this line and go to next point
      if(line[0].distanceXY(p) < 10){
        line.push_back(p);
        addedToExisting = true;
        break;
      }
    }
    if(addedToExisting) continue;
    
    // If the point was not added to any line, create a new line for it
    vector<Point> line;
    line.push_back(p);
    sort(line.begin(), line.end(),[](Point i1, Point i2){return i1.z < i2.z;});
    pointsByLines.push_back(line);
  }
  
  return pointsByLines;
}
