//
//  Helix.cpp
//
//  Created by Jeremi Niedziela on 18/12/2018.
//

#include "Helix.hpp"

Helix::Helix(Point *_origin, shared_ptr<Point> _momentum, int _charge,
             int _nCycles, double _thickness, double _zRegularityTolerance) :
origin(_origin),
momentum(_momentum),
charge(_charge),
nCycles(_nCycles),
thickness(_thickness),
zRegularityTolerance(_zRegularityTolerance)
{
  radius = GetRadiusInMagField(momentum->x, momentum->y, solenoidField);
  slope = momentum->GetVectorSlopeC();
  
  tShift = acos(-origin->x/sqrt(pow(origin->x,2)+pow(origin->y,2)));
  tMin = -tShift;
  tMax = nCycles*2*TMath::Pi();
  tStep = 0.01;
  
  ShiftByVector();
}

Helix::Helix(double _slope, Circle circle,
             int _nCycles, double _thickness, double _zRegularityTolerance) :
slope(_slope),
nCycles(_nCycles),
thickness(_thickness),
zRegularityTolerance(_zRegularityTolerance)
{
  radius    = circle.R;
  tShift    = circle.tShift;
  origin    = new Point(circle.x, circle.y, circle.z + tShift*slope);
  momentum  = shared_ptr<Point>(new Point(circle.px, circle.py, 0.0));
}

void Helix::Print()
{
  cout<<"R:"<<radius<<"\tc:"<<slope<<"\toffset:("<<origin->x<<","<<origin->y<<","<<origin->z<<")\tnPoints:"<<nPoints<<"\tnPionPoints:"<<nPionPoints;
  cout<<"\tpz:"<<momentum->z<<endl;
  
}

void Helix::ShiftByVector()
{
  // take a vector perpendicular to the pion's momentum vector
  Point v = Point(charge * -momentum->y,charge * momentum->x, momentum->z);
  
  double vTransverseLength = sqrt(v.x*v.x+v.y*v.y);
  tShift = acos(-v.x/vTransverseLength);
  
  double scale = radius/vTransverseLength;
  
  v.x *= scale;
  v.y *= scale;
  
  origin->x += v.x;
  origin->y += v.y;
  origin->z += tShift*slope;
}

vector<Point> Helix::GetPointsHittingSilicon()
{
  vector<Point> points;
  
  double dh = sqrt(origin->x*origin->x+origin->y*origin->y);
  double Rl, C, delta;
  double x1,y1,x2,y2,z1,z2,t1,t2;
  
  for(int iLayer=0;iLayer<4/*nLayers*/;iLayer++){
    Rl = layerR[iLayer];
    C = (Rl*Rl+dh*dh-radius*radius)/2.;
    
    delta = 4*origin->y*origin->y*C*C-4*dh*dh*(C*C-Rl*Rl*origin->x*origin->x);
    if(delta < 0) continue;
    
    y1 = (2*origin->y*C+sqrt(delta))/(2*dh*dh);
    y2 = (2*origin->y*C-sqrt(delta))/(2*dh*dh);
    
    x1 = (C-y1*origin->y)/origin->x;
    x2 = (C-y2*origin->y)/origin->x;
    
    t1 = atan2(y1-origin->y,x1-origin->x);
    t2 = atan2(y2-origin->y,x2-origin->x);
    
    for(int n=0;n<nCycles;n++){
      if(n>0 || t1 > -tShift){
        z1 = origin->z + slope*(t1+n*2*TMath::Pi());
        points.push_back(Point(x1, y1, z1));
      }
      if(n>0 || t2 > -tShift){
        z2 = origin->z + slope*(t2+n*2*TMath::Pi());
        points.push_back(Point(x2, y2, z2));
      }
    }
    
  }
  return points;
}



void Helix::SetPoints(const vector<Point> &_points)
{
  nPoints = 0;
  nPionPoints = 0;
  points.clear();
  
  for(Point p : _points){
    Point q = GetClosestPoint(p);
    if(p.distance(q) < thickness){
      nPoints++;
      if(p.isPionHit) nPionPoints++;
      points.push_back(p);
    }
  }
  
  CalculateNregularPoints();
}

double Helix::GetChi2()
{
  double chi2 = 0;
  for(Point p : points){
    Point q = GetClosestPoint(p);
    chi2 += pow(p.distance(q), 2);
  }
  return chi2 / nPoints;
}

Point Helix::GetClosestPoint(Point p)
{
  double t = atan2(p.y-origin->y, p.x-origin->x);
  
  double x = radius*cos(t) + origin->x;
  double y = radius*sin(t) + origin->y;
  double z = slope*t      + origin->z;
  
  double absC = fabs(slope);
  
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

void Helix::CalculateNregularPoints()
{
  vector<vector<Point>> pointsByLine = SplitPointsIntoLines();
  vector<double> possibleDistances;
  
  for(auto line : pointsByLine){
    bool first=true;
    for(auto p : line){
      if(first){first=false;continue;}
      possibleDistances.push_back(line[0].distance(p));
    }
  }
  
  nRegularPoints = 0;
  
  for(double testingDistance : possibleDistances){
    int nPointsForDistance = 0;
    
    for(auto line : pointsByLine){
      for(int n=0;n<20;n++){ // should just go till the edge of the pixel barrel or the most distant hit
        for(auto q : line){
          if(fabs(q.distance(line[0])-n*testingDistance) < zRegularityTolerance){
            nPointsForDistance++;
          }
        }
      }
    }
    
    if(nPointsForDistance > nRegularPoints) nRegularPoints = nPointsForDistance;
  }
}
