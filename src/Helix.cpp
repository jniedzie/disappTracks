//
//  Helix.cpp
//
//  Created by Jeremi Niedziela on 18/12/2018.
//

#include "Helix.hpp"

Helix::Helix(const unique_ptr<Point> &_origin,
             const unique_ptr<Point> &_momentum,
             int _charge,
             int _nCycles, double _thickness, double _zRegularityTolerance) :
origin(make_unique<Point>(*_origin)),
momentum(make_unique<Point>(*_momentum)),
charge(_charge),
nCycles(_nCycles),
thickness(_thickness),
zRegularityTolerance(_zRegularityTolerance)
{
  radius = GetRadiusInMagField(momentum->GetX(), momentum->GetY(), solenoidField);
  slope = momentum->GetVectorSlopeC();
  
  tMax = nCycles*2*TMath::Pi();
  tStep = 0.01;
  
  // take a vector perpendicular to the pion's momentum vector
  Point v = Point(charge * -momentum->GetY(),charge * momentum->GetX(), momentum->GetZ());
  
  double vTransverseLength = sqrt(pow(v.GetX(),2)+pow(v.GetY(),2));
  tShift = acos(-v.GetX()/vTransverseLength);
  double scale = radius/vTransverseLength;
  
  origin->SetX(origin->GetX() + scale*v.GetX());
  origin->SetY(origin->GetY() + scale*v.GetY());
  origin->SetZ(origin->GetZ() + tShift*slope);
}

Helix::Helix(double _slope, const unique_ptr<Circle> &_circle,
             int _nCycles, double _thickness, double _zRegularityTolerance) :
slope(_slope),
nCycles(_nCycles),
thickness(_thickness),
zRegularityTolerance(_zRegularityTolerance)
{
  radius    = _circle->GetRadius();
  tShift    = _circle->GetToffset();
  origin    = make_unique<Point>(_circle->GetCenter()->GetX(),
                                 _circle->GetCenter()->GetY(),
                                 _circle->GetCenter()->GetZ() + tShift*slope);
  momentum  = _circle->GetMomentum();
}

void Helix::Print()
{
  cout<<"R:"<<radius<<"\tc:"<<slope<<"\toffset:("<<origin->GetX()<<","<<origin->GetY()<<","<<origin->GetZ()<<")\tnPoints:"<<points.size()<<"\tnPionPoints:"<<nPionPoints;
  cout<<"\tpz:"<<momentum->GetZ()<<endl;
  
}

vector<Point> Helix::GetPointsHittingSilicon()
{
  vector<Point> points;
  
  double dh = sqrt(pow(origin->GetX(),2)+pow(origin->GetY(),2));
  double Rl, C, delta;
  double x1,y1,x2,y2,z1,z2,t1,t2;
  
  for(int iLayer=0;iLayer<4/*nLayers*/;iLayer++){
    Rl = layerR[iLayer];
    C = (Rl*Rl+dh*dh-radius*radius)/2.;
    
    delta = 4*origin->GetY()*origin->GetY()*C*C-4*dh*dh*(C*C-Rl*Rl*origin->GetX()*origin->GetX());
    if(delta < 0) continue;
    
    y1 = (2*origin->GetY()*C+sqrt(delta))/(2*dh*dh);
    y2 = (2*origin->GetY()*C-sqrt(delta))/(2*dh*dh);
    
    x1 = (C-y1*origin->GetY())/origin->GetX();
    x2 = (C-y2*origin->GetY())/origin->GetX();
    
    t1 = atan2(y1-origin->GetY(),x1-origin->GetX());
    t2 = atan2(y2-origin->GetY(),x2-origin->GetX());
    
    for(int n=0;n<nCycles;n++){
      if(n>0 || t1 > -tShift){
        z1 = origin->GetZ() + slope*(t1+n*2*TMath::Pi());
        points.push_back(Point(x1, y1, z1));
      }
      if(n>0 || t2 > -tShift){
        z2 = origin->GetZ() + slope*(t2+n*2*TMath::Pi());
        points.push_back(Point(x2, y2, z2));
      }
    }
    
  }
  return points;
}

void Helix::SetPoints(const vector<Point> &_points)
{
  nPionPoints = 0;
  points.clear();
  
  for(Point p : _points){
    Point q = GetClosestPoint(p);
    if(p.distance(q) < thickness){
      if(p.IsPionHit()) nPionPoints++;
      points.push_back(p);
    }
  }
  CalculateNregularPoints();
}

double Helix::GetChi2()
{
  double chi2 = 0;
  for(Point p : points){
    chi2 += pow(p.distance(GetClosestPoint(p)), 2);
  }
  return chi2 / points.size();
}

Point Helix::GetClosestPoint(Point p)
{
  double t = atan2(p.GetY()-origin->GetY(), p.GetX()-origin->GetX());
  
  double x = radius*cos(t) + origin->GetX();
  double y = radius*sin(t) + origin->GetY();
  double z = slope*t       + origin->GetZ();
  
  double absC = fabs(slope);
  
  while(fabs(z-p.GetZ()) >= absC*2*TMath::Pi()){
    if(z < p.GetZ())  z += absC*2*TMath::Pi();
    else              z -= absC*2*TMath::Pi();
  }
  
  double currentDistanceZ = fabs(p.GetZ()-z);
  
  if(fabs(p.GetZ()-(z-absC*2*TMath::Pi())) < currentDistanceZ){
    z -= absC*2*TMath::Pi();
    currentDistanceZ = fabs(p.GetZ()-z);
  }
  else if(fabs(p.GetZ()-(z+absC*2*TMath::Pi())) < currentDistanceZ){
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
    sort(line.begin(), line.end(),[](Point p1, Point p2){return p1.GetZ() < p2.GetZ();});
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
