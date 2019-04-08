//
//  Helix.cpp
//
//  Created by Jeremi Niedziela on 18/12/2018.
//

#include "Helix.hpp"

Helix::Helix(const Point &_origin,
             const unique_ptr<Point> &_momentum,
             int _charge) :
vertex(make_unique<Point>(_origin)),
origin(_origin),
momentum(make_unique<Point>(*_momentum)),
charge(_charge),
pointsProcessor(make_unique<PointsProcessor>()),
circleProcessor(make_unique<CircleProcessor>())
{
  radius = GetRadiusInMagField(momentum->GetX(), momentum->GetY(), solenoidField);
  slope = radius * charge * momentum->GetVectorSlopeC();
  slopeAbs = fabs(slope);
  tStep = 0.01;
  
  // take a vector perpendicular to the pion's momentum vector
  Point v = Point(charge * momentum->GetY(),charge * -momentum->GetX(), 0.0);
  const double scale = radius/sqrt(pow(v.GetX(),2)+pow(v.GetY(),2));
  
  v.SetX(scale * v.GetX());
  v.SetY(scale * v.GetY());
  
  origin.SetX(origin.GetX() + v.GetX());
  origin.SetY(origin.GetY() + v.GetY());
 
  if(momentum->GetZ() > 0){
    if(charge > 0) tShift = TMath::Pi() - atan2(-v.GetX(), v.GetY());
    if(charge < 0) tShift = TMath::Pi() - atan2(-v.GetY(), v.GetX());
  }
  if(momentum->GetZ() < 0){
    if(charge > 0) tShift = TMath::Pi() - atan2(-v.GetY(), v.GetX());
    if(charge < 0) tShift = TMath::Pi() - atan2(-v.GetX(), v.GetY());
  }
  origin.SetZ(origin.GetZ() - fabs(tShift)*fabs(slope));
  
  tMax = GetNcycles()*2*TMath::Pi();
}

Helix::Helix(const Track &_track, const Point &p1, const Point &p2) :
origin(Point(0,0,0)),
momentum(make_unique<Point>(0,0,0)),
circleProcessor(make_unique<CircleProcessor>())
{
  int nLayers = _track.GetNtrackerLayers();
  
  Lmin = layerR[nLayers-1];
  Lmax = layerR[nLayers];
  
  // approximated position of the decay vertex
  vertex = make_unique<Point>((Lmin-Lmax)/2 * cos(_track.GetPhi()),
                              (Lmin-Lmax)/2 * sin(_track.GetPhi()),
                              (Lmin-Lmax)/2 / tan(_track.GetTheta()));
  
  PointsTriplet points = {make_shared<Point>(*vertex), make_shared<Point>(p1), make_shared<Point>(p2)};
  auto circle = circleProcessor->GetCircleFromTriplet(points);
  
  origin = circle->GetCenter();
  
  double t1 = circle->GetPointAngle(make_shared<Point>(p1));
  double t2 = circle->GetPointAngle(make_shared<Point>(p2));
  
  double b1 = (p2.GetZ()+p2.GetZerr() - Lmin/tan(_track.GetTheta())) / (t1*t2);
  double b2 = (p2.GetZ()-p2.GetZerr() - Lmin/tan(_track.GetTheta())) / (t1*t2);
  double b3 = (p2.GetZ()+p2.GetZerr() - Lmax/tan(_track.GetTheta())) / (t1*t2);
  double b4 = (p2.GetZ()-p2.GetZerr() - Lmax/tan(_track.GetTheta())) / (t1*t2);
  
  bmin = min(min(min(b1, b2), b3), b4);
  bmax = max(max(max(b1, b2), b3), b4);
  
  double s0_1 = (p1.GetZ()+p1.GetZerr() - (p2.GetZ()+p2.GetZerr()) + bmin*(t1*t1 - t2*t2)) / (t1 - t2);
  double s0_2 = (p1.GetZ()+p1.GetZerr() - (p2.GetZ()-p2.GetZerr()) + bmin*(t1*t1 - t2*t2)) / (t1 - t2);
  double s0_3 = (p1.GetZ()-p1.GetZerr() - (p2.GetZ()+p2.GetZerr()) + bmin*(t1*t1 - t2*t2)) / (t1 - t2);
  double s0_4 = (p1.GetZ()-p1.GetZerr() - (p2.GetZ()-p2.GetZerr()) + bmin*(t1*t1 - t2*t2)) / (t1 - t2);
  
  double s0_5 = (p1.GetZ()+p1.GetZerr() - (p2.GetZ()+p2.GetZerr()) + bmax*(t1*t1 - t2*t2)) / (t1 - t2);
  double s0_6 = (p1.GetZ()+p1.GetZerr() - (p2.GetZ()-p2.GetZerr()) + bmax*(t1*t1 - t2*t2)) / (t1 - t2);
  double s0_7 = (p1.GetZ()-p1.GetZerr() - (p2.GetZ()+p2.GetZerr()) + bmax*(t1*t1 - t2*t2)) / (t1 - t2);
  double s0_8 = (p1.GetZ()-p1.GetZerr() - (p2.GetZ()-p2.GetZerr()) + bmax*(t1*t1 - t2*t2)) / (t1 - t2);
  
  vector<double> s0 = { s0_1, s0_2, s0_3, s0_4, s0_5, s0_6, s0_7, s0_8 };
  
  s0min = *min_element(s0.begin(), s0.end());
  s0max = *max_element(s0.begin(), s0.end());
  
  amin = ((p1.GetX() - Lmin*cos(_track.GetPhi()))*(cos(t1)-cos(t2))) / (sin(t1)*cos(t2)*(t1-t2));
  amax = ((p1.GetX() - Lmax*cos(_track.GetPhi()))*(cos(t1)-cos(t2))) / (sin(t1)*cos(t2)*(t1-t2));
  
  if(amin > amax) swap(amin, amax);
  
  R0min = ((p1.GetX() - p2.GetX()) + amin*(t1*cos(t1)-t2*cos(t2))) / (cos(t1)-cos(t2));
  R0max = ((p1.GetX() - p2.GetX()) + amax*(t1*cos(t1)-t2*cos(t2))) / (cos(t1)-cos(t2));
  
  if(R0min > R0max) swap(R0min, R0max);
}

void Helix::Print()
{
  cout<<"\tVertex:("<<vertex->GetX()<<","<<vertex->GetY()<<","<<vertex->GetZ()<<")\n";
  cout<<"\tOrigin:("<<origin.GetX()<<","<<origin.GetY()<<","<<origin.GetZ()<<")\n";
  cout<<"\tMomentum:("<<momentum->GetX()<<","<<momentum->GetY()<<","<<momentum->GetZ()<<")\n";
  cout<<"\tCharge: "<<charge<<"\tR:"<<radius<<"\tc:"<<slope<<"\n";
  cout<<"\tnPoints:"<<points.size()<<"\tnPionPoints:"<<nPionPoints<<"\tnRegularPoints:"<<nRegularPoints<<"\n";
  cout<<"\tL:"<<Lmin<<" -- "<<Lmax<<endl;
  cout<<"\tb:"<<bmin<<" -- "<<bmax<<endl;
  cout<<"\ts0:"<<s0min<<" -- "<<s0max<<endl;
  cout<<"\ta:"<<amin<<" -- "<<amax<<endl;
  cout<<"\tR:"<<R0min<<" -- "<<R0max<<endl;
}

void Helix::SetPoints(const vector<shared_ptr<Point>> &_points)
{
  nPionPoints = 0;
  points.clear();
  
  for(auto &p : _points){
    Point q = GetClosestPoint(*p);
    if(pointsProcessor->distance(*p,q) < config->helixThickness){
      if(p->IsPionHit()) nPionPoints++;
      points.push_back(p);
    }
  }
}

double Helix::GetChi2()
{
  double chi2 = 0;
  for(auto &p : points){
    chi2 += pow(pointsProcessor->distance(*p, GetClosestPoint(*p)), 2);
  }
  return chi2 / points.size();
}

Point Helix::GetClosestPoint(Point p)
{
  int zSign = sgn(momentum->GetZ());
  
//  if(   (zSign > 0 && p.GetZ() < origin->GetZ())
//     || (zSign < 0 && p.GetZ() > origin->GetZ())){
//    return Point(origin->GetX(), origin->GetY(), origin->GetZ());
//  }
  double t = atan2((p.GetY()-origin.GetY()),(p.GetX()-origin.GetX()));
  
  if(charge < 0) t = TMath::Pi()/2. - t;
  
  double x = origin.GetX();
  double y = origin.GetY();
  double z = origin.GetZ() + slopeAbs*t;
  
  if(charge > 0){
    x += radius*cos(t);
    y += radius*sin(t);
  }
  else{
    x += radius*sin(t);
    y += radius*cos(t);
  }
  
  int nCycles = round(fabs(p.GetZ() - z) / (slopeAbs * 2 * TMath::Pi()));
  z += zSign * nCycles * slopeAbs * 2 * TMath::Pi();

  return Point(x,y,z);
}
