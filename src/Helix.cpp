//
//  Helix.cpp
//
//  Created by Jeremi Niedziela on 18/12/2018.
//

#include "Helix.hpp"

Helix::Helix(const unique_ptr<Point> &_origin,
             const unique_ptr<Point> &_momentum,
             int _charge) :
points(make_shared<vector<Point>>()),
vertex(make_unique<Point>(*_origin)),
origin(make_unique<Point>(*_origin)),
momentum(make_unique<Point>(*_momentum)),
charge(_charge),
pointsProcessor(make_unique<PointsProcessor>())
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
  
  origin->SetX(origin->GetX() + v.GetX());
  origin->SetY(origin->GetY() + v.GetY());
 
  if(momentum->GetZ() > 0){
    if(charge > 0) tShift = TMath::Pi() - atan2(-v.GetX(), v.GetY());
    if(charge < 0) tShift = TMath::Pi() - atan2(-v.GetY(), v.GetX());
  }
  if(momentum->GetZ() < 0){
    if(charge > 0) tShift = TMath::Pi() - atan2(-v.GetY(), v.GetX());
    if(charge < 0) tShift = TMath::Pi() - atan2(-v.GetX(), v.GetY());
  }
  origin->SetZ(origin->GetZ() - fabs(tShift)*fabs(slope));
  
  tMax = GetNcycles()*2*TMath::Pi();
}

void Helix::Print()
{
  cout<<"\tVertex:("<<vertex->GetX()<<","<<vertex->GetY()<<","<<vertex->GetZ()<<")\n";
  cout<<"\tOrigin:("<<origin->GetX()<<","<<origin->GetY()<<","<<origin->GetZ()<<")\n";
  cout<<"\tMomentum:("<<momentum->GetX()<<","<<momentum->GetY()<<","<<momentum->GetZ()<<")\n";
  cout<<"\tCharge: "<<charge<<"\tR:"<<radius<<"\tc:"<<slope<<"\n";
  cout<<"\tnPoints:"<<points->size()<<"\tnPionPoints:"<<nPionPoints<<"\tnRegularPoints:"<<nRegularPoints<<"\n";
}

void Helix::SetPoints(const shared_ptr<vector<Point>> _points)
{
  nPionPoints = 0;
  points->clear();
  
  for(Point p : *_points){
    Point q = GetClosestPoint(p);
    if(pointsProcessor->distance(p,q) < config->helixThickness){
      if(p.IsPionHit()) nPionPoints++;
      points->push_back(p);
    }
  }
}

double Helix::GetChi2()
{
  double chi2 = 0;
  for(Point p : *points){
    chi2 += pow(pointsProcessor->distance(p, GetClosestPoint(p)), 2);
  }
  return chi2 / points->size();
}

Point Helix::GetClosestPoint(Point p)
{
  int zSign = sgn(momentum->GetZ());
  
//  if(   (zSign > 0 && p.GetZ() < origin->GetZ())
//     || (zSign < 0 && p.GetZ() > origin->GetZ())){
//    return Point(origin->GetX(), origin->GetY(), origin->GetZ());
//  }
  double t = atan2((p.GetY()-origin->GetY()),(p.GetX()-origin->GetX()));
  
  if(charge < 0) t = TMath::Pi()/2. - t;
  
  double x = origin->GetX();
  double y = origin->GetY();
  double z = origin->GetZ() + slopeAbs*t;
  
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
