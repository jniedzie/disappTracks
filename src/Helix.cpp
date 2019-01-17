//
//  Helix.cpp
//
//  Created by Jeremi Niedziela on 18/12/2018.
//

#include "Helix.hpp"

Helix::Helix(const unique_ptr<Point> &_origin,
             const unique_ptr<Point> &_momentum,
             int _charge, shared_ptr<FitterConfig> _config) :
origin(make_unique<Point>(*_origin)),
momentum(make_unique<Point>(*_momentum)),
charge(_charge),
config(_config),
pointsProcessor(make_unique<PointsProcessor>())
{
  radius = GetRadiusInMagField(momentum->GetX(), momentum->GetY(), solenoidField);
  slope = charge*momentum->GetVectorSlopeC();
  slopeAbs = fabs(slope);
  tStep = 0.01;
  
  // take a vector perpendicular to the pion's momentum vector
  const Point v = Point(charge * -momentum->GetY(),charge * momentum->GetX(), 0.0);
  const double scale = radius/sqrt(pow(v.GetX(),2)+pow(v.GetY(),2));
  
  origin->SetX(origin->GetX() + scale*v.GetX());
  origin->SetY(origin->GetY() + scale*v.GetY());
 
  tShift = -atan2(charge * v.GetY(), -charge * v.GetX());
  
  origin->SetZ(origin->GetZ() - tShift*slopeAbs);
  tMax = GetNcycles()*2*TMath::Pi();
}

Helix::Helix(const unique_ptr<Circle> &_circle, double _pz, int _charge) :
charge(_charge),
pointsProcessor(make_unique<PointsProcessor>())
{
  tStep = 0.01;
  config    = _circle->GetConfig();
  radius    = _circle->GetRadius();
  tShift    = _circle->GetToffset();
  
  momentum  = _circle->GetMomentum();
  momentum->SetZ(_pz);
  
  slope = charge * momentum->GetVectorSlopeC();
  slopeAbs = fabs(slope);
  
  origin = make_unique<Point>(_circle->GetCenter()->GetX(),
                              _circle->GetCenter()->GetY(),
                              _circle->GetCenter()->GetZ() - tShift*slopeAbs);
  
  tMax = GetNcycles()*2*TMath::Pi();
}

void Helix::Print()
{
  cout<<"\tOrigin:("<<origin->GetX()<<","<<origin->GetY()<<","<<origin->GetZ()<<")\t";
  cout<<"Momentum:("<<momentum->GetX()<<","<<momentum->GetY()<<","<<momentum->GetZ()<<")\n";
cout<<"\tR:"<<radius<<"\tc:"<<slope<<"\tnPoints:"<<points.size()<<"\tnPionPoints:"<<nPionPoints<<"\tnRegularPoints:"<<nRegularPoints<<"\n";
}

vector<Point> Helix::GetPointsHittingSilicon()
{
  vector<Point> points;
  
  double dh = sqrt(pow(origin->GetX(),2)+pow(origin->GetY(),2));
  double Rl, C, delta;
  double x1,y1,x2,y2,z1,z2,t1,t2;
  
  for(int iLayer=0;iLayer<nPixelLayers;iLayer++){
    Rl = layerR[iLayer];
    C = (Rl*Rl+dh*dh-radius*radius)/2.;
    
    delta = 4*origin->GetY()*origin->GetY()*C*C-4*dh*dh*(C*C-Rl*Rl*origin->GetX()*origin->GetX());
    if(delta < 0) continue;
    
    y1 = (2*origin->GetY()*C+sqrt(delta))/(2*dh*dh);
    y2 = (2*origin->GetY()*C-sqrt(delta))/(2*dh*dh);
    
    x1 = (C-y1*origin->GetY())/origin->GetX();
    x2 = (C-y2*origin->GetY())/origin->GetX();
    
    t1 = atan2((y1-origin->GetY()),(x1-origin->GetX()));
    t2 = atan2(y2-origin->GetY(),x2-origin->GetX());
    
    if(charge < 0){
      t1 = TMath::Pi()/2. - t1;
      t2 = TMath::Pi()/2. - t2;
    }
    
    double nCycles = fabs(GetNcycles());
    int signZ = sgn(momentum->GetZ());
    
    for(double n=0;n<nCycles;n+=1){
      if(n>0
         || (signZ > 0 && t1 > tShift)
         || (signZ < 0 && t1 < tShift)
         ){
        z1 = origin->GetZ() + slopeAbs*(t1 + signZ*n*2*TMath::Pi());
        
        points.push_back(Point(x1, y1, z1));
      }
      if(n>0
         || (signZ > 0 && t2 > tShift)
         || (signZ < 0 && t2 < tShift)
         ){
        z2 = origin->GetZ() + slopeAbs*(t2 + signZ*n*2*TMath::Pi());
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
    if(pointsProcessor->distance(p,q) < config->GetHelixThickness()){
      if(p.IsPionHit()) nPionPoints++;
      points.push_back(p);
    }
  }
}

double Helix::GetChi2()
{
  double chi2 = 0;
  for(Point p : points){
    chi2 += pow(pointsProcessor->distance(p, GetClosestPoint(p)), 2);
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

void Helix::CalculateNregularPoints(int limit)
{
  vector<vector<Point>> pointsByLine = pointsProcessor->SplitPointsIntoLines(points, config->GetLinesToleranceForRegularity());
  vector<double> possibleDistances;
  set<double> possibleDistancesSet;
  nRegularPoints = 0;
  int nPointsForDistance;
  double zRegularityTolerance = config->GetZregularityTolerance();
  bool found;
  double testingDistance;
  
  for(auto line : pointsByLine){
    int iPoint;
    for(iPoint=0; iPoint < (int)line.size()-1; iPoint++){
      testingDistance = pointsProcessor->distance(line[iPoint], line[iPoint+1]);
      found = false;
      for(double dd : possibleDistances){
        if(fabs(testingDistance-dd) < zRegularityTolerance){found = true;break;}
      }
      if(found) continue;
      
      possibleDistances.push_back(testingDistance);
      nPointsForDistance = 0;
      
      for(auto line2 : pointsByLine){
        for(int i=0;i<(int)line2.size();i++){
          if(std::abs(pointsProcessor->distance(line2[0], line2[i])-i*testingDistance) < zRegularityTolerance)  nPointsForDistance++;
        }
      }
      if(nPointsForDistance > nRegularPoints){
        nRegularPoints = nPointsForDistance;
        if(nRegularPoints > limit) return;
      }
    }
  }
  
}

vector<int> Helix::AreHelicesIdentical(const unique_ptr<Helix> &h1, const unique_ptr<Helix> &h2)
{
  vector<int> reasons;
  shared_ptr<FitterConfig> config = h1->config;
  
  if(fabs(h1->GetOrigin()->GetX() - h2->GetOrigin()->GetX()) > config->GetToleranceX()) reasons.push_back(1);
  if(fabs(h1->GetOrigin()->GetY() - h2->GetOrigin()->GetY()) > config->GetToleranceY()) reasons.push_back(2);
  if(fabs(h1->GetOrigin()->GetZ() - h2->GetOrigin()->GetZ()) > config->GetToleranceZ()) reasons.push_back(3);
  if(fabs(h1->GetMomentum()->GetX() - h2->GetMomentum()->GetX()) > config->GetTolerancePx()) reasons.push_back(4);
  if(fabs(h1->GetMomentum()->GetY() - h2->GetMomentum()->GetY()) > config->GetTolerancePy()) reasons.push_back(5);
  if(fabs(h1->GetMomentum()->GetZ() - h2->GetMomentum()->GetZ()) > config->GetTolerancePz()) reasons.push_back(6);
  if(h1->GetCharge() != h2->GetCharge()) reasons.push_back(7);
  
  return reasons;
}
