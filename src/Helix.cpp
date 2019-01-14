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
config(_config)
{
  radius = GetRadiusInMagField(momentum->GetX(), momentum->GetY(), solenoidField);
  slope = momentum->GetVectorSlopeC();
  
  tMax = GetNcycles()*2*TMath::Pi();
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

Helix::Helix(double _slope, const unique_ptr<Circle> &_circle, shared_ptr<FitterConfig> _config) :
slope(_slope),
config(_config)
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
  cout<<"R:"<<radius<<"\tc:"<<slope<<"\toffset:("<<origin->GetX()<<","<<origin->GetY()<<","<<origin->GetZ()<<")\t";
  cout<<"nPoints:"<<points.size()<<"\tnPionPoints:"<<nPionPoints<<"\tnRegularPoints:"<<nRegularPoints<<"\t";
  cout<<"pz:"<<momentum->GetZ()<<endl;
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
    
    for(int n=0;n<GetNcycles();n++){
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
    if(p.distance(q) < config->GetHelixThickness()){
      if(p.IsPionHit()) nPionPoints++;
      points.push_back(p);
    }
  }
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
  

  int nCycles = floor( (p.GetZ() - z) / (slope * 2 * TMath::Pi()));
  z += nCycles * slope * 2 * TMath::Pi();
  
  
  double absC = fabs(slope);
  double currentDistanceZ = fabs(p.GetZ()-z);
  
  if(fabs(p.GetZ()-(z+slope*2*TMath::Pi())) < currentDistanceZ){
    z += absC*2*TMath::Pi();
  }

  return Point(x,y,z);
  
  // This is an attempt to calculate real distance between point and helix:
  //
//  TF1 *fun = new TF1("fun","[0]*sin(x)*([1]-[2]) + [0]*cos(x)*([3]-[4]) + [5]*([6]-[7])-pow([5],2)*x",-TMath::Pi(),TMath::Pi());
//  fun->SetParameter(0, radius);
//  fun->SetParameter(1, origin->GetX());
//  fun->SetParameter(4, origin->GetY());
//  fun->SetParameter(5, slope);
//  fun->SetParameter(7, origin->GetZ());
  //  fun->SetParameter(2, p.GetX());
  //  fun->SetParameter(3, p.GetY());
  //
  //  double pz = p.GetZ();
  ////  cout<<"pz before:"<<pz-origin->GetZ()<<"\t";
  //
  //  int nCycles = floor((p.GetZ() - origin->GetZ())/(slope * 2*TMath::Pi()));
  //  pz -= nCycles * slope * 2*TMath::Pi();
  ////  cout<<"after:"<< (pz-origin->GetZ())/(slope * 2*TMath::Pi())<<endl;
  //
  //  fun->SetParameter(6, pz);
  //
  //  double t1 = fun->GetX(0.0, -TMath::Pi(), 0.0);
  //  double t2 = fun->GetX(0.0, 0.0, TMath::Pi());
  //
  //  double x = radius*cos(t1) + origin->GetX();
  //  double y = radius*sin(t1) + origin->GetY();
  //  double z = slope*t1       + origin->GetZ();
  //
  ////  int nCycles = floor((p.GetZ() - origin->GetZ())/(slope * 2*TMath::Pi()));
  //  z += nCycles * slope * 2*TMath::Pi();
  //
  //  Point q(x,y,z);
  //
  //  if(p.distance(q) > 1.0){
  //    x = radius*cos(t2) + origin->GetX();
  //    y = radius*sin(t2) + origin->GetY();
  //    z = slope*t2       + origin->GetZ();
  //    q = Point(x,y,z);
  //
  //    if(p.distance(q) > 1.0){
  //      cout<<"Point P:"; p.Print();
  //      cout<<"\tPoint Q:"; q.Print();
  //      cout<<"\tdistance:"<<p.distance(q)<<endl;
  //    }
  //  }
}

void Helix::CalculateNregularPoints(int limit)
{
  vector<vector<Point>> pointsByLine = Point::SplitPointsIntoLines(points, config->GetLinesToleranceForRegularity());
  vector<double> possibleDistances;
  set<double> possibleDistancesSet;
  nRegularPoints = 0;
  int nPointsForDistance;
  double zRegularityTolerance = config->GetZregularityTolerance();
  bool first, found;
  double testingDistance;
  
  for(auto line : pointsByLine){
    first=true;
    int iPoint;
    for(iPoint=0; iPoint < line.size()-1; iPoint++){
      testingDistance = line[iPoint].distance(line[iPoint+1]);
      found = false;
      for(double dd : possibleDistances){
        if(fabs(testingDistance-dd) < zRegularityTolerance){found = true;break;}
      }
      if(found) continue;
      
      possibleDistances.push_back(testingDistance);
      nPointsForDistance = 0;
      
      for(auto line2 : pointsByLine){
        for(int i=0;i<line2.size();i++){
          if(std::abs(line2[0].distance(line2[i])-i*testingDistance) < zRegularityTolerance)  nPointsForDistance++;
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
  
  return reasons;
}
