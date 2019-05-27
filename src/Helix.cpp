//
//  Helix.cpp
//
//  Created by Jeremi Niedziela on 18/12/2018.
//

#include "Helix.hpp"
#include "Fitter.hpp"

Helix::Helix() :
origin(0,0,0),
eventVertex(0,0,0)
{
  
}

Helix::Helix(const Point &_origin,
             const unique_ptr<Point> &_momentum,
             int _charge) :
vertex(make_unique<Point>(_origin)),
origin(_origin),
momentum(make_unique<Point>(*_momentum)),
charge(_charge),
eventVertex(0,0,0)
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

Helix::Helix(const Helix &h) :
iCycles(h.iCycles),
isFinished(h.isFinished),
seedID(h.seedID),
points(h.points),
tShift(h.tShift),
tMax(h.tMax),
tStep(h.tStep),
nRegularPoints(h.nRegularPoints),
nPionPoints(h.nPionPoints),
vertex(make_unique<Point>(*h.vertex)),
origin(h.origin),
momentum(make_unique<Point>(*h.momentum)),
track(h.track),
radius(h.radius),
slope(h.slope),
slopeAbs(h.slopeAbs),
charge(h.charge),
eventVertex(h.eventVertex),
helixParams(h.helixParams),
chi2(h.chi2),
increasing(h.increasing),
shouldRefit(h.shouldRefit)
{
  uniqueID = h.uniqueID;
//  uniqueID = reinterpret_cast<uint64_t>(this);
}

Helix& Helix::operator=(const Helix &h)
{
  origin = h.origin;
  momentum = make_unique<Point>(*h.momentum);
  charge = h.charge;
  
  uniqueID = h.uniqueID;
//  uniqueID = reinterpret_cast<uint64_t>(this);

  iCycles        = h.iCycles;
  isFinished     = h.isFinished;
  seedID         = h.seedID;
  points         = h.points;
  tShift         = h.tShift;
  tMax           = h.tMax;
  tStep          = h.tStep;
  nRegularPoints = h.nRegularPoints;
  nPionPoints    = h.nPionPoints;
  vertex         = make_unique<Point>(*h.vertex);
  track          = h.track;
  radius         = h.radius;
  slope          = h.slope;
  slopeAbs       = h.slopeAbs;
  eventVertex    = h.eventVertex;
  
  helixParams    = h.helixParams;
  chi2           = h.chi2;
  increasing     = h.increasing;
  shouldRefit    = h.shouldRefit;
 
  return *this;
}

bool Helix::operator==(const Helix &h)
{
  if(points.size() == h.points.size() &&
     fabs(helixParams.R0 - h.helixParams.R0) < 0.00001 &&
     fabs(helixParams.s0 - h.helixParams.s0) < 0.00001 &&
     fabs(helixParams.a - h.helixParams.a) < 0.00001 &&
     fabs(helixParams.b - h.helixParams.b) < 0.00001 &&
     origin == h.origin &&
     vertex == h.vertex) return true;
  
  return false;
}

void Helix::Print()
{
  cout<<"\tseedID: "<<seedID<<endl;
  cout<<"\tuniqueID: "<<uniqueID<<endl;
  cout<<"\tVertex:("<<vertex->GetX()<<","<<vertex->GetY()<<","<<vertex->GetZ()<<")\n";
  cout<<"\tOrigin:("<<origin.GetX()<<","<<origin.GetY()<<","<<origin.GetZ()<<")\n";
  cout<<"\tMomentum:("<<momentum->GetX()<<","<<momentum->GetY()<<","<<momentum->GetZ()<<")\n";
  cout<<"\tCharge: "<<charge<<"\tR:"<<radius<<"\tc:"<<slope<<"\n";
  cout<<"\ts0: "<<helixParams.s0<<"\tb: "<<helixParams.b<<"\tR0: "<<helixParams.R0<<"\t a: "<<helixParams.a<<endl;
  cout<<"\tnPoints:"<<points.size()<<"\tnPionPoints:"<<nPionPoints<<"\tnRegularPoints:"<<nRegularPoints<<"\n";
  cout<<"\tt min:"<<tShift<<"\tt max:"<<tMax<<endl;
  cout<<"\tchi2:"<<chi2<<endl;
}

void Helix::SetPoints(const vector<shared_ptr<Point>> &_points)
{
  nPionPoints = 0;
  points.clear();
  
  for(auto &p : _points){
    Point q = GetClosestPoint(*p);
    if(pointsProcessor.distance(*p,q) < config.helixThickness){
      if(p->IsPionHit()) nPionPoints++;
      points.push_back(p);
    }
  }
}

Point Helix::GetClosestPoint(const Point &p) const
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


//------------------------------------------------------------------------------------------------
// Annother approach
//

Helix::Helix(const HelixParams &_params,
             const Point &_decayVertex,
             const Point &_origin,
             const vector<shared_ptr<Point>> &_points,
             const Track &_track) :
origin(_origin),
vertex(make_unique<Point>(_decayVertex)),
helixParams(_params),
track(_track)
{
  seedID = uniqueID = reinterpret_cast<uint64_t>(this);
  charge = track.GetCharge();
 
  points.push_back(make_shared<Point>(_decayVertex));
  for(auto &p : _points) points.push_back(p);
  
  // this only gives approximate direction of the momentum vector
  momentum = make_unique<Point>(points[1]->GetX() - vertex->GetX(),
                                points[1]->GetY() - vertex->GetY(),
                                points[1]->GetZ() - vertex->GetZ());
  
  
  tShift  = points.front()->GetT();
  tMax    = points.back()->GetT();
  tStep   = 0.01;
  iCycles =0;
}

void Helix::AddPoint(const shared_ptr<Point> &point)
{
  points.push_back(point);
  
  double t;
  if(charge < 0){
    t = atan2(point->GetY() - origin.GetY(), point->GetX() - origin.GetX());
    // FIX this:
    while(t < tMax) t += 2*TMath::Pi();
  }
  else{
    t = TMath::Pi()/2 - atan2(point->GetX() - origin.GetX(), point->GetY() - origin.GetY());
    // FIX this:
    while(t > tMax) t -= 2*TMath::Pi();
  }
  
  point->SetT(t);
  tMax = t;
}

void Helix::UpdateOrigin(const Point &_origin)
{
  origin = _origin;
  
  for(auto &p : points){
    double t;
    if(charge < 0)  t = atan2(p->GetY() - origin.GetY(), p->GetX() - origin.GetX());
    else            t = TMath::Pi()/2 - atan2(p->GetX() - origin.GetX(), p->GetY() - origin.GetY());
    p->SetT(t);
  }
  
  sort(points.begin(), points.end(), PointsProcessor::ComparePointByLayer());
//  sort(points.begin(), points.end(), PointsProcessor::ComparePointByT());

  // this should be fixed:
  for(int iPoint=0; iPoint<points.size()-1; iPoint++){
    if(charge < 0){
      while(points[iPoint+1]->GetT() < points[iPoint]->GetT()){
        points[iPoint+1]->SetT(points[iPoint+1]->GetT() + 2*TMath::Pi());
      }
    }
    else{
      while(points[iPoint+1]->GetT() > points[iPoint]->GetT()){
        points[iPoint+1]->SetT(points[iPoint+1]->GetT() - 2*TMath::Pi());
      }
    }
  }
  
  tShift = points.front()->GetT();
  tMax   = points.back()->GetT();
}
