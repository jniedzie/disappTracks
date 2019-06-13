//
//  Helix.cpp
//
//  Created by Jeremi Niedziela on 18/12/2018.
//

#include "Helix.hpp"
#include "Fitter.hpp"

Helix::Helix(const Point &_origin,
             const unique_ptr<Point> &_momentum,
             int _charge) :
vertex(make_unique<Point>(_origin)),
origin(_origin),
momentum(make_unique<Point>(*_momentum)),
charge(_charge),
iCycles(0),
isFinished(false),
helixParams(HelixParams()),
chi2(inf),
shouldRefit(false),
nMissingHits(0),
nMissingHitsInRow(0),
isPreviousHitMissing(false),
firstTurningPointIndex(-1)
{
  seedID = uniqueID = reinterpret_cast<uint64_t>(this);
  
  helixParams.R0 = GetRadiusInMagField(momentum->GetX(), momentum->GetY(), solenoidField);
  helixParams.s0 = GetRadius(0) * charge * momentum->GetVectorSlopeC();
  
  // take a vector perpendicular to the pion's momentum vector
  Point v = Point(charge * momentum->GetY(),charge * -momentum->GetX(), 0.0);
  const double scale = GetRadius(0)/sqrt(pow(v.GetX(),2)+pow(v.GetY(),2));
  
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
  origin.SetZ(origin.GetZ() - fabs(tShift)*fabs(GetSlope(0)));
  
  tMax = GetNcycles()*2*TMath::Pi();
  tStep = 0.01;
}

Helix::Helix(const HelixParams &_params,
             const Point &_decayVertex,
             const Point &_origin,
             const vector<shared_ptr<Point>> &_points,
             const Track &_track) :
origin(_origin),
vertex(make_unique<Point>(_decayVertex)),
helixParams(_params),
track(_track),
iCycles(0),
isFinished(false),
chi2(inf),
shouldRefit(false),
nMissingHits(0),
nMissingHitsInRow(0),
isPreviousHitMissing(false),
firstTurningPointIndex(-1)
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
}

Helix::Helix(const Helix &h) :
iCycles(h.iCycles),
isFinished(h.isFinished),
seedID(h.seedID),
points(h.points),
tShift(h.tShift),
tMax(h.tMax),
tStep(h.tStep),
vertex(make_unique<Point>(*h.vertex)),
origin(h.origin),
momentum(make_unique<Point>(*h.momentum)),
track(h.track),
charge(h.charge),
helixParams(h.helixParams),
chi2(h.chi2),
shouldRefit(h.shouldRefit),
nMissingHits(h.nMissingHits),
nMissingHitsInRow(h.nMissingHitsInRow),
isPreviousHitMissing(h.isPreviousHitMissing),
firstTurningPointIndex(h.firstTurningPointIndex)
{
  uniqueID = reinterpret_cast<uint64_t>(this);
}

Helix& Helix::operator=(const Helix &h)
{
  origin = h.origin;
  momentum = make_unique<Point>(*h.momentum);
  charge = h.charge;
  
  uniqueID = reinterpret_cast<uint64_t>(this);

  iCycles        = h.iCycles;
  isFinished     = h.isFinished;
  seedID         = h.seedID;
  points         = h.points;
  tShift         = h.tShift;
  tMax           = h.tMax;
  tStep          = h.tStep;
  vertex         = make_unique<Point>(*h.vertex);
  track          = h.track;
  helixParams    = h.helixParams;
  chi2           = h.chi2;
  shouldRefit    = h.shouldRefit;
  firstTurningPointIndex = h.firstTurningPointIndex;
 
  nMissingHits         = h.nMissingHits;
  nMissingHitsInRow    = h.nMissingHitsInRow;
  isPreviousHitMissing = h.isPreviousHitMissing;
  
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
  cout<<"\tCharge: "<<charge<<"\n";
  cout<<"\ts0: "<<helixParams.s0<<"\tb: "<<helixParams.b<<"\tR0: "<<helixParams.R0<<"\t a: "<<helixParams.a<<endl;
  cout<<"\tnPoints:"<<points.size()<<"\tnMissingHits: "<<nMissingHits<<"\n";
  cout<<"\tt min:"<<tShift<<"\tt max:"<<tMax<<endl;
  cout<<"\tchi2:"<<chi2<<endl;
}

void Helix::SetVertex(const Point &_vertex)
{
  points[0] = make_shared<Point>(_vertex);
  vertex    = make_unique<Point>(_vertex);
}

double Helix::GetNcycles() const
{
  return sgn(momentum->GetZ())*((sgn(momentum->GetZ())*trackerZsize) - origin.GetZ())/(fabs(GetSlope(0))*2*TMath::Pi());
}

double Helix::GetRadius(double t) const
{
  return (helixParams.R0 - helixParams.a*t);
}

double Helix::GetSlope(double t) const
{
  return (helixParams.s0 - helixParams.b*t);
}

void Helix::AddPoint(const shared_ptr<Point> &point)
{
  shared_ptr<Point> previousPoint = points.back();
  
  int previousPointLayer = previousPoint->GetLayer();
  
  double t = pointsProcessor.GetTforPoint(*point, origin, charge);
  int thisPointLayer = point->GetLayer();
  
  if((thisPointLayer != previousPointLayer) || (points.size()==firstTurningPointIndex)){
    if(charge < 0) while(t < tMax) t += 2*TMath::Pi();
    else           while(t > tMax) t -= 2*TMath::Pi();
  }
  else{
    double previousT = previousPoint->GetT();
         if(fabs(t-previousT) > fabs(t+2*TMath::Pi()-previousT)) t += 2*TMath::Pi();
    else if(fabs(t-previousT) > fabs(t-2*TMath::Pi()-previousT)) t -= 2*TMath::Pi();
  }
  point->SetT(t);
  tMax = t;
  
  points.push_back(point);
}

void Helix::UpdateOrigin(const Point &_origin)
{
  origin = _origin;
  
  pointsProcessor.SetPointsT(points, origin, charge);
  
  tShift = points.front()->GetT();
  tMax   = points.back()->GetT();
}

void Helix::IncreaseMissingHits(){
  nMissingHits++;
  if(isPreviousHitMissing) nMissingHitsInRow++;
  else                     nMissingHitsInRow=1;
  
  isPreviousHitMissing = true;
}

vector<shared_ptr<Point>> Helix::GetLastPoints() const
{
  vector<shared_ptr<Point>> resultPoints;
  int lastPointLayer = points.back()->GetLayer();
  
  for(int iPoint=0; iPoint<points.size(); iPoint++){
    auto point = points[iPoint];
    
    if((point->GetLayer() == lastPointLayer) &&
       (iPoint>=firstTurningPointIndex)){
      resultPoints.push_back(point);
    }
  }
  
  return resultPoints;
}

vector<shared_ptr<Point>> Helix::GetSecontToLastPoints() const
{
  vector<shared_ptr<Point>> resultPoints;
  
  int lastPointLayer = points.back()->GetLayer();
  
  for(int iPoint=0; iPoint<points.size(); iPoint++){
    auto point = points[iPoint];
    
    if(firstTurningPointIndex<0){
      if((point->GetLayer() == lastPointLayer-1)) resultPoints.push_back(point);
    }
    else if(iPoint<firstTurningPointIndex){
      if((point->GetLayer() == lastPointLayer)) resultPoints.push_back(point);
    }
    else{
      if((point->GetLayer() == lastPointLayer+1)) resultPoints.push_back(point);
    }
  }
  
  return resultPoints;
}
