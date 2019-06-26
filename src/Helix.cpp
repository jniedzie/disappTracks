//
//  Helix.cpp
//
//  Created by Jeremi Niedziela on 18/12/2018.
//

#include "Helix.hpp"
#include "Fitter.hpp"

Helix::Helix(const Point &_decayVertex,
             const unique_ptr<Point> &_momentum,
             int _charge) :
origin(_decayVertex),
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
  points.push_back(make_shared<Point>(_decayVertex));
  
  seedID = uniqueID = reinterpret_cast<uint64_t>(this);
  
  helixParams.R0 = GetRadiusInMagField(momentum->GetX(), momentum->GetY(), solenoidField);
  helixParams.s0 = GetRadius(0) * charge * momentum->GetVectorSlopeC();
  helixParams.a = 0;
  helixParams.b = 0;
  
  // take a vector perpendicular to the pion's momentum vector
  Point v = Point(charge * momentum->GetY(),charge * -momentum->GetX(), 0.0);
  const double scale = GetRadius(0)/sqrt(pow(v.GetX(),2)+pow(v.GetY(),2));
  
  v.SetX(scale * v.GetX());
  v.SetY(scale * v.GetY());
  
  origin.SetX(origin.GetX() + v.GetX());
  origin.SetY(origin.GetY() + v.GetY());
 
  double tShift=0;
  
  if(momentum->GetZ() > 0){
    if(charge > 0) tShift = TMath::Pi() - atan2(-v.GetX(), v.GetY());
    if(charge < 0) tShift = TMath::Pi() - atan2(-v.GetY(), v.GetX());
  }
  if(momentum->GetZ() < 0){
    if(charge > 0) tShift = TMath::Pi() - atan2(-v.GetY(), v.GetX());
    if(charge < 0) tShift = TMath::Pi() - atan2(-v.GetX(), v.GetY());
  }
  origin.SetZ(origin.GetZ() - fabs(tShift)*fabs(GetSlope(0)));
  
  points.front()->SetT(tShift);
  
  double tMax = GetNcycles()*2*TMath::Pi();
  
  auto lastPoint = make_shared<Point>(origin.GetX() + GetRadius(tMax)*cos(tMax),
                                      origin.GetY() + GetRadius(tMax)*sin(tMax),
                                      -charge*origin.GetZ() + GetSlope(tMax)*tMax);
  
  lastPoint->SetT(tMax);
  lastPoint->SetLayer(20);
  points.push_back(lastPoint);
  
  tStep = 0.01;
}

Helix::Helix(const HelixParams &_params,
             const Point &_decayVertex,
             const Point &_origin,
             const vector<shared_ptr<Point>> &_points,
             const Track &_track) :
origin(_origin),
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
  for(auto &p : _points) points.push_back(make_shared<Point>(*p));
  
  // this only gives approximate direction of the momentum vector
  momentum = make_unique<Point>(points[1]->GetX() - points[0]->GetX(),
                                points[1]->GetY() - points[0]->GetY(),
                                points[1]->GetZ() - points[0]->GetZ());
  
  tStep   = 0.01;
}

Helix::Helix(const Helix &h) :
iCycles(h.iCycles),
isFinished(h.isFinished),
seedID(h.seedID),
points(h.points),
tStep(h.tStep),
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
  tStep          = h.tStep;
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
     origin == h.origin) return true;
  
  return false;
}

void Helix::Print()
{
  cout<<"\tseedID: "<<seedID<<endl;
  cout<<"\tuniqueID: "<<uniqueID<<endl;
  if(points.size()>0){ cout<<"\tVertex: "; points.front()->Print(); cout<<"\n";}
  cout<<"\tOrigin:("<<origin.GetX()<<","<<origin.GetY()<<","<<origin.GetZ()<<")\n";
  cout<<"\tMomentum:("<<momentum->GetX()<<","<<momentum->GetY()<<","<<momentum->GetZ()<<")\n";
  cout<<"\tCharge: "<<charge<<"\n";
  cout<<"\ts0: "<<helixParams.s0<<"\tb: "<<helixParams.b<<"\tR0: "<<helixParams.R0<<"\t a: "<<helixParams.a<<endl;
  cout<<"\tnPoints:"<<points.size()<<"\tnLayers: "<<GetNlayers()<<"\tnMissingHits: "<<nMissingHits<<"\n";
  cout<<"\tt min:"<<GetTmin()<<"\tt max:"<<GetTmax()<<endl;
  cout<<"\tchi2:"<<chi2<<endl;
}

void Helix::SetVertex(const Point &_vertex)
{
  points[0] = make_shared<Point>(_vertex);
}

double Helix::GetNcycles() const
{
  return sgn(momentum->GetZ())*((sgn(momentum->GetZ())*trackerZsize) - origin.GetZ())/(fabs(GetSlope(0))*2*TMath::Pi());
}

double Helix::GetRadius(double t) const
{
  if(config.expRadiusFunction) return helixParams.R0 * exp(-helixParams.a*t);
  return (helixParams.R0 - helixParams.a*t);
}

double Helix::GetSlope(double t) const
{
  if(config.expRadiusFunction) return helixParams.s0 * exp(-helixParams.b*t);
  return (helixParams.s0 - helixParams.b*t);
}

void Helix::AddPoint(const shared_ptr<Point> &point)
{
  double t = pointsProcessor.GetTforPoint(*point, origin, charge);
  double tMax = GetTmax();
  
  while(fabs(t-tMax) > fabs(t+2*TMath::Pi()-tMax)) t += 2*TMath::Pi();
  while(fabs(t-tMax) > fabs(t-2*TMath::Pi()-tMax)) t -= 2*TMath::Pi();
  
  point->SetT(t);
  points.push_back(make_shared<Point>(*point));
}

void Helix::UpdateOrigin(const Point &_origin)
{
  origin = _origin;
  pointsProcessor.SetPointsT(points, origin, charge);
}

void Helix::IncreaseMissingHits(){
  nMissingHits++;
  if(isPreviousHitMissing) nMissingHitsInRow++;
  else                     nMissingHitsInRow=1;
  
  isPreviousHitMissing = true;
}

void Helix::RemoveLastPoint()
{
  points.pop_back();
}

vector<shared_ptr<Point>> Helix::GetLastPoints() const
{
  vector<shared_ptr<Point>> resultPoints;
  int lastPointLayer = points.back()->GetLayer();
  
  for(int iPoint=0; iPoint<points.size(); iPoint++){
    auto point = points[iPoint];
    
    if((point->GetLayer() == lastPointLayer) &&
       (iPoint>=firstTurningPointIndex)){
      resultPoints.push_back(make_shared<Point>(*point));
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
      if((point->GetLayer() == lastPointLayer-1)) resultPoints.push_back(make_shared<Point>(*point));
    }
    else if(iPoint<firstTurningPointIndex){
      if((point->GetLayer() == lastPointLayer)) resultPoints.push_back(make_shared<Point>(*point));
    }
    else{
      if((point->GetLayer() == lastPointLayer+1)) resultPoints.push_back(make_shared<Point>(*point));
    }
  }
  
  return resultPoints;
}

uint Helix::GetNlayers() const
{
  unordered_set<int> layers;
  for(int iPoint=0; iPoint<points.size(); iPoint++){
    int layer = points[iPoint]->GetLayer();
    if(layer > 0){
      if(iPoint < firstTurningPointIndex) layers.insert( layer);
      else                                layers.insert(-layer);
    }
  }
  return layers.size();
}

double Helix::GetTmax() const
{
  vector<shared_ptr<Point>> lastPoints = GetLastPoints();
  if(lastPoints.size()==1) return lastPoints.front()->GetT();
  
  double extremeT=charge*inf;
  
  for(auto &point : lastPoints){
    if(charge > 0 && point->GetT() < extremeT) extremeT = point->GetT();
    if(charge < 0 && point->GetT() > extremeT) extremeT = point->GetT();
  }

  return extremeT;
}
  
