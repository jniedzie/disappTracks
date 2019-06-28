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
  seedID = uniqueID = reinterpret_cast<uint64_t>(this);
  
  helixParams.R0 = GetRadiusInMagField(momentum->GetX(), momentum->GetY(), solenoidField);
  helixParams.s0 = helixParams.R0 * charge * momentum->GetVectorSlopeC();
  helixParams.a = 0;
  helixParams.b = 0;
  
  // take a vector perpendicular to the pion's momentum vector
  Point v = Point(charge * momentum->GetY(),charge * -momentum->GetX(), 0.0);
  const double scale = helixParams.R0/sqrt(pow(v.GetX(),2)+pow(v.GetY(),2));
  
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
  origin.SetZ(origin.GetZ() - fabs(tShift)*fabs(helixParams.s0));
  
  points.push_back(make_shared<Point>(_decayVertex));
  pointsT.push_back(tShift);
  
  double tMax = GetNcycles()*2*TMath::Pi();
  
  auto lastPoint = make_shared<Point>(origin.GetX() + GetRadius(tMax)*cos(tMax),
                                      origin.GetY() + GetRadius(tMax)*sin(tMax),
                                      -charge*origin.GetZ() + GetSlope(tMax)*tMax);
  
  lastPoint->SetLayer(20);
  
  points.push_back(lastPoint);
  pointsT.push_back(tMax);
  
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
  pointsT.push_back(pointsProcessor.GetTforPoint(_decayVertex, origin, charge));
  
  for(auto &p : _points){
    points.push_back(p);
    pointsT.push_back(pointsProcessor.GetTforPoint(*p, origin, charge));
  }
  
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
pointsT(h.pointsT),
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
  pointsT        = h.pointsT;
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
  return sgn(momentum->GetZ())*((sgn(momentum->GetZ())*trackerZsize) - origin.GetZ())/(fabs(GetSlope(GetTmin()))*2*TMath::Pi());
}

double Helix::GetRadius(double t) const
{
  return GetRofT(helixParams.R0, helixParams.a, GetTmin(), t, charge);
}

double Helix::GetSlope(double t) const
{
  return GetSofT(helixParams.s0, helixParams.b, GetTmin(), t, charge);
}

void Helix::AddPoint(const shared_ptr<Point> &point)
{
  double t = pointsProcessor.GetTforPoint(*point, origin, charge);
  double tMax = GetTmax();
  
  if(charge > 0){ if(t > tMax) t -= 2*TMath::Pi(); }
  else          { if(t < tMax) t += 2*TMath::Pi(); }

  points.push_back(point);
  pointsT.push_back(t);
}

void Helix::UpdateOrigin(const Point &_origin)
{
  origin = _origin;
  
  for(int iPoint=0; iPoint<points.size(); iPoint++){
    auto point = points[iPoint];
    double t = pointsProcessor.GetTforPoint(*point, origin, charge);
    if(iPoint!=0){
      double previousT = pointsT[iPoint-1];
      while(fabs(t+2*TMath::Pi()-previousT) < fabs(t-previousT)) t += 2*TMath::Pi();
      while(fabs(t-2*TMath::Pi()-previousT) < fabs(t-previousT)) t -= 2*TMath::Pi();
    }
    pointsT[iPoint] = t;
  }
}

void Helix::IncreaseMissingHits(){
  nMissingHits++;
  if(isPreviousHitMissing) nMissingHitsInRow++;
  else                     nMissingHitsInRow=1;
  
  isPreviousHitMissing = true;
}

void Helix::RemoveLastPoint()
{
  if(points.back()->GetSubDetName()=="missing"){
    nMissingHits--;
    nMissingHitsInRow--;
  }
  points.pop_back();
  pointsT.pop_back();
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

vector<size_t> Helix::GetLastPointsIndices() const
{
  vector<size_t> resultIndices;
  int lastPointLayer = points.back()->GetLayer();
  
  for(int iPoint=0; iPoint<points.size(); iPoint++){
    if((points[iPoint]->GetLayer() == lastPointLayer) &&
       (iPoint>=firstTurningPointIndex)){
      resultIndices.push_back(iPoint);
    }
  }
  return resultIndices;
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

size_t Helix::GetNlayers() const
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
  vector<size_t> lastPointsIndices = GetLastPointsIndices();
  if(lastPointsIndices.size()==1) return pointsT[lastPointsIndices.front()];
  
  double extremeT=charge*inf;
  
  for(size_t index : lastPointsIndices){
    if(charge > 0 && pointsT[index] < extremeT) extremeT = pointsT[index];
    if(charge < 0 && pointsT[index] > extremeT) extremeT = pointsT[index];
  }

  return extremeT;
}

void Helix::SortPointsByT(bool inverted)
{
  // initialize original index locations
  vector<size_t> idx(points.size());
  iota(idx.begin(), idx.end(), 0);
  
  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&](size_t i1, size_t i2) {return inverted ? pointsT[i1] > pointsT[i2] : pointsT[i1] < pointsT[i2];});
  
  vector<shared_ptr<Point>> sortedPoints;
  vector<double> sortedPointsT;
  
  for(auto iter : idx){
    sortedPoints.push_back(points[iter]);
    sortedPointsT.push_back(pointsT[iter]);
  }
  points = sortedPoints;
  pointsT = sortedPointsT;
}

void Helix::SetPointsAndSortByT(const vector<shared_ptr<Point>> &_points, const vector<double> &_pointsT)
{
  points = _points;
  pointsT = _pointsT;
  SortPointsByT(charge > 0);
}
