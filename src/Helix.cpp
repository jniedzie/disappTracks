//  Helix.cpp
//
//  Created by Jeremi Niedziela on 18/12/2018.

#include "Helix.hpp"
#include "Fitter.hpp"

Helix::Helix(const Point &_decayVertex,
             const Point &_momentum,
             int _charge) :
origin(_decayVertex),
firstTurningPointIndex(-1),
momentum(_momentum),
charge(_charge),
helixParams(HelixParams()),
tStep(0.01),
isFinished(false),
shouldRefit(false),
chi2(inf),
nMissingHits(0),
nMissingHitsInRow(0),
isPreviousHitMissing(false),
nRecLayers(0),
nRecHits(0)
{
  seedID = uniqueID = reinterpret_cast<uint64_t>(this);
  
  helixParams.R0 = GetRadiusInMagField(momentum.GetX(), momentum.GetY(), solenoidField);
  helixParams.s0 = helixParams.R0 * charge * momentum.GetVectorSlopeC();
  helixParams.a = 0;
  helixParams.b = 0;
  
  // take a vector perpendicular to the pion's momentum vector
  Point v = Point(charge * momentum.GetY(),charge * -momentum.GetX(), 0.0);
  const double scale = helixParams.R0/sqrt(pow(v.GetX(),2)+pow(v.GetY(),2));
  
  v.SetX(scale * v.GetX());
  v.SetY(scale * v.GetY());
  
  origin.SetX(origin.GetX() + v.GetX());
  origin.SetY(origin.GetY() + v.GetY());
 
  double tShift=0;
  
  if(momentum.GetZ() > 0){
    if(charge > 0) tShift = TMath::Pi() - atan2(-v.GetX(), v.GetY());
    if(charge < 0) tShift = TMath::Pi() - atan2(-v.GetY(), v.GetX());
  }
  if(momentum.GetZ() < 0){
    if(charge > 0) tShift = TMath::Pi() - atan2(-v.GetY(), v.GetX());
    if(charge < 0) tShift = TMath::Pi() - atan2(-v.GetX(), v.GetY());
  }
  origin.SetZ(origin.GetZ() - fabs(tShift)*fabs(helixParams.s0));
  
  points.push_back(make_shared<Point>(_decayVertex));
  secondToLastPoints.push_back(make_shared<Point>(_decayVertex));
  pointsT.push_back(tShift);
  
  double tMax = GetNcycles()*2*TMath::Pi();
  
  auto lastPoint = make_shared<Point>(origin.GetX() + GetRadius(tMax)*cos(tMax),
                                      origin.GetY() + GetRadius(tMax)*sin(tMax),
                                      -charge*origin.GetZ() + GetSlope(tMax)*tMax);
  
  lastPoint->SetLayer(20);
  
  points.push_back(lastPoint);
  lastPoints.push_back(lastPoint);
  pointsT.push_back(tMax);
}

Helix::Helix(const HelixParams &_params,
             const Point &_decayVertex,
             const Point &_origin,
             int _charge) :
origin(_origin),
firstTurningPointIndex(-1),
momentum(Point(0, 0, 0)),
charge(_charge),
helixParams(_params),
tStep(0.01),
isFinished(false),
shouldRefit(false),
chi2(inf),
nMissingHits(0),
nMissingHitsInRow(0),
isPreviousHitMissing(false),
nRecLayers(0),
nRecHits(0)
{
  seedID = uniqueID = reinterpret_cast<uint64_t>(this);
  
  points.push_back(make_shared<Point>(_decayVertex));
  lastPoints.push_back(make_shared<Point>(_decayVertex));
  pointsT.push_back(pointsProcessor.GetTforPoint(_decayVertex, origin, charge));
}

Helix::Helix(const ROOT::Fit::FitResult &result,
             const Track &track,
             const Point &eventVertex,
             int _charge) :
charge(_charge),
firstTurningPointIndex(-1),
momentum(Point(0, 0, 0)),
tStep(0.01),
isFinished(false),
shouldRefit(false),
chi2(inf),
nMissingHits(0),
nMissingHitsInRow(0),
isPreviousHitMissing(false)
{
  helixParams = HelixParams(result.GetParams()[0],
                            result.GetParams()[1],
                            result.GetParams()[2],
                            result.GetParams()[3]);
  
  double L  = result.GetParams()[4];
  Point _decayVertex = pointsProcessor.GetPointOnTrack(L, track, eventVertex);
  
  seedID = uniqueID = reinterpret_cast<uint64_t>(this);
  
  points.push_back(make_shared<Point>(_decayVertex));
  lastPoints.push_back(make_shared<Point>(_decayVertex));
  pointsT.push_back(pointsProcessor.GetTforPoint(_decayVertex, origin, charge));
  
  double x0 = result.GetParams()[5];
  double y0 = result.GetParams()[6];
  double z0 = result.GetParams()[7];
  
  origin = Point(x0, y0, z0);
  UpdateOrigin(origin);
  seedChi2 = chi2 = result.MinFcnValue();
}


Helix::Helix(const Helix &h) :
isFinished(h.isFinished),
seedID(h.seedID),
points(h.points),
lastPoints(h.lastPoints),
secondToLastPoints(h.secondToLastPoints),
thirdToLastPoints(h.thirdToLastPoints),
pointsT(h.pointsT),
tStep(h.tStep),
origin(h.origin),
momentum(h.momentum),
charge(h.charge),
helixParams(h.helixParams),
chi2(h.chi2),
shouldRefit(h.shouldRefit),
nMissingHits(h.nMissingHits),
nMissingHitsInRow(h.nMissingHitsInRow),
isPreviousHitMissing(h.isPreviousHitMissing),
firstTurningPointIndex(h.firstTurningPointIndex),
nRecLayers(h.nRecLayers),
nRecHits(h.nRecHits)
{
  uniqueID = reinterpret_cast<uint64_t>(this);
}

Helix& Helix::operator=(const Helix &h)
{
  origin   = h.origin;
  momentum = h.momentum;
  charge   = h.charge;
  
  uniqueID = reinterpret_cast<uint64_t>(this);

  isFinished     = h.isFinished;
  seedID         = h.seedID;
  points         = h.points;
  lastPoints     = h.lastPoints;
  secondToLastPoints = h.secondToLastPoints;
  thirdToLastPoints = h.thirdToLastPoints;
  pointsT        = h.pointsT;
  tStep          = h.tStep;
  helixParams    = h.helixParams;
  chi2           = h.chi2;
  shouldRefit    = h.shouldRefit;
  firstTurningPointIndex = h.firstTurningPointIndex;
 
  nMissingHits         = h.nMissingHits;
  nMissingHitsInRow    = h.nMissingHitsInRow;
  isPreviousHitMissing = h.isPreviousHitMissing;
  
  nRecLayers  = h.nRecLayers;
  nRecHits    = h.nRecHits;
  
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
  cout<<"\tMomentum:("<<momentum.GetX()<<","<<momentum.GetY()<<","<<momentum.GetZ()<<")\n";
  cout<<"\tCharge: "<<charge<<"\n";
  cout<<"\ts0: "<<helixParams.s0<<"\tb: "<<helixParams.b<<"\tR0: "<<helixParams.R0<<"\t a: "<<helixParams.a<<endl;
  cout<<"\tnPoints:"<<points.size()<<"\tnLayers: "<<GetNlayers()<<"\tnMissingHits: "<<nMissingHits<<"\n";
  cout<<"\tt min:"<<GetTmin()<<"\tt max:"<<GetTmax()<<endl;
  cout<<"\tchi2:"<<chi2<<endl;
}

double Helix::GetNcycles() const
{
  return sgn(momentum.GetZ())*((sgn(momentum.GetZ())*trackerZsize) - origin.GetZ())/(fabs(GetSlope(GetTmin()))*2*TMath::Pi());
}

double Helix::GetRadius(double t) const
{
  return GetRofT(helixParams.R0, helixParams.a, GetTmin(), t, charge);
}

double Helix::GetSlope(double t) const
{
  return GetSofT(helixParams.s0, helixParams.b, GetTmin(), t, charge);
}

void Helix::AddPoint(const shared_ptr<Point> &point, vector<size_t> *_lastPointIndices)
{
  double t = pointsProcessor.GetTforPoint(*point, origin, charge);
  double tMax = GetTmax(_lastPointIndices);
  
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

void Helix::RemoveLastPoints()
{
  if(points.back()->GetSubDetName()=="missing"){
    nMissingHits--;
    nMissingHitsInRow--;
  }
  for(auto &p : lastPoints){
    for(size_t i=0; i<points.size();){
      if(*p == *points[i]){
        points.erase(points.begin()+i);
        pointsT.erase(pointsT.begin()+i);
      }
      else i++;
    }
  }
  
  lastPoints = secondToLastPoints;
  secondToLastPoints = thirdToLastPoints;
}

vector<size_t> Helix::GetLastPointsIndices() const
{
  vector<size_t> resultIndices;

  for(int iPoint=0; iPoint<points.size(); iPoint++){
    for(auto &p : lastPoints){
      if(*p == *points[iPoint]) resultIndices.push_back(iPoint);
    }
  }
  return resultIndices;
}

void Helix::SetLastPoints(Points _points)
{
  vector<size_t> lastPointsIndices = GetLastPointsIndices();
  thirdToLastPoints = secondToLastPoints;
  secondToLastPoints = lastPoints;
  lastPoints.clear();
  
  for(auto &p : _points){
    AddPoint(p, &lastPointsIndices);
    lastPoints.push_back(p);
  }
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
    int disk = points[iPoint]->GetDisk();
    if(disk >= 0){
      layers.insert(disk);
    }
  }
  return layers.size();
}

double Helix::GetTmax(vector<size_t> *_lastPointIndices) const
{
  vector<size_t> lastPointsIndices;
  
  if(_lastPointIndices) lastPointsIndices = *_lastPointIndices;
  else lastPointsIndices = GetLastPointsIndices();
  
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
  
  Points sortedPoints;
  vector<double> sortedPointsT;
  
  for(auto iter : idx){
    sortedPoints.push_back(points[iter]);
    sortedPointsT.push_back(pointsT[iter]);
  }
  points = sortedPoints;
  pointsT = sortedPointsT;
}

void Helix::SetPointsAndSortByT(const Points &_points, const vector<double> &_pointsT)
{
  points = _points;
  pointsT = _pointsT;
  SortPointsByT(charge > 0);
  
  lastPoints.clear();
  
  auto lastPoint = points.back();
  if(lastPoint->IsEndcapHit()){
    int lastPointDisk = lastPoint->GetDisk();
    for(auto &p : points){
      if(p->GetDisk() == lastPointDisk) lastPoints.push_back(p);
    }
  }
  else{
    int lastPointLayer = lastPoint->GetLayer();
    for(auto &p : points){
      if(p->GetLayer() == lastPointLayer) lastPoints.push_back(p);
    }
  }
}
