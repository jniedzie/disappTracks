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
slope_valmin(h.slope_valmin),
slope_valmax(h.slope_valmax),
radius_valmin(h.radius_valmin),
radius_valmax(h.radius_valmax),
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
helixParams(h.helixParams)
{
  uniqueID = reinterpret_cast<uint64_t>(this);
  
  for(int i=0;i<kNhelixParams;i++){
    params[i] = h.params[i];
  }
}

Helix Helix::operator=(const Helix &h)
{
  Helix result(h.origin, h.momentum, h.charge);
  
  result.uniqueID = reinterpret_cast<uint64_t>(this);

  result.iCycles        = h.iCycles;
  result.isFinished     = h.isFinished;
  result.slope_valmin   = h.slope_valmin;
  result.slope_valmax   = h.slope_valmax;
  result.radius_valmin  = h.radius_valmin;
  result.radius_valmax  = h.radius_valmax;
  result.seedID         = h.seedID;
  result.points         = h.points;
  result.tShift         = h.tShift;
  result.tMax           = h.tMax;
  result.tStep          = h.tStep;
  result.nRegularPoints = h.nRegularPoints;
  result.nPionPoints    = h.nPionPoints;
  result.vertex         = make_unique<Point>(*h.vertex);
  result.track          = h.track;
  result.radius         = h.radius;
  result.slope          = h.slope;
  result.slopeAbs       = h.slopeAbs;
  result.eventVertex    = h.eventVertex;
  result.helixParams    = h.helixParams;
  
  for(int i=0;i<kNhelixParams;i++){
    result.params[i]  = h.params[i];
  }
  
  return result;
}

Helix::Helix(const Track &_track, const Point &p1, const Point &p2, const Point &_eventVertex) :
origin(Point(0,0,0)),
eventVertex(_eventVertex)
{
  seedID = uniqueID = reinterpret_cast<uint64_t>(this);
  
  track = _track;
  charge = track.GetCharge();
  int nLayers = track.GetNtrackerLayers();
  double Lmin = layerR[nLayers-1];
  double Lmax = layerR[nLayers];
  
  // approximated position of the decay vertex
  vertex = make_unique<Point>((Lmin+Lmax)/2. * cos(track.GetPhi())    + 10*eventVertex.GetX(),
                              (Lmin+Lmax)/2. * sin(track.GetPhi())    + 10*eventVertex.GetY(),
                              (Lmin+Lmax)/2. / tan(track.GetTheta())  + 10*eventVertex.GetZ(),
                              0.0, "",
                              fabs( (Lmax-Lmin)/2. * cos(track.GetPhi()) ),
                              fabs( (Lmax-Lmin)/2. * sin(track.GetPhi()) ),
                              fabs( (Lmax-Lmin)/2. / tan(track.GetTheta()) )
                              );

  points.push_back(make_shared<Point>(*vertex));
  points.push_back(make_shared<Point>(p1));
  points.push_back(make_shared<Point>(p2));
  
  auto circle = circleProcessor.GetCircleFromTriplet(points);
  // this is an approximation, but seems to be ok. It would be perfect to analytically solve the real
  // equations' set for first 3 points on the helix and calculate origin precisely...
  origin = circle->GetCenter();
  
  double t0 = atan2(points[0]->GetY() - origin.GetY(), points[0]->GetX() - origin.GetX());
  double t1 = atan2(points[1]->GetY() - origin.GetY(), points[1]->GetX() - origin.GetX());
  while(t1 < t0) t1 += 2*TMath::Pi();
  double t2 = atan2(points[2]->GetY() - origin.GetY(), points[2]->GetX() - origin.GetX());
  while(t2 < t1) t2 += 2*TMath::Pi();
  
  points[0]->SetT(t0);
  points[1]->SetT(t1);
  points[2]->SetT(t2);
  
  // this only gives approximate direction of the momentum vector
  momentum = make_unique<Point>(p1.GetX() - vertex->GetX(),
                                p1.GetY() - vertex->GetY(),
                                p1.GetZ() - vertex->GetZ());
  
  
  
  vector<Point> p0_variants;
  p0_variants.push_back(Point(Lmin * cos(track.GetPhi())   + 10*eventVertex.GetX(),
                              Lmin * sin(track.GetPhi())   + 10*eventVertex.GetY(),
                              Lmin / tan(track.GetTheta()) + 10*eventVertex.GetZ()));
  
  p0_variants.push_back(Point(Lmax * cos(track.GetPhi())   + 10*eventVertex.GetX(),
                              Lmax * sin(track.GetPhi())   + 10*eventVertex.GetY(),
                              Lmax / tan(track.GetTheta()) + 10*eventVertex.GetZ()));
  
  p0_variants[0].SetT(t0);
  p0_variants[1].SetT(t0);
  
  vector<Point> p1_variants;
  vector<Point> p2_variants;
  
  vector<vector<int>> signs = { { 1, 1, 1},
                                { 1, 1,-1},
                                { 1,-1, 1},
                                { 1,-1,-1},
                                {-1, 1, 1},
                                {-1, 1,-1},
                                {-1,-1, 1},
                                {-1,-1,-1}
  };
  
  for(int i=0;i<8;i++){
    p1_variants.push_back(Point(p1));
    p1_variants[i].SetX(p1.GetX() + signs[i][0] * p1.GetXerr());
    p1_variants[i].SetY(p1.GetY() + signs[i][1] * p1.GetYerr());
    p1_variants[i].SetZ(p1.GetZ() + signs[i][2] * p1.GetZerr());
    p1_variants[i].SetT(t1);
    
    p2_variants.push_back(Point(p2));
    p2_variants[i].SetX(p2.GetX() + signs[i][0] * p2.GetXerr());
    p2_variants[i].SetY(p2.GetY() + signs[i][1] * p2.GetYerr());
    p2_variants[i].SetZ(p2.GetZ() + signs[i][2] * p2.GetZerr());
    p1_variants[i].SetT(t2);
  }
  
  for(Point p0 : p0_variants){
    for(Point p1 : p1_variants){
      for(Point p2 : p2_variants){
        CalcAndUpateHelixParams(p0, p1, p2);
      }
    }
  }

  tStep = 0.01;
  iCycles=0;
  
}

bool Helix::ExtendByPoint(const Point &point)
{
  // make sure that this point is not yet on the helix
  for(auto &p : points){
    if(*p==point) return false;
  }
  
  Point p2 = point;
  double t2 = atan2(p2.GetY() - origin.GetY(), p2.GetX() - origin.GetX());
  while(t2 < GetLastPoint()->GetT()) t2 += 2*TMath::Pi();
  p2.SetT(t2);
  
  vector<Point> limitPoints(kNhelixParams);
  
  for(int iParamNum = 0; iParamNum<kNhelixParams; iParamNum++){
    EHelixParams iParam = static_cast<EHelixParams>(iParamNum);
    Point paramOrigin = GetOrigin(iParam);
    
    // for given helix params, minimizing of maximizing radius or slope, find angle for the new point
    // and calculate limiting point in that plane
    double t = atan2(p2.GetY() - paramOrigin.GetY(), p2.GetX() - paramOrigin.GetX());
    while(t < params[iParam].tMax) t += 2*TMath::Pi();
    
    limitPoints[iParam].SetX(paramOrigin.GetX() + GetRadius(t, iParam) * cos(t));
    limitPoints[iParam].SetY(paramOrigin.GetY() + GetRadius(t, iParam) * sin(t));
    limitPoints[iParam].SetZ(paramOrigin.GetZ() + GetSlope(t, iParam)  * t);
  }
  
  // then, check if new point is between those 4 limiting points (within it's errors)
  double minX=inf, maxX=-inf;
  double minY=inf, maxY=-inf;
  double minZ=inf, maxZ=-inf;
  
  for(Point p : limitPoints){
    if(p.GetX() < minX) minX = p.GetX();
    if(p.GetX() > maxX) maxX = p.GetX();
    if(p.GetY() < minY) minY = p.GetY();
    if(p.GetY() > maxY) maxY = p.GetY();
    if(p.GetZ() < minZ) minZ = p.GetZ();
    if(p.GetZ() > maxZ) maxZ = p.GetZ();
  }
  
  if(   p2.GetX() + p2.GetXerr() + config.helixThickness < minX
     || p2.GetX() - p2.GetXerr() - config.helixThickness > maxX
     || p2.GetY() + p2.GetYerr() + config.helixThickness < minY
     || p2.GetY() - p2.GetYerr() - config.helixThickness > maxY
     || p2.GetZ() + p2.GetZerr() + config.helixThickness < minZ
     || p2.GetZ() - p2.GetZerr() - config.helixThickness > maxZ
     ){
    return false;
  }
  // New idea - if this point is within limiting helices, take this point and all previous ones
  // and try to fit a helix to them minimizing/maximizing slope and radius at inf, giving a huge
  // penalty for missing any of the points by more then its error.
  
  auto fitter = Fitter();
  
  HelixParams newHelixParams[kNhelixParams];
  
  for(int iParam=0; iParam<kNhelixParams; iParam++){
    newHelixParams[iParam] = fitter.FitHelixParams(points, p2, origin, track, eventVertex, (EHelixParams)iParam);
//    newHelixParams[iParam].tShift = params[iParam].tShift;
//    newHelixParams[iParam].tMax   = t2;
//    newHelixParams[iParam].zShift = params[iParam].zShift;
    params[iParam] = newHelixParams[iParam];
  }

  
  /*
  // if the new point it ok, trim its errors to fit within the limiting helices
  TrimPointToLimits(p2, minX, maxX, minY, maxY, minZ, maxZ);
  
  // next, take vertex, the last point on the helix and the new point and narrow down the limits
  Point p0 = *vertex;
  Point p1 = *points[points.size()-1];
  
  double Lmin = layerR[track.GetNtrackerLayers()-1];
  double Lmax = layerR[track.GetNtrackerLayers()];
  
  vector<Point> p0_variants;
  p0_variants.push_back(Point(Lmin * cos(track.GetPhi())   + 10*eventVertex.GetX(),
                              Lmin * sin(track.GetPhi())   + 10*eventVertex.GetY(),
                              Lmin / tan(track.GetTheta()) + 10*eventVertex.GetZ()));
  
  p0_variants.push_back(Point(Lmax * cos(track.GetPhi())   + 10*eventVertex.GetX(),
                              Lmax * sin(track.GetPhi())   + 10*eventVertex.GetY(),
                              Lmax / tan(track.GetTheta()) + 10*eventVertex.GetZ()));
  
  p0_variants[0].SetT(vertex->GetT());
  p0_variants[1].SetT(vertex->GetT());
  
  vector<Point> p1_variants;
  vector<Point> p2_variants;
  
  vector<vector<int>> signs = { { 1, 1, 1},
                                { 1, 1,-1},
                                { 1,-1, 1},
                                { 1,-1,-1},
                                {-1, 1, 1},
                                {-1, 1,-1},
                                {-1,-1, 1},
                                {-1,-1,-1}
  };
  
  for(int i=0;i<8;i++){
    p1_variants.push_back(Point(p1));
    p1_variants[i].SetX(p1.GetX() + signs[i][0] * p1.GetXerr());
    p1_variants[i].SetY(p1.GetY() + signs[i][1] * p1.GetYerr());
    p1_variants[i].SetZ(p1.GetZ() + signs[i][2] * p1.GetZerr());
    
    p2_variants.push_back(Point(p2));
    p2_variants[i].SetX(p2.GetX() + signs[i][0] * p2.GetXerr());
    p2_variants[i].SetY(p2.GetY() + signs[i][1] * p2.GetYerr());
    p2_variants[i].SetZ(p2.GetZ() + signs[i][2] * p2.GetZerr());
  }
  
  vector<HelixParams> newParams;
  
  for(Point p0 : p0_variants){
    for(Point p1 : p1_variants){
      for(Point p2 : p2_variants){
        newParams.push_back(CalcHelixParams(p0, p1, p2));
      }
    }
  }
 
  double new_slopeMin = inf;
  double new_slopeMax = -inf;
  double new_radiusMin = inf;
  double new_radiusMax = -inf;
  
  HelixParams newLimits[kNhelixParams];
  
  double testT = 100*TMath::Pi();
  
  for(HelixParams par : newParams){
    
    
    
    double valR = par.R0 - par.a*testT;
    double valS = par.s0 - par.b*testT;
    
    if(valS < new_slopeMin){new_slopeMin = valS; newLimits[kMinS] = par;}
    if(valS > new_slopeMax){new_slopeMax = valS; newLimits[kMaxS] = par;}
    
    if(valR < new_radiusMin){new_radiusMin = valR; newLimits[kMinR] = par;}
    if(valR > new_radiusMax){new_radiusMax = valR; newLimits[kMaxR] = par;}
  }
  
  // Check if new limits overlap with the current ones
  if(   new_slopeMin > slope_valmax
     || new_slopeMax < slope_valmin
     || new_radiusMin > radius_valmax
     || new_radiusMax < radius_valmin
     ){
    cout<<"new limits don't overlap with the current ones. Failed to add the point"<<endl;
    return false;
  }
  
  // Try to narrow down the limits
  if( new_slopeMin > slope_valmin){   slope_valmin  = new_slopeMin;  params[kMinS] = newLimits[kMinS];}
  if( new_slopeMax < slope_valmax){ // new limit is better at the end of helix
    // check if helix doesn't start to diverge at already existing points
    // and if it does, bring params to values that will make helix go through
    // all the points
    
    double z0 = origin.GetZ();
    double t2 = p2.GetT();
    double z2 = z0 + (newLimits[kMaxS].s0 - newLimits[kMaxS].b*t2)*t2;
    cout<<"Val before:"<<z2<<endl;
    
    for(int iPoint=1; iPoint<points.size(); iPoint++){
      auto testPoint = points[iPoint];
      
      double t = testPoint->GetT();
      double testZ = origin.GetZ() + (newLimits[kMaxS].s0 - newLimits[kMaxS].b*t)*t; // z(t_i) according to updated limits
      
      double zMax = testPoint->GetZ() + testPoint->GetZerr();
      double zMin = testPoint->GetZ() - testPoint->GetZerr();
      
      if(testZ > zMax){
        newLimits[kMaxS].b  = zMax/(t*(t2-t)) - z2/(t2*(t2-t)) - z0/(t*t2);
        newLimits[kMaxS].s0 = (z2-z0)/t2 + newLimits[kMaxS].b*t2;
      }
      if(testZ < zMin){
        newLimits[kMaxS].b  = zMin/(t*(t2-t)) - z2/(t2*(t2-t)) - z0/(t*t2);
        newLimits[kMaxS].s0 = (z2-z0)/t2 + newLimits[kMaxS].b*t2;
      }
      break;
    }
    z2 = z0 + (newLimits[kMaxS].s0 - newLimits[kMaxS].b*t2)*t2;
    cout<<"Val after:"<<z2<<endl;
    
    newLimits[kMaxS].zShift = z2 - (newLimits[kMaxS].s0-newLimits[kMaxS].b*t2)*t2;
    slope_valmax  = newLimits[kMaxS].s0 - newLimits[kMaxS].b*testT;
    params[kMaxS] = newLimits[kMaxS];
  }
  if( new_radiusMin > radius_valmin){ radius_valmin = new_radiusMin; params[kMinR] = newLimits[kMinR];}
  if( new_radiusMax < radius_valmax){ radius_valmax = new_radiusMax; params[kMaxR] = newLimits[kMaxR];}
  
  // If we reached this stage, it means that the new point is a good candidate and we add them to the vector
  // of points on this helix
   
   // Finally, it's time to trim all previous points to the new limits
   // And update limits that were based on the old errors of points
   
  */
   
  points.push_back(make_shared<Point>(p2));

  return true;
}

HelixParams Helix::CalcHelixParams(const Point &p0, const Point &p1, const Point &p2)
{
  double t0 = atan2(p0.GetY() - origin.GetY(), p0.GetX() - origin.GetX());
  double t1 = atan2(p1.GetY() - origin.GetY(), p1.GetX() - origin.GetX());
  while(t1 < t0) t1 += 2*TMath::Pi();
  double t2 = atan2(p2.GetY() - origin.GetY(), p2.GetX() - origin.GetX());
  while(t2 < t1) t2 += 2*TMath::Pi();
  
  double b_num = (p0.GetZ()-p2.GetZ())*(t0-t1) - (p0.GetZ()-p1.GetZ())*(t0-t2);
  double b_den = (t2*t2-t0*t0)*(t0-t1) - (t1*t1-t0*t0)*(t0-t2);
  double b = b_num/b_den;
  double s0 = ( p0.GetZ() - p1.GetZ() - b*(t1*t1-t0*t0)) / (t0-t1);
  
  double a_num = (p0.GetX()-p2.GetX())*(cos(t0)-cos(t1)) - (p0.GetX()-p1.GetX())*(cos(t0)-cos(t2));
  double a_den = (cos(t0)-cos(t1))*(t2*cos(t2)-t0*cos(t0)) - (cos(t0)-cos(t2))*(t1*cos(t1)-t0*cos(t0));
  double a = a_num/a_den;
  double R0 = ( p0.GetX() - p1.GetX() - a*(t1*cos(t1)-t0*cos(t0))) / (cos(t0)-cos(t1));
  
  HelixParams params(R0, a, s0, b);
  
  params.tShift     = t0;
  params.tSecondMax = t1;
  params.tMax       = t2;
  params.zShift = p2.GetZ() - (params.s0-params.b*t2)*t2;
  
  return params;
}

void Helix::CalcAndUpateHelixParams(const Point &p0, const Point &p1, const Point &p2)
{
  HelixParams par = CalcHelixParams(p0, p1, p2);
  
  double testT = 100*TMath::Pi();
  
  double valR = par.R0 - par.a*testT;
  double valS = par.s0 - par.b*testT;
  
  if(valS < slope_valmin){slope_valmin = valS; params[kMinS] = par;}
  if(valS > slope_valmax){slope_valmax = valS; params[kMaxS] = par;}
  
  if(valR < radius_valmin){radius_valmin = valR; params[kMinR] = par;}
  if(valR > radius_valmax){radius_valmax = valR; params[kMaxR] = par;}
}

void Helix::Print()
{
  cout<<"\tseedID: "<<seedID<<endl;
  cout<<"\tuniqueID: "<<uniqueID<<endl;
  cout<<"\tVertex:("<<vertex->GetX()<<","<<vertex->GetY()<<","<<vertex->GetZ()<<")\n";
  cout<<"\tOrigin:("<<origin.GetX()<<","<<origin.GetY()<<","<<origin.GetZ()<<")\n";
  cout<<"\tMomentum:("<<momentum->GetX()<<","<<momentum->GetY()<<","<<momentum->GetZ()<<")\n";
  cout<<"\tCharge: "<<charge<<"\tR:"<<radius<<"\tc:"<<slope<<"\n";
  cout<<"\tnPoints:"<<points.size()<<"\tnPionPoints:"<<nPionPoints<<"\tnRegularPoints:"<<nRegularPoints<<"\n";
  cout<<"\tt min:"<<tShift<<"\tt max:"<<tMax<<endl;
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

double Helix::GetChi2() const
{
  double chi2 = 0;
  for(auto &p : points){
    chi2 += pow(pointsProcessor.distance(*p, GetClosestPoint(*p)), 2);
  }
  return chi2 / points.size();
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

void Helix::TrimPointToLimits(Point &point,
                              double minX, double maxX,
                              double minY, double maxY,
                              double minZ, double maxZ)
{
  double upperX = min(point.GetX() + point.GetXerr(), maxX);
  double lowerX = max(point.GetX() - point.GetXerr(), minX);
  
  double upperY = min(point.GetY() + point.GetYerr(), maxY);
  double lowerY = max(point.GetY() - point.GetYerr(), minY);
  
  double upperZ = min(point.GetZ() + point.GetZerr(), maxZ);
  double lowerZ = max(point.GetZ() - point.GetZerr(), minZ);
  
  point = Point((upperX+lowerX)/2., (upperY+lowerY)/2., (upperZ+lowerZ)/2.,
                point.GetValue(), point.GetSubDetName(),
                (upperX-lowerX)/2., (upperY-lowerY)/2., (upperZ-lowerZ)/2.,
                point.GetT());
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
  
  // approximated position of the decay vertex
 
  points.push_back(make_shared<Point>(_decayVertex));
  for(auto &p : _points) points.push_back(p);
  
  
  // this only gives approximate direction of the momentum vector
  momentum = make_unique<Point>(points[0]->GetX() - vertex->GetX(),
                                points[0]->GetY() - vertex->GetY(),
                                points[0]->GetZ() - vertex->GetZ());
  
  
  tShift  = points[0]->GetT();
  tMax    = points.back()->GetT();
  tStep   = 0.01;
  iCycles =0;
  
}
