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

Helix::Helix(const Helix &h) :
Lmin(h.Lmin),
Lmax(h.Lmax),
bmin(h.bmin),
bmax(h.bmax),
s0min(h.s0min),
s0max(h.s0max),
amin(h.amin),
amax(h.amax),
R0min(h.R0min),
R0max(h.R0max),
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
pointsProcessor(make_unique<PointsProcessor>()),
circleProcessor(make_unique<CircleProcessor>())
{
  uniqueID = reinterpret_cast<uint64_t>(this);
}

Helix::Helix(const Track &_track, const Point &p1, const Point &p2, const Point &_eventVertex) :
origin(Point(0,0,0)),
circleProcessor(make_unique<CircleProcessor>())
{
  seedID = uniqueID = reinterpret_cast<uint64_t>(this);
  
  track = _track;
  
  int nLayers = track.GetNtrackerLayers();
  
  Lmin = layerR[nLayers-1];
  Lmax = layerR[nLayers];
  
  // approximated position of the decay vertex
  vertex = make_unique<Point>((Lmin+Lmax)/2. * cos(track.GetPhi())    + 10*_eventVertex.GetX(),
                              (Lmin+Lmax)/2. * sin(track.GetPhi())    + 10*_eventVertex.GetY(),
                              (Lmin+Lmax)/2. / tan(track.GetTheta())  + 10*_eventVertex.GetZ(),
                              0.0, "",
                              fabs( (Lmax-Lmin)/2. * cos(track.GetPhi()) ),
                              fabs( (Lmax-Lmin)/2. * sin(track.GetPhi()) ),
                              fabs( (Lmax-Lmin)/2. / tan(track.GetTheta()) )
                              );

  points.push_back(make_shared<Point>(*vertex));
  points.push_back(make_shared<Point>(p1));
  points.push_back(make_shared<Point>(p2));
  
  auto circle = circleProcessor->GetCircleFromTriplet(points);
  
  origin = circle->GetCenter();
  
  double t0 = atan2(vertex->GetY() - origin.GetY(), vertex->GetX() - origin.GetX());
  
  double t1 = atan2(p1.GetY() - origin.GetY(), p1.GetX() - origin.GetX());
  while(t1 < t0) t1 += 2*TMath::Pi();
  
  double t2 = atan2(p2.GetY() - origin.GetY(), p2.GetX() - origin.GetX());
  while(t2 < t1) t2 += 2*TMath::Pi();
  
  tShift = t0;
  
  // this only gives approximate direction of the momentum vector
  momentum = make_unique<Point>(p1.GetX() - vertex->GetX(),
                                p1.GetY() - vertex->GetY(),
                                p1.GetZ() - vertex->GetZ());
  
  charge = track.GetCharge();
  
  double z0_min = vertex->GetZ()-vertex->GetZerr() - config.helixThickness;
  double z0_max = vertex->GetZ()+vertex->GetZerr() + config.helixThickness;
  
  double z1_min = p1.GetZ()-p1.GetZerr() - config.helixThickness;
  double z1_max = p1.GetZ()+p1.GetZerr() + config.helixThickness;
  
  double z2_min = p2.GetZ()-p2.GetZerr() - config.helixThickness;
  double z2_max = p2.GetZ()+p2.GetZerr() + config.helixThickness;
  
  // calculate all possible pairs of b -- s0 values and find extreme ones:
  CalcAndUpdateSlopeVars(z0_min, t0, z1_min, t1, z2_min, t2);
  CalcAndUpdateSlopeVars(z0_max, t0, z1_min, t1, z2_min, t2);
  CalcAndUpdateSlopeVars(z0_min, t0, z1_max, t1, z2_min, t2);
  CalcAndUpdateSlopeVars(z0_min, t0, z1_min, t1, z2_max, t2);
  CalcAndUpdateSlopeVars(z0_max, t0, z1_max, t1, z2_max, t2);
  CalcAndUpdateSlopeVars(z0_min, t0, z1_max, t1, z2_max, t2);
  CalcAndUpdateSlopeVars(z0_max, t0, z1_min, t1, z2_max, t2);
  CalcAndUpdateSlopeVars(z0_max, t0, z1_max, t1, z2_min, t2);

  // make sure that ranges are in the correct order:
  if(s0min > s0max) swap(s0min, s0max);
  if(bmin > bmax) swap(bmin, bmax);
  
  origin.SetZ((z0_min+z0_max)/2. - GetSlope(t0)*t0);
  
  double x0_min = vertex->GetX() - vertex->GetXerr() - config.helixThickness;
  double x0_max = vertex->GetX() + vertex->GetXerr() + config.helixThickness;
  
  double x1_min = p1.GetX()-p1.GetXerr()-config.helixThickness;
  double x1_max = p1.GetX()+p1.GetXerr()+config.helixThickness;
  
  double x2_min = p2.GetX()-p2.GetXerr()-config.helixThickness;
  double x2_max = p2.GetX()+p2.GetXerr()+config.helixThickness;
  
  CalcAndUpdateRadiiVars(x0_min, t0, x1_min, t1, x2_min, t2);
  CalcAndUpdateRadiiVars(x0_max, t0, x1_min, t1, x2_min, t2);
  CalcAndUpdateRadiiVars(x0_min, t0, x1_max, t1, x2_min, t2);
  CalcAndUpdateRadiiVars(x0_min, t0, x1_min, t1, x2_max, t2);
  CalcAndUpdateRadiiVars(x0_max, t0, x1_max, t1, x2_max, t2);
  CalcAndUpdateRadiiVars(x0_min, t0, x1_max, t1, x2_max, t2);
  CalcAndUpdateRadiiVars(x0_max, t0, x1_min, t1, x2_max, t2);
  CalcAndUpdateRadiiVars(x0_max, t0, x1_max, t1, x2_min, t2);
  
  if(R0min > R0max) swap(R0min, R0max);
  if(amin > amax) swap(amin, amax);
  
  tMax = t2;
  tStep = 0.01;
  iCycles=0;
}

pair<double, double> Helix::CalcSlopeVars(double z0, double t0, double z1, double t1, double z2, double t2)
{
  double b_num = (z0-z2)*(t0-t1) - (z0-z1)*(t0-t2);
  double b_den = (t2*t2-t0*t0)*(t0-t1) - (t1*t1-t0*t0)*(t0-t2);
  double b = b_num/b_den;
  double s0 = ( z0 - z1 - b*(t1*t1-t0*t0)) / (t0-t1);
  
  return make_pair(s0, b);
}

pair<double, double> Helix::CalcRadiiVars(double x0, double t0, double x1, double t1, double x2, double t2)
{
  double a_num = (x0-x2)*(cos(t0)-cos(t1)) - (x0-x1)*(cos(t0)-cos(t2));
  double a_den = (cos(t0)-cos(t1))*(t2*cos(t2)-t0*cos(t0)) - (cos(t0)-cos(t2))*(t1*cos(t1)-t0*cos(t0));
  double a = a_num/a_den;
  double R0 = ( x0 - x1 - a*(t1*cos(t1)-t0*cos(t0))) / (cos(t0)-cos(t1));
  
  return make_pair(R0, a);
}

void Helix::CalcAndUpdateSlopeVars(double z0, double t0, double z1, double t1, double z2, double t2)
{
  const auto &[s0, b] = CalcSlopeVars(z0, t0, z1, t1, z2, t2);
  
  double testT = 2*TMath::Pi();
  double val = s0 - b*testT;
  
  if(val < slope_valmin){ slope_valmin = val; s0min = s0; bmin = b; }
  if(val > slope_valmax){ slope_valmax = val; s0max = s0; bmax = b; }
}

void Helix::CalcAndUpdateRadiiVars(double x0, double t0, double x1, double t1, double x2, double t2)
{
  const auto &[R0, a] = CalcRadiiVars(x0, t0, x1, t1, x2, t2);
  
  double testT = 2*TMath::Pi();
  double val = R0 - a*testT;
  
  if(val < radius_valmin){ radius_valmin = val; R0min = R0; amin = a; }
  if(val > radius_valmax){ radius_valmax = val; R0max = R0; amax = a; }
}



bool Helix::ExtendByPoint(const Point &point)
{
  // make sure that this point is not yet on the helix
  for(auto &p : points){
    if(*p==point) return false;
  }
  
  
  Point p0 = *points[0];
  Point p1 = *points[points.size()-1];
  Point p2 = point;
  
  double t0 = tShift;
  double t1 = tMax;
  double t2 = atan2(p2.GetY() - origin.GetY(), p2.GetX() - origin.GetX());
  
  while(t2 < t1) t2 += 2*TMath::Pi();

  tMax = t2;

  // vertex
  double z0_min = p0.GetZ()-p0.GetZerr() - config.helixThickness;
  double z0_max = p0.GetZ()+p0.GetZerr() + config.helixThickness;
  
  // last point on the helix
  double z1_min = p1.GetZ()-p1.GetZerr() - config.helixThickness;
  double z1_max = p1.GetZ()+p1.GetZerr() + config.helixThickness;
  
  // the new candidate point
  double z2_min = p2.GetZ()-p2.GetZerr() - config.helixThickness;
  double z2_max = p2.GetZ()+p2.GetZerr() + config.helixThickness;
  
  vector<pair<double, double>> s0_b;
  
  s0_b.push_back(CalcSlopeVars(z0_min, t0, z1_min, t1, z2_min, t2));
  s0_b.push_back(CalcSlopeVars(z0_max, t0, z1_min, t1, z2_min, t2));
  s0_b.push_back(CalcSlopeVars(z0_min, t0, z1_max, t1, z2_min, t2));
  s0_b.push_back(CalcSlopeVars(z0_min, t0, z1_min, t1, z2_max, t2));
  s0_b.push_back(CalcSlopeVars(z0_max, t0, z1_max, t1, z2_max, t2));
  s0_b.push_back(CalcSlopeVars(z0_min, t0, z1_max, t1, z2_max, t2));
  s0_b.push_back(CalcSlopeVars(z0_max, t0, z1_min, t1, z2_max, t2));
  s0_b.push_back(CalcSlopeVars(z0_max, t0, z1_max, t1, z2_min, t2));
  
  double testT = 2*TMath::Pi();
  
  double new_s0min = inf;
  double new_s0max = -inf;
  
  double new_bmin = inf;
  double new_bmax = -inf;
  
  double new_valmin = inf;
  double new_valmax = -inf;
  
  for(auto &[s0, b] : s0_b){
    double val = s0 - b*testT;
    
    if(val < new_valmin){ new_valmin = val; new_s0min = s0; new_bmin = b; }
    if(val > new_valmax){ new_valmax = val; new_s0max = s0; new_bmax = b; }
  }
  
  if(new_s0min > new_s0max) swap(new_s0min, new_s0max);
  if(new_bmin > new_bmax) swap(new_bmin, new_bmax);
  
  // if there's no overlap between new limits and the current ones:
  if(new_s0max < s0min || new_s0min > s0max ||
     new_bmax < bmin || new_bmin > bmax
     ){
    return false;
  }
  
  if(new_s0max < s0max) s0max = new_s0max;
  if(new_s0min > s0min) s0min = new_s0min;
  
  if(new_bmax < bmax) bmax = new_bmax;
  if(new_bmin > bmin) bmin = new_bmin;
  
  // vertex
  double x0_min = p0.GetX()-p0.GetXerr() - config.helixThickness;
  double x0_max = p0.GetX()+p0.GetXerr() + config.helixThickness;
  
  // last point on the helix
  double x1_min = p1.GetX()-p1.GetXerr() - config.helixThickness;
  double x1_max = p1.GetX()+p1.GetXerr() + config.helixThickness;
  
  // the new candidate point
  double x2_min = p2.GetX()-p2.GetXerr() - config.helixThickness;
  double x2_max = p2.GetX()+p2.GetXerr() + config.helixThickness;
  
  vector<pair<double, double>> R0_a;
  
  R0_a.push_back(CalcRadiiVars(x0_min, t0, x1_min, t1, x2_min, t2));
  R0_a.push_back(CalcRadiiVars(x0_max, t0, x1_min, t1, x2_min, t2));
  R0_a.push_back(CalcRadiiVars(x0_min, t0, x1_max, t1, x2_min, t2));
  R0_a.push_back(CalcRadiiVars(x0_min, t0, x1_min, t1, x2_max, t2));
  R0_a.push_back(CalcRadiiVars(x0_max, t0, x1_max, t1, x2_max, t2));
  R0_a.push_back(CalcRadiiVars(x0_min, t0, x1_max, t1, x2_max, t2));
  R0_a.push_back(CalcRadiiVars(x0_max, t0, x1_min, t1, x2_max, t2));
  R0_a.push_back(CalcRadiiVars(x0_max, t0, x1_max, t1, x2_min, t2));
  
  double new_R0min = inf;
  double new_R0max = -inf;
  
  double new_amin = inf;
  double new_amax = -inf;
  
  double new_radius_valmin = inf;
  double new_radius_valmax = -inf;
  
  testT = 2*TMath::Pi();
  
  for(auto &[R0, a] : R0_a){
    double val = R0 - a*testT;
    
    if(val < new_radius_valmin){ new_radius_valmin = val; new_R0min = R0; new_amin = a; }
    if(val > new_radius_valmax){ new_radius_valmax = val; new_R0max = R0; new_amax = a; }
  }
  
  if(new_R0min > new_R0max) swap(new_R0min, new_R0max);
  if(new_amin > new_amax) swap(new_amin, new_amax);
  
  // if there's no overlap between new limits and the current ones:
  if(new_R0max < R0min || new_R0min > R0max ||
     new_amax < amin || new_amin > amax
     ){
    return false;
  }
  
  if(new_R0max < R0max) R0max = new_R0max;
  if(new_R0min > R0min) R0min = new_R0min;
  
  if(new_amax < amax) amax = new_amax;
  if(new_amin > amin) amin = new_amin;
  
  // adjust position of the origin
  double y0_min = p0.GetY()-p0.GetYerr() - config.helixThickness;
  double y0_max = p0.GetY()+p0.GetYerr() + config.helixThickness;
  
  origin.SetX((x0_min+x0_max)/2. - GetRadius(t0)*cos(t0));
  origin.SetY((y0_min+y0_max)/2. - GetRadius(t0)*sin(t0));
  origin.SetZ((z0_min+z0_max)/2. - GetSlope(t0)*t0);
  
  points.push_back(make_shared<Point>(point));

  return true;
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
  cout<<"\tL:"<<Lmin<<" -- "<<Lmax<<endl;
  cout<<"\tb:"<<bmin<<" -- "<<bmax<<endl;
  cout<<"\ts0:"<<s0min<<" -- "<<s0max<<endl;
  cout<<"\ta:"<<amin<<" -- "<<amax<<endl;
  cout<<"\tR:"<<R0min<<" -- "<<R0max<<endl;
  cout<<"\tt min:"<<tShift<<"\tt max:"<<tMax<<endl;
}

void Helix::SetPoints(const vector<shared_ptr<Point>> &_points)
{
  nPionPoints = 0;
  points.clear();
  
  for(auto &p : _points){
    Point q = GetClosestPoint(*p);
    if(pointsProcessor->distance(*p,q) < config.helixThickness){
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
