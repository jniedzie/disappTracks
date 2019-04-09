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
R0(h.R0),
//R0min(h.R0min),
//R0max(h.R0max),
iCycles(h.iCycles),
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

}

Helix::Helix(const Track &_track, const Point &p1, const Point &p2, const Point &_eventVertex) :
origin(Point(0,0,0)),
circleProcessor(make_unique<CircleProcessor>())
{
  track = _track;
  
  int nLayers = track.GetNtrackerLayers();
  
  Lmin = layerR[nLayers-1];
  Lmax = layerR[nLayers];
  
  // approximated position of the decay vertex
  vertex = make_unique<Point>((Lmin+Lmax)/2. * cos(track.GetPhi())    + 10*_eventVertex.GetX(),
                              (Lmin+Lmax)/2. * sin(track.GetPhi())    + 10*_eventVertex.GetY(),
                              (Lmin+Lmax)/2. / tan(track.GetTheta())  + 10*_eventVertex.GetZ());
  
  points.push_back(make_shared<Point>(*vertex));
  points.push_back(make_shared<Point>(p1));
  points.push_back(make_shared<Point>(p2));
  
  auto circle = circleProcessor->GetCircleFromTriplet(points);
  
  origin = circle->GetCenter();
  
  double t0 = atan2(vertex->GetY() - origin.GetY(), vertex->GetX() - origin.GetX());
  double t1 = atan2(p1.GetY() - origin.GetY(), p1.GetX() - origin.GetX());
  double t2 = atan2(p2.GetY() - origin.GetY(), p2.GetX() - origin.GetX());
  tShift = t0;
  
  // this only gives approximate direction of the momentum vector
  momentum = make_unique<Point>(p1.GetX() - vertex->GetX(),
                                p1.GetY() - vertex->GetY(),
                                p1.GetZ() - vertex->GetZ());
  
  charge = track.GetCharge();

  double z0_min = Lmin/tan(track.GetTheta()) + 10*_eventVertex.GetZ();
  double z0_max = Lmax/tan(track.GetTheta()) + 10*_eventVertex.GetZ();
  
  double z1_min = p1.GetZ()-p1.GetZerr();
  double z1_max = p1.GetZ()+p1.GetZerr();
  
  double z2_min = p2.GetZ()-p2.GetZerr();
  double z2_max = p2.GetZ()+p2.GetZerr();
  
  // calculate all possible pairs of b -- s0 values and find extreme ones:
  CalcAndUpdateSlopeVars(z0_min, t0, z1_min, t1, z2_min, t2);
  CalcAndUpdateSlopeVars(z0_max, t0, z1_min, t1, z2_min, t2);
  CalcAndUpdateSlopeVars(z0_min, t0, z1_max, t1, z2_min, t2);
  CalcAndUpdateSlopeVars(z0_min, t0, z1_min, t1, z2_max, t2);
  CalcAndUpdateSlopeVars(z0_max, t0, z1_max, t1, z2_max, t2);
  CalcAndUpdateSlopeVars(z0_min, t0, z1_max, t1, z2_max, t2);
  CalcAndUpdateSlopeVars(z0_max, t0, z1_min, t1, z2_max, t2);
  CalcAndUpdateSlopeVars(z0_max, t0, z1_max, t1, z2_min, t2);

  origin.SetZ((z0_min+z0_max)/2. - GetSlope(t0)*t0);
  
  // at the beginning, we cannot tell how much the radius will shrink and what are the limits
  // in the initial radius value
  // average should be zero/initial radius, but any value other than zero should improve those limits
  amin = -100000;
  amax =  100000;
  
  R0 = circle->GetRadius();
  
  tMax = t2;
  tStep = 0.005;
  iCycles=0;
}

void Helix::CalcAndUpdateSlopeVars(double z0, double t0, double z1, double t1, double z2, double t2)
{
  double b_num = (z0-z2)*(t0-t1) - (z0-z1)*(t0-t2);
  double b_den = (t2*t2-t0*t0)*(t0-t1) - (t1*t1-t0*t0)*(t0-t2);
  double b = b_num/b_den;
  double s0 = ( z0 - z1 - b*(t1*t1-t0*t0)) / (t0-t1);
  
  double testT = 2*TMath::Pi();
  double val = s0 - b*testT;
  
  if(val < valmin){ valmin = val; s0min = s0; bmin = b; }
  if(val > valmax){ valmax = val; s0max = s0; bmax = b; }
}

bool Helix::ExtendByPoint(const Point &point)
{
  // make sure that this point is not yet on the helix
  for(auto &p : points){
    if(*p==point) return false;
  }
  
  Point lastPoint = *points.back();
  
  double t0 = atan2(lastPoint.GetY() - origin.GetY(), lastPoint.GetX() - origin.GetX());
  double t  = atan2(point.GetY() - origin.GetY(), point.GetX() - origin.GetX());
  
  t0  += iCycles * 2*TMath::Pi();
  t   += iCycles * 2*TMath::Pi();
  
  if( t < t0){
    iCycles++;
    t += 2*TMath::Pi();
  }
  
  double b1 = (point.GetZ()+point.GetZerr() - Lmin/tan(track.GetTheta())) / (t0*t);
  double b2 = (point.GetZ()-point.GetZerr() - Lmin/tan(track.GetTheta())) / (t0*t);
  double b3 = (point.GetZ()+point.GetZerr() - Lmax/tan(track.GetTheta())) / (t0*t);
  double b4 = (point.GetZ()-point.GetZerr() - Lmax/tan(track.GetTheta())) / (t0*t);
  
  double bmin_new = min(min(min(b1, b2), b3), b4);
  double bmax_new = max(max(max(b1, b2), b3), b4);
  
  // slope cannot get larger with time:
  if(bmax_new < 0) bmax_new = 0;
  if(bmin_new < 0) bmin_new = 0;
  
  // if new value of b is outside of the allowed range, don't use this point
  if(bmax_new < bmin || bmin_new > bmax) return false;
  
  // if new point narrows down b range, update
  if(bmax_new < bmax) bmax = bmax_new;
  if(bmin_new > bmin) bmin = bmin_new;
  
  double s0_1 = (lastPoint.GetZ()+lastPoint.GetZerr() - (point.GetZ()+point.GetZerr()) + bmin*(t0*t0 - t*t)) / (t0 - t);
  double s0_2 = (lastPoint.GetZ()+lastPoint.GetZerr() - (point.GetZ()-point.GetZerr()) + bmin*(t0*t0 - t*t)) / (t0 - t);
  double s0_3 = (lastPoint.GetZ()-lastPoint.GetZerr() - (point.GetZ()+point.GetZerr()) + bmin*(t0*t0 - t*t)) / (t0 - t);
  double s0_4 = (lastPoint.GetZ()-lastPoint.GetZerr() - (point.GetZ()-point.GetZerr()) + bmin*(t0*t0 - t*t)) / (t0 - t);
  
  double s0_5 = (lastPoint.GetZ()+lastPoint.GetZerr() - (point.GetZ()+point.GetZerr()) + bmax*(t0*t0 - t*t)) / (t0 - t);
  double s0_6 = (lastPoint.GetZ()+lastPoint.GetZerr() - (point.GetZ()-point.GetZerr()) + bmax*(t0*t0 - t*t)) / (t0 - t);
  double s0_7 = (lastPoint.GetZ()-lastPoint.GetZerr() - (point.GetZ()+point.GetZerr()) + bmax*(t0*t0 - t*t)) / (t0 - t);
  double s0_8 = (lastPoint.GetZ()-lastPoint.GetZerr() - (point.GetZ()-point.GetZerr()) + bmax*(t0*t0 - t*t)) / (t0 - t);
  
  vector<double> s0 = { s0_1, s0_2, s0_3, s0_4, s0_5, s0_6, s0_7, s0_8 };
  
  double s0min_new = *min_element(s0.begin(), s0.end());
  double s0max_new = *max_element(s0.begin(), s0.end());
  
  // if new value of s0 is outside of the allowed range, don't use this point
  if(s0max_new < s0min || s0min_new > s0max) return false;
  
  // if new point narrows down s0 range, update
  if(s0max_new < s0max) s0max = s0max_new;
  if(s0min_new > s0min) s0min = s0min_new;
  
  double amin_new = ((lastPoint.GetX() - Lmin*cos(track.GetPhi()))*(cos(t0)-cos(t))) / (sin(t0)*cos(t)*(t0-t));
  double amax_new = ((lastPoint.GetX() - Lmax*cos(track.GetPhi()))*(cos(t0)-cos(t))) / (sin(t0)*cos(t)*(t0-t));
  
  if(amin_new > amax_new) swap(amin_new, amax_new);
  
  // radius cannot get larger with time
//  if(amax_new < 0) amax_new = 0;
//  if(amin_new < 0) amin_new = 0;
  
  // if new value of a is outside of the allowed range, don't use this point
  if(amax_new < amin || amin_new > amax) return false;
  
  // if new point narrows down a range, update
  if(amax_new < amax) amax = amax_new;
  if(amin_new > amin) amin = amin_new;
  
//  double R0min_new = ((lastPoint.GetX() - point.GetX()) + amin*(t0*cos(t0)-t*cos(t))) / (cos(t0)-cos(t));
//  double R0max_new = ((lastPoint.GetX() - point.GetX()) + amax*(t0*cos(t0)-t*cos(t))) / (cos(t0)-cos(t));
  
//  if(R0min_new > R0max_new) swap(R0min_new, R0max_new);
  
//  if(R0max_new < R0min || R0min_new > R0max) return false;
  
  // if new point narrows down R0 range, update
//  if(R0max_new < R0max) R0max = R0max_new;
//  if(R0min_new > R0min) R0min = R0min_new;
  
  points.push_back(make_shared<Point>(point));
  
  tMax = t;
  
  return true;
}

void Helix::Print()
{
  cout<<"\tVertex:("<<vertex->GetX()<<","<<vertex->GetY()<<","<<vertex->GetZ()<<")\n";
  cout<<"\tOrigin:("<<origin.GetX()<<","<<origin.GetY()<<","<<origin.GetZ()<<")\n";
  cout<<"\tMomentum:("<<momentum->GetX()<<","<<momentum->GetY()<<","<<momentum->GetZ()<<")\n";
  cout<<"\tCharge: "<<charge<<"\tR:"<<radius<<"\tc:"<<slope<<"\n";
  cout<<"\tnPoints:"<<points.size()<<"\tnPionPoints:"<<nPionPoints<<"\tnRegularPoints:"<<nRegularPoints<<"\n";
  cout<<"\tL:"<<Lmin<<" -- "<<Lmax<<endl;
  cout<<"\tb:"<<bmin<<" -- "<<bmax<<endl;
  cout<<"\ts0:"<<s0min<<" -- "<<s0max<<endl;
  cout<<"\ta:"<<amin<<" -- "<<amax<<endl;
  cout<<"\tR0:"<<R0<<endl;
//  cout<<"\tR:"<<R0min<<" -- "<<R0max<<endl;
  cout<<"\tt min:"<<tShift<<"\tt max:"<<tMax<<endl;
}

void Helix::SetPoints(const vector<shared_ptr<Point>> &_points)
{
  nPionPoints = 0;
  points.clear();
  
  for(auto &p : _points){
    Point q = GetClosestPoint(*p);
    if(pointsProcessor->distance(*p,q) < config->helixThickness){
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
