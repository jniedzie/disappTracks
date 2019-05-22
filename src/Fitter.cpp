//
//  Fitter.cpp
//
//  Created by Jeremi Niedziela on 14/12/2018.
//

#include "Fitter.hpp"

Fitter::Fitter() :
eventVertex(Point(0, 0, 0))
{

}

Fitter::~Fitter()
{
  
}

vector<Helix> Fitter::FitHelices(const vector<shared_ptr<Point>> &_points,
                                 const Track &_track,
                                 const Point &_eventVertex)
{
  points      = _points;
  track       = _track;
  eventVertex = _eventVertex;
  
  double maxChi2 = 1.0;
  
  double middleSeedHitMaxDeltaPhi = 1.5;
  double middleSeedHitMaxDeltaZ = 200;
  
  double lastSeedHitMaxDeltaPhi = 1.5;
  double lastSeedHitMaxDeltaZ = 100;
  
//  pointsProcessor.FilterNearbyPoints(points, 50);
  
  vector<vector<shared_ptr<Point>>> pointsByLayer = pointsProcessor.SortByLayer(points);
  
  // find possible middle and last seeds' points
  int trackLayers = track.GetNtrackerLayers();
  vector<shared_ptr<Point>> possibleMiddlePoints = pointsByLayer[trackLayers];
  vector<shared_ptr<Point>> possibleLastPoints   = pointsByLayer[trackLayers+1];
  
  // build seeds for all possible combinations of points
  vector<Helix> fittedHelices;
  
  cout<<"Fitter -- looking for seeds"<<endl;
  
  Point trackPointMin = pointsProcessor.GetPointOnTrack(layerR[track.GetNtrackerLayers()-1],  track, eventVertex);
  Point trackPointMax = pointsProcessor.GetPointOnTrack(layerR[track.GetNtrackerLayers()],    track, eventVertex);
  
  // Count how many seeds we'll need to test
  int nPairs=0;
  for(auto &middlePoint : possibleMiddlePoints){
    
    // Make sure that middle hit can be reached from track without crossing other layers of tracker (assuming straight track!!)
    double middleHitDeltaPhi = min(pointsProcessor.GetPointingAngleXY(Point(0,0,0), trackPointMin, *middlePoint),
                                   pointsProcessor.GetPointingAngleXY(Point(0,0,0), trackPointMax, *middlePoint));
    
    if(middleHitDeltaPhi > middleSeedHitMaxDeltaPhi) continue;
    
    double middleHitDeltaZ = min( fabs(middlePoint->GetZ() - trackPointMin.GetZ()),
                                  fabs(middlePoint->GetZ() - trackPointMax.GetZ()));
    
    if(middleHitDeltaZ > middleSeedHitMaxDeltaZ) continue;
    
    for(auto &lastPoint : possibleLastPoints){
      double lastHitDeltaPhi = min(pointsProcessor.GetPointingAngleXY(trackPointMin, *middlePoint, *lastPoint),
                                   pointsProcessor.GetPointingAngleXY(trackPointMax, *middlePoint, *lastPoint));
      
      if(lastHitDeltaPhi > lastSeedHitMaxDeltaPhi) continue;
      
      double lastPointDeltaZ = fabs(middlePoint->GetZ() - lastPoint->GetZ());
      
      if(lastPointDeltaZ > lastSeedHitMaxDeltaZ) continue;
      
      nPairs++;
    }
  }
  cout<<"N pairs to test:"<<nPairs<<endl;
  int nTested=0;
  for(auto &middlePoint : possibleMiddlePoints){
    
    double middleHitDeltaPhi = min(pointsProcessor.GetPointingAngleXY(Point(0,0,0), trackPointMin, *middlePoint),
                                   pointsProcessor.GetPointingAngleXY(Point(0,0,0), trackPointMax, *middlePoint));
    
    if(middleHitDeltaPhi > middleSeedHitMaxDeltaPhi) continue;
    
    double middleHitDeltaZ = min( fabs(middlePoint->GetZ() - trackPointMin.GetZ()),
                                 fabs(middlePoint->GetZ() - trackPointMax.GetZ()));
    
    if(middleHitDeltaZ > middleSeedHitMaxDeltaZ) continue;
    
    for(auto &lastPoint : possibleLastPoints){
      double lastHitDeltaPhi = min(pointsProcessor.GetPointingAngleXY(trackPointMin, *middlePoint, *lastPoint),
                                   pointsProcessor.GetPointingAngleXY(trackPointMax, *middlePoint, *lastPoint));
      
      if(lastHitDeltaPhi > lastSeedHitMaxDeltaPhi) continue;
      
      double lastPointDeltaZ = fabs(middlePoint->GetZ() - lastPoint->GetZ());
      
      if(lastPointDeltaZ > lastSeedHitMaxDeltaZ) continue;
      
      vector<shared_ptr<Point>> points = { middlePoint, lastPoint };
      unique_ptr<Helix> helixPos = FitSeed(points,  1);
      unique_ptr<Helix> helixNeg = FitSeed(points, -1);
      
      vector<unique_ptr<Helix>> goodHelices;
      
      if(helixPos){
        if(helixPos->chi2 < maxChi2) goodHelices.push_back(move(helixPos));
      }
      if(helixNeg){
        if(helixNeg->chi2 < maxChi2) goodHelices.push_back(move(helixNeg));
      }
      
      for(auto &helix : goodHelices){
        helix->increasing = true; // add decreasing later
        fittedHelices.push_back(*helix);
      }
      nTested++;
      if(nTested%100 == 0) cout<<"Tested: "<<nTested<<"/"<<nPairs<<endl;
    }
  }
  cout<<"Fitter -- found "<<fittedHelices.size()<<" valid seeds"<<endl;
  cout<<"Fitter -- extending seeds"<<endl;
  ExtendSeeds(fittedHelices, pointsByLayer, maxChi2);
  cout<<"Fitter -- merging overlapping helices. Before: "<<fittedHelices.size();
  MergeHelices(fittedHelices);
  cout<<"\tafter:"<<fittedHelices.size()<<endl;
  
  vector<Helix> longHelices;
  
  // Remove helices that even after merging have only 3 hits
  cout<<"Fitter -- removing very short helices"<<endl;
  for(int iHelix=0; iHelix<fittedHelices.size(); iHelix++){
    if(fittedHelices[iHelix].GetPoints().size() > 4){
      longHelices.push_back(fittedHelices[iHelix]);
    }
  }
  
  cout<<"Fitter -- refitting surviving helices ("<<longHelices.size()<<")"<<endl;
  for(auto &helix : longHelices){
    if(helix.shouldRefit) RefitHelix(helix);
  }
  
  cout<<"Fitter -- done"<<endl;
  return longHelices;
}

unique_ptr<Helix> Fitter::FitSeed(const vector<shared_ptr<Point>> &points, int charge)
{
  double Lmin = layerR[track.GetNtrackerLayers()-1];
  double Lmax = layerR[track.GetNtrackerLayers()];
  
  auto fitter = GetSeedFitter(range<double>(Lmin, Lmax));
  
  auto chi2Function = [&](const double *par) {
    double R0 = par[0];
    double a  = par[1];
    double s0 = par[2];
    double b  = par[3];
    double L  = par[4];
    double x0 = par[5];
    double y0 = par[6];
    double z0 = par[7];
    
    // First add distance to the vertex
    auto vertex = make_shared<Point>(pointsProcessor.GetPointOnTrack(L, track, eventVertex));
    auto pointsTriplet = {vertex, points[0], points[1] };
    
    double t, distX, distY, distZ;
    double f=0;
    
    // Then add distances to all other points
    for(auto &p : pointsTriplet){
      if(charge < 0) t = atan2(p->GetY() - y0, p->GetX() - x0);
      else           t = atan2(p->GetX() - x0, p->GetY() - y0);
      
      
      // find helix point for this point's t
      double x,y,z;
      if(charge < 0){
        x = x0 + (R0 - a*t)*cos(t);
        y = y0 + (R0 - a*t)*sin(t);
        z = z0 + (s0 - b*t)*t;
      }
      else{
        x = x0 + (R0 - a*t)*sin(t);
        y = y0 + (R0 - a*t)*cos(t);
        z = z0 + (s0 - b*t)*t;
      }
      
      // calculate distance between helix and point's boundary (taking into account its errors)
      distX = pow(x-p->GetX(), 2);
      distY = pow(y-p->GetY(), 2);
      distZ = pow(z-p->GetZ(), 2);
      
      distX /= p->GetXerr() > 0 ? pow(p->GetXerr(), 2) : fabs(p->GetX());
      distY /= p->GetYerr() > 0 ? pow(p->GetYerr(), 2) : fabs(p->GetY());
      distZ /= p->GetZerr() > 0 ? pow(p->GetZerr(), 2) : fabs(p->GetZ());
      
      f += distX + distY + distZ;
    }
    return f;
  };
  

  int nPar=8;
  auto fitFunction = ROOT::Math::Functor(chi2Function, nPar);
  double pStart[nPar];
  for(int i=0; i<nPar; i++) pStart[i] = fitter->Config().ParSettings(i).Value();
  fitter->SetFCN(fitFunction, pStart);
  
  unique_ptr<Helix> resultHelix = nullptr;
  
  if(fitter->FitFCN()) {
    auto result = fitter->Result();
    
    // Build helix from fitters output
    HelixParams resultParams;
    resultParams.R0 = result.GetParams()[0];
    resultParams.a  = result.GetParams()[1];
    resultParams.s0 = result.GetParams()[2];
    resultParams.b  = result.GetParams()[3];
    
    double L  = result.GetParams()[4];
    double x0 = result.GetParams()[5];
    double y0 = result.GetParams()[6];
    double z0 = result.GetParams()[7];
    
    Point origin(x0, y0, z0);
    Point vertex = pointsProcessor.GetPointOnTrack(L, track, eventVertex);
    if(charge < 0)  vertex.SetT(atan2(vertex.GetY() - origin.GetY(), vertex.GetX() - origin.GetX()));
    else            vertex.SetT(atan2(vertex.GetX() - origin.GetX(), vertex.GetY() - origin.GetY()));
    
    for(auto &p : points){
      if(charge < 0)  p->SetT(atan2(p->GetY() - y0, p->GetX() - x0));
      else            p->SetT(atan2(p->GetX() - x0, p->GetY() - y0));
    }
    
    // Check if ordering of the points is correct (third point in a cone relative to vertex+first point)
    double phi = pointsProcessor.GetPointingAngle(vertex, *points[0], *points[1]);
    
    if(phi < TMath::Pi()/2.){
      resultHelix = make_unique<Helix>(resultParams, vertex, origin, points, track);
      resultHelix->chi2 = result.MinFcnValue();
      resultHelix->SetCharge(charge);
    }
  }
  
  return resultHelix;
}

void Fitter::ExtendSeeds(vector<Helix> &helices,
                         const vector<vector<shared_ptr<Point>>> &pointsByLayer,
                         double maxChi2)
{
  double nextPointMaxDeltaPhi = 0.7;
  double nextPointMaxDeltaZ = 200;
  
  bool finished;
  int nSteps=0;
  do{
    cout<<"Performing step "<<nSteps<<endl;
    cout<<"Helices to extend: "<<helices.size()<<endl;
    finished = true;
    vector<Helix> nextStepHelices;
    
    // for all helices from previous step
    for(Helix &helix : helices){
      
      if(helix.isFinished){
        nextStepHelices.push_back(helix);
      }
      else{
        // Find points that could extend this helix
        int lastPointLayer = helix.GetLastPoint()->GetLayer();
        vector<shared_ptr<Point>> possiblePoints;
        
        if(lastPointLayer+1 < pointsByLayer.size() && lastPointLayer-1 >= 0){
          if(helix.increasing)  possiblePoints = pointsByLayer[lastPointLayer+1];
          else                  possiblePoints = pointsByLayer[lastPointLayer-1];
        }
        
        vector<Helix> extendedHelices;
        
        // try to extend by all possible points
        for(auto &point : possiblePoints){
          
          // Check if new point is within some cone
          size_t nHelixPoints = helix.GetPoints().size();
          double deltaPhi = pointsProcessor.GetPointingAngleXY(*helix.GetPoints()[nHelixPoints-2],
                                                               *helix.GetPoints()[nHelixPoints-1],
                                                               *point);
          
          if(deltaPhi > nextPointMaxDeltaPhi) continue;
          
          double deltaZ   = fabs(helix.GetPoints()[nHelixPoints-1]->GetZ() - point->GetZ());
          
          if(deltaZ > nextPointMaxDeltaZ) continue;
          
          /// Extend helix by the new point and refit its params
          Helix helixCopy(helix);
          helixCopy.AddPoint(point);
          RefitHelix(helixCopy);
          
          // check if chi2 is small enough
          if(helixCopy.chi2 > maxChi2) continue;
          
          extendedHelices.push_back(helixCopy);
        }
        // if it was possible to extend the helix
        if(extendedHelices.size() != 0){
          nextStepHelices.insert(nextStepHelices.end(),
                                 extendedHelices.begin(),
                                 extendedHelices.end());
          finished = false;
        }
        else if(!helix.isFinished){ // if helix could not be extended
          helix.isFinished = true;
          nextStepHelices.push_back(helix);
        }
      }
    }
    
    helices.clear();
    for(auto &h : nextStepHelices) helices.push_back(h);
    /*
    for(auto &h : nextStepHelices){
      if(nSteps+3 <= minNhits){
        if(h.GetNpoints() >= nSteps+3) helices.push_back(h);
      }
      else{
        if(h.GetNpoints() >= minNhits) helices.push_back(h);
      }
    }
    */
    nSteps++;
  }
  while(!finished);
}

void Fitter::MergeHelices(vector<Helix> &helices)
{
  // Merge helices that are very similar to each other
  bool merged=true;
  
  while(merged){
    merged=false;
    
    for(int iHelix1=0; iHelix1<helices.size(); iHelix1++){
      Helix &helix1 = helices[iHelix1];
      
      for(int iHelix2=iHelix1+1; iHelix2<helices.size(); iHelix2++){
        Helix &helix2 = helices[iHelix2];
        
        vector<shared_ptr<Point>> points1 = helix1.GetPoints();
        vector<shared_ptr<Point>> points2 = helix2.GetPoints();
        
        sort(points1.begin(), points1.end());
        sort(points2.begin(), points2.end());
        
        vector<shared_ptr<Point>> samePoints;
        set_intersection(points1.begin(), points1.end(),
                         points2.begin(), points2.end(),
                         back_inserter(samePoints));
        
        double samePointsFraction = samePoints.size()/(double)points1.size();

        // merging
        if(samePointsFraction > 0.7){
          // remove second helix
          helices.erase(helices.begin() + iHelix2);
          
          // update first helix
          unordered_set<shared_ptr<Point>> uniquePoints;
          for(auto &p : points1) uniquePoints.insert(p);
          for(auto &p : points2) uniquePoints.insert(p);
          vector<shared_ptr<Point>> allPoints(uniquePoints.begin(), uniquePoints.end());
          helix1.ReplacePoints(allPoints);
          helix1.shouldRefit = true;
          
          merged=true;
          break;
        }
      }
      
      if(merged) break;
    }
  }
}

void Fitter::RefitHelix(Helix &helix)
{
  vector<shared_ptr<Point>> helixPoints = helix.GetPoints();
  double ht = config.helixThickness;
  
  auto chi2Function = [&](const double *par) {
    double R0 = par[0];
    double a  = par[1];
    double s0 = par[2];
    double b  = par[3];
    double L  = par[4];
    double x0 = par[5];
    double y0 = par[6];
    double z0 = par[7];
    
    // First add distance to the new vertex
    auto vertex = make_shared<Point>(pointsProcessor.GetPointOnTrack(L, track, eventVertex));
    auto points = helixPoints;
    points.push_back(vertex);
    double t, distX, distY, distZ;
    double f=0;
    
    // Then add distances to all other points
//    bool first=true;
    for(auto &p : points){
      double x, y, z;
      
      // find helix point for this point's t
      if(helix.GetCharge() < 0){
        t = atan2(p->GetY() - y0, p->GetX() - x0);
        x = x0 + (R0 - a*t)*cos(t);
        y = y0 + (R0 - a*t)*sin(t);
        z = z0 + (s0 - b*t)*t;
      }
      else{
        t = atan2(p->GetX() - x0, p->GetY() - y0);
        x = x0 + (R0 - a*t)*sin(t);
        y = y0 + (R0 - a*t)*cos(t);
        z = z0 + (s0 - b*t)*t;
      }
      
      // calculate distance between helix and point's boundary (taking into account its errors)
      distX = distY = distZ = 0;
      
      if(fabs(x-p->GetX()) > p->GetXerr()+ht){
        double distX_1 = x - (p->GetX() + p->GetXerr() + ht);
        double distX_2 = x - (p->GetX() - p->GetXerr() - ht);
        distX = min(pow(distX_1, 2), pow(distX_2, 2));
      }
      if(fabs(y-p->GetY()) > p->GetYerr()+ht){
        double distY_1 = y - (p->GetY() + p->GetYerr() + ht);
        double distY_2 = y - (p->GetY() - p->GetYerr() - ht);
        distY = min(pow(distY_1, 2), pow(distY_2, 2));
      }
      if(fabs(z-p->GetZ()) > p->GetZerr()+ht){
        double distZ_1 = z - (p->GetZ() + p->GetZerr() + ht);
        double distZ_2 = z - (p->GetZ() - p->GetZerr() - ht);
        distZ = min(pow(distZ_1, 2), pow(distZ_2, 2));
      }
      
      distX /= p->GetXerr() > 0 ? pow(p->GetXerr(), 2) : fabs(p->GetX());
      distY /= p->GetYerr() > 0 ? pow(p->GetYerr(), 2) : fabs(p->GetY());
      distZ /= p->GetZerr() > 0 ? pow(p->GetZerr(), 2) : fabs(p->GetZ());
      
      f += distX + distY + distZ;
    }
    
    return f/(3*helixPoints.size()+8);
  };
  
  ROOT::Fit::Fitter *fitter = new ROOT::Fit::Fitter();
  int nPar = 8;
  ROOT::Math::Functor fitFunction = ROOT::Math::Functor(chi2Function, nPar);
  double pStart[nPar];
  fitter->SetFCN(fitFunction, pStart);
  
  double Lmin = layerR[track.GetNtrackerLayers()-1];
  double Lmax = layerR[track.GetNtrackerLayers()];
  
  SetParameter(fitter, 0, "R0", helix.helixParams.R0,
               GetRadiusInMagField(config.minPx, config.minPy, solenoidField),
               GetRadiusInMagField(config.maxPx, config.maxPy, solenoidField));
  SetParameter(fitter, 1, "a" , helix.helixParams.a, 0, 10000);
  SetParameter(fitter, 2, "s0", helix.helixParams.s0, -1000, 1000);
  SetParameter(fitter, 3, "b" , helix.helixParams.b, -10000, 0);
  SetParameter(fitter, 4, "L" , (helix.GetVertex()->GetX()-10*eventVertex.GetX())/cos(track.GetPhi()),
               Lmin, Lmax);
  SetParameter(fitter, 5, "x0", helix.GetOrigin().GetX(), -layerR[nLayers-1] , layerR[nLayers-1]);
  SetParameter(fitter, 6, "y0", helix.GetOrigin().GetY(), -layerR[nLayers-1] , layerR[nLayers-1]);
  SetParameter(fitter, 7, "z0", helix.GetOrigin().GetZ(), -pixelBarrelZsize  , pixelBarrelZsize);
  
  for(int i=0; i<nPar; i++) pStart[i] = fitter->Config().ParSettings(i).Value();
  
//  if(fitter->FitFCN()) {
  
  helix.SetCharge(1);
  fitter->FitFCN();
  auto resultPos = fitter->Result();
  double chi2pos = resultPos.MinFcnValue();
  
  helix.SetCharge(-1);
  fitter->FitFCN();
  auto resultNeg = fitter->Result();
  double chi2neg = resultNeg.MinFcnValue();
  
  auto result = (chi2pos < chi2neg) ? resultPos : resultNeg;
  
  if(chi2pos < chi2neg){
    helix.SetCharge(1);
    helix.chi2 = chi2pos;
  }
  else{
    helix.SetCharge(-1);
    helix.chi2 = chi2neg;
  }
  
  double L  = result.GetParams()[4];
  Point vertex = pointsProcessor.GetPointOnTrack(L, track, eventVertex);
  
  helix.SetVertex(vertex);
  
  double x0 = result.GetParams()[5];
  double y0 = result.GetParams()[6];
  double z0 = result.GetParams()[7];
  Point origin(x0, y0, z0);
  helix.UpdateOrigin(origin);
  
  HelixParams resultParams;
  resultParams.R0 = result.GetParams()[0];
  resultParams.a  = result.GetParams()[1];
  resultParams.s0 = result.GetParams()[2];
  resultParams.b  = result.GetParams()[3];
  helix.helixParams = resultParams;
  
}

ROOT::Fit::Fitter* Fitter::GetSeedFitter(range<double> rangeL)
{
  ROOT::Fit::Fitter *fitter = new ROOT::Fit::Fitter();
  
  // This is a stupid hack to be able to set fit parameters before actually setting a fitting function
  // Params we want to set only once, but function will change in each iteration of the loop, because
  // it captures loop iterators.
  auto f = [&](const double*) {return 0;};
  int nPar = 8;
  ROOT::Math::Functor fitFunction = ROOT::Math::Functor(f, nPar);
  double pStart[nPar];
  fitter->SetFCN(fitFunction, pStart);
  
  double minR = GetRadiusInMagField(config.minPx, config.minPy, solenoidField);
  double maxR = GetRadiusInMagField(config.maxPx, config.maxPy, solenoidField);
  
  SetParameter(fitter, 0, "R0", (minR+maxR)/2., 0, 10000);
  SetParameter(fitter, 1, "a" , (minR+maxR)/2., 0, 10000);
  SetParameter(fitter, 2, "s0", 0, -1000, 1000);
  SetParameter(fitter, 3, "b" , 0, -10000, 0);
  SetParameter(fitter, 4, "L" , (rangeL.GetMin()+rangeL.GetMax())/2., rangeL.GetMin(), rangeL.GetMax());
  SetParameter(fitter, 5, "x0", 0, -layerR[nLayers-1] , layerR[nLayers-1]);
  SetParameter(fitter, 6, "y0", 0, -layerR[nLayers-1] , layerR[nLayers-1]);
  SetParameter(fitter, 7, "z0", 0, -pixelBarrelZsize  , pixelBarrelZsize);
  
  return fitter;
}

void Fitter::SetParameter(ROOT::Fit::Fitter *fitter, int i, string name, double start, double min, double max, bool fix)
{
  fitter->Config().ParSettings(i).SetName(name);
  fitter->Config().ParSettings(i).SetValue(start);
  fitter->Config().ParSettings(i).SetLimits((min < max) ? min : max,(min < max) ? max : min);
  fitter->Config().ParSettings(i).SetStepSize(0.0001);
  if(fix) fitter->Config().ParSettings(i).Fix();
}

void Fitter::FixParameter(ROOT::Fit::Fitter *fitter, int i, string name, double val)
{
  fitter->Config().ParSettings(i).SetName(name);
  fitter->Config().ParSettings(i).SetValue(val);
  fitter->Config().ParSettings(i).Fix();
}
