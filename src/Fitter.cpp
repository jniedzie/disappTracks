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

  vector<vector<shared_ptr<Point>>> pointsByLayer = pointsProcessor.SortByLayer(points);
  
  // Find seeds
  cout<<"Looking for seeds..."<<endl;
  vector<Helix> fittedHelices = GetSeeds(pointsByLayer);
  cout<<"Number of valid seeds: "<<fittedHelices.size()<<endl;

  int iter=0;
  for(auto &helix : fittedHelices){
    helix.uniqueID = iter++;
  }
  
  // Extend seeds
  cout<<"Extending seeds..."<<endl;
  ExtendSeeds(fittedHelices, pointsByLayer);
  cout<<"Candidates found: "<<fittedHelices.size()<<endl;
  
  // Merge similar candidates
  cout<<"Merging overlapping helices...";
  while(MergeHelices2(fittedHelices));
  cout<<" merged down to: "<<fittedHelices.size()<<endl;
  
  // Remove helices that are too short
  cout<<"Removing very short helices...";
  vector<Helix> longHelices;
  for(int iHelix=0; iHelix<fittedHelices.size(); iHelix++){
//    if(fittedHelices[iHelix].GetPoints().size() > 7){
      longHelices.push_back(fittedHelices[iHelix]);
//    }
  }
  cout<<" long helices: "<<longHelices.size()<<endl;
  
  cout<<"Refitting surviving helices...";
  for(auto &helix : longHelices){
    if(helix.shouldRefit) RefitHelix(helix);
  }
  cout<<" done."<<endl;
  
  return longHelices;
}

vector<Helix> Fitter::GetSeeds(vector<vector<shared_ptr<Point>>> pointsByLayer)
{
  double maxChi2 = 3.0;
  
  double middleSeedHitMaxDeltaPhi = 1.5;
  double middleSeedHitMaxDeltaZ = 200;
  
  double lastSeedHitMaxDeltaPhi = 1.5;
  double lastSeedHitMaxDeltaZ = 100;
  
  // find possible middle and last seeds' points
  int trackLayers = track.GetNtrackerLayers();
  cout<<"Track layers:"<<trackLayers<<endl;
  vector<shared_ptr<Point>> possibleMiddlePoints = pointsByLayer[trackLayers];
  vector<shared_ptr<Point>> possibleLastPoints   = pointsByLayer[trackLayers+1];
  
  vector<Helix> seeds;
  
  Point trackPointMin = pointsProcessor.GetPointOnTrack(layerR[track.GetNtrackerLayers()-1],  track, eventVertex);
  Point trackPointMax = pointsProcessor.GetPointOnTrack(layerR[track.GetNtrackerLayers()],    track, eventVertex);
  int nPairs=0;
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
      
      nPairs++;
      auto points = { middlePoint, lastPoint };
      auto helix = FitSeed(points,  track.GetCharge());
      
      if(helix){
        helix->increasing = true; // add decreasing later
        if(helix->chi2 < maxChi2) seeds.push_back(*helix);
        else{
          cout<<"Seed chi2 out of limits"<<endl;
        }
      }
    }
  }
  cout<<"Tested pairs: "<<nPairs<<endl;
  return seeds;
}

unique_ptr<Helix> Fitter::FitSeed(const vector<shared_ptr<Point>> &points, int charge)
{
  double Lmin = layerRanges[track.GetNtrackerLayers()-1].GetMax();
  double Lmax = layerRanges[track.GetNtrackerLayers()].GetMin();
  
  auto fitter = GetSeedFitter(range<double>(Lmin, Lmax));
  double ht = 0.0;//config.helixThickness;
  
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
//    cout<<endl;
    // Then add distances to all other points
    for(auto &p : pointsTriplet){
      if(charge < 0) t = atan2(p->GetY() - y0, p->GetX() - x0);
      else           t = TMath::Pi()/2 - atan2(p->GetX() - x0, p->GetY() - y0);
      
//      cout<<"t:"<<t<<endl;
      
      // find helix point for this point's t
      double x,y,z;
      x = x0 + (R0 - a*t)*cos(t);
      y = y0 + (R0 - a*t)*sin(t);
      z = z0 + (s0 - b*t)*t;
      
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
//      distX = pow(x-p->GetX(), 2);
//      distY = pow(y-p->GetY(), 2);
//      distZ = pow(z-p->GetZ(), 2);
      
      distX /= fabs(p->GetX());
      distY /= fabs(p->GetY());
      distZ /= fabs(p->GetZ());
      
//      distX /= p->GetXerr() > 0 ? pow(p->GetXerr(), 2) : fabs(p->GetX());
//      distY /= p->GetYerr() > 0 ? pow(p->GetYerr(), 2) : fabs(p->GetY());
//      distZ /= p->GetZerr() > 0 ? pow(p->GetZerr(), 2) : fabs(p->GetZ());
   
      f += distX + distY + distZ;
    }
    return f/(3*3+8);
//    return f;
  };
  

  int nPar=8;
  auto fitFunction = ROOT::Math::Functor(chi2Function, nPar);
  double pStart[nPar];
  for(int i=0; i<nPar; i++) pStart[i] = fitter->Config().ParSettings(i).Value();
  fitter->SetFCN(fitFunction, pStart);
  
  unique_ptr<Helix> resultHelix = nullptr;
  
  fitter->FitFCN();
  auto result = fitter->Result();
  
  HelixParams resultParams;
  resultParams.R0 = result.GetParams()[0];
  resultParams.a  = result.GetParams()[1];
  resultParams.s0 = result.GetParams()[2];
  resultParams.b  = result.GetParams()[3];
  
  double L  = result.GetParams()[4];
  Point vertex = pointsProcessor.GetPointOnTrack(L, track, eventVertex);
  
  double x0 = result.GetParams()[5];
  double y0 = result.GetParams()[6];
  double z0 = result.GetParams()[7];

  Point origin(x0, y0, z0);
  
  resultHelix = make_unique<Helix>(resultParams, vertex, origin, points, track);
  resultHelix->SetCharge(charge);
  resultHelix->UpdateOrigin(origin);
  resultHelix->chi2 = result.MinFcnValue();
  
  cout<<"Final chi2:"<<resultHelix->chi2<<endl;
  
  return resultHelix;
}

void Fitter::ExtendSeeds(vector<Helix> &helices,
                         const vector<vector<shared_ptr<Point>>> &pointsByLayer)
{
  double maxChi2 = 0.1;
  double nextPointMaxDeltaPhi = 0.7;
  double nextPointMaxDeltaZ = 200;
  
  bool finished;
  int nSteps=0;
  do{
    cout<<"Helices before "<<nSteps<<" step: "<<helices.size()<<endl;
    
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
          
          double deltaZ = fabs(helix.GetPoints()[nHelixPoints-1]->GetZ() - point->GetZ());
          
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
    nSteps++;
  }
  while(!finished);
}

void Fitter::MergeHelices(vector<Helix> &helices)
{
  // Merge helices that are very similar to each other
  bool merged = true;
  int nPasses = 0;
  
  while(merged){
    merged=false;
//    cout<<"Merging pass "<<nPasses<<", n helices: "<<helices.size()<<endl;
    
    for(int iHelix1=0; iHelix1<helices.size(); iHelix1++){
      Helix &helix1 = helices[iHelix1];
      vector<shared_ptr<Point>> points1 = helix1.GetPoints();
      sort(points1.begin(), points1.end());
      
      for(int iHelix2=iHelix1+1; iHelix2<helices.size(); iHelix2++){
        
        Helix &helix2 = helices[iHelix2];
        vector<shared_ptr<Point>> points2 = helix2.GetPoints();
        sort(points2.begin(), points2.end());
        
        vector<shared_ptr<Point>> samePoints;
        set_intersection(points1.begin(), points1.end(),
                         points2.begin(), points2.end(),
                         back_inserter(samePoints));
        
        double samePointsFraction = samePoints.size()/(double)points1.size();
        int nDifferentPoints = max(points1.size()-samePoints.size(),
                                   points2.size()-samePoints.size());
        
//        if(samePointsFraction > 0.8){
        if(nDifferentPoints < 3){
          // update first helix
          unordered_set<shared_ptr<Point>> uniquePoints;
          for(auto &p : points1) uniquePoints.insert(p);
          for(auto &p : points2) uniquePoints.insert(p);
          vector<shared_ptr<Point>> allPoints(uniquePoints.begin(), uniquePoints.end());
          helix1.ReplacePoints(allPoints);
          helix1.shouldRefit = true;
          
          // remove second helix
          helices.erase(helices.begin() + iHelix2);
          iHelix2--;
          merged=true;
//          break;
        }
      }
      
//      if(merged) break;
    }
    nPasses++;
  }
}

bool Fitter::MergeHelices2(vector<Helix> &helices)
{
  // Merge helices that are very similar to each other
  bool merged = false;
  
  vector<pair<Helix, vector<int>>> helixLinks;
  
  for(auto helix : helices){
    helixLinks.push_back(make_pair(helix, vector<int>()));
  }
  
  // build links between different helices
  for(int iHelix1=0; iHelix1<helixLinks.size(); iHelix1++){
    Helix &helix1 = helixLinks[iHelix1].first;
    vector<shared_ptr<Point>> points1 = helix1.GetPoints();
    sort(points1.begin(), points1.end());
      
    for(int iHelix2=iHelix1+1; iHelix2<helixLinks.size(); iHelix2++){
      Helix &helix2 = helixLinks[iHelix2].first;
      vector<shared_ptr<Point>> points2 = helix2.GetPoints();
      sort(points2.begin(), points2.end());
        
      vector<shared_ptr<Point>> samePoints;
      set_intersection(points1.begin(), points1.end(),
                       points2.begin(), points2.end(),
                       back_inserter(samePoints));
        
      
      size_t nDifferentPoints = max(points1.size()-samePoints.size(),
                                    points2.size()-samePoints.size());
        
      
      if(nDifferentPoints < 3){
        helixLinks[iHelix1].second.push_back(iHelix2);
        merged = true;
      }
    }
  }
  
  // Merge all linked helices
  vector<int> alreadyUsedHelices;
  vector<int> toRemove;
  
  for(int iHelix=0; iHelix<helixLinks.size(); iHelix++){
    Helix &helix1 = helixLinks[iHelix].first;
    vector<shared_ptr<Point>> points1 = helix1.GetPoints();
    
    
    unordered_set<shared_ptr<Point>> uniquePoints;
    for(auto &p : points1) uniquePoints.insert(p);
    
    alreadyUsedHelices.push_back(iHelix);
    
    for(auto iChild : helixLinks[iHelix].second){
      if(find(alreadyUsedHelices.begin(), alreadyUsedHelices.end(), iChild) != alreadyUsedHelices.end()) continue;
      
      alreadyUsedHelices.push_back(iChild);
      
      Helix &helix2 = helices[iChild];
      vector<shared_ptr<Point>> points2 = helix2.GetPoints();
      for(auto &p : points2) uniquePoints.insert(p);
      toRemove.push_back(iChild);
    }
    vector<shared_ptr<Point>> allPoints(uniquePoints.begin(), uniquePoints.end());
    helix1.ReplacePoints(allPoints);
    helix1.shouldRefit = true;
  }
  
  // Store all helices that were not merged into another helix
  vector<Helix> goodHelices;
  for(int iHelix=0; iHelix<helixLinks.size(); iHelix++){
    if(find(toRemove.begin(), toRemove.end(), iHelix) != toRemove.end()) continue;
    goodHelices.push_back(helixLinks[iHelix].first);
  }
  helices = goodHelices;
  
  return merged;
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
    points[0] = vertex;
    
    double t, distX, distY, distZ;
    double f=0;
    
    // Then add distances to all other points
    for(auto &p : points){
      double x, y, z;
      
      // find helix point for this point's t
      if(helix.GetCharge() < 0) t = atan2(p->GetY() - y0, p->GetX() - x0);
      else                      t = TMath::Pi()/2 - atan2(p->GetX() - x0, p->GetY() - y0);
      x = x0 + (R0 - a*t)*cos(t);
      y = y0 + (R0 - a*t)*sin(t);
      z = z0 + (s0 - b*t)*t;

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
      if(fabs(z-p->GetZ()) > p->GetZerr()){
        double distZ_1 = z - (p->GetZ() + p->GetZerr());
        double distZ_2 = z - (p->GetZ() - p->GetZerr());
        distZ = min(pow(distZ_1, 2), pow(distZ_2, 2));
      }
      
      distX /= fabs(p->GetX());
      distY /= fabs(p->GetY());
      distZ /= fabs(p->GetZ());
      
//      distX /= p->GetXerr() > 0 ? pow(p->GetXerr(), 2) : fabs(p->GetX());
//      distY /= p->GetYerr() > 0 ? pow(p->GetYerr(), 2) : fabs(p->GetY());
//      distZ /= p->GetZerr() > 0 ? pow(p->GetZerr(), 2) : fabs(p->GetZ());
      
      f += distX + distY + distZ;
    }
    return f/(3*helixPoints.size()+8);
  };
  
  ROOT::Fit::Fitter *fitter = new ROOT::Fit::Fitter();
  int nPar = 8;
  ROOT::Math::Functor fitFunction = ROOT::Math::Functor(chi2Function, nPar);
  double pStart[nPar];
  fitter->SetFCN(fitFunction, pStart);
  
  double Lmin = layerRanges[track.GetNtrackerLayers()-1].GetMax();
  double Lmax = layerRanges[track.GetNtrackerLayers()].GetMin();
  
  SetParameter(fitter, 0, "R0", helix.helixParams.R0,
               GetRadiusInMagField(config.minPx, config.minPy, solenoidField),
               GetRadiusInMagField(config.maxPx, config.maxPy, solenoidField));
  SetParameter(fitter, 1, "a" , helix.helixParams.a, 0, 10000);
  SetParameter(fitter, 2, "s0", helix.helixParams.s0, -10000, 10000);
  SetParameter(fitter, 3, "b" , helix.helixParams.b, -10000, 10000);
  SetParameter(fitter, 4, "L" , (helix.GetVertex()->GetX()-10*eventVertex.GetX())/cos(track.GetPhi()),
               Lmin, Lmax);
  SetParameter(fitter, 5, "x0", helix.GetOrigin().GetX(), -layerR[nLayers-1] , layerR[nLayers-1]);
  SetParameter(fitter, 6, "y0", helix.GetOrigin().GetY(), -layerR[nLayers-1] , layerR[nLayers-1]);
  SetParameter(fitter, 7, "z0", helix.GetOrigin().GetZ(), -pixelBarrelZsize  , pixelBarrelZsize);
  
  for(int i=0; i<nPar; i++) pStart[i] = fitter->Config().ParSettings(i).Value();
  
//  if(fitter->FitFCN()) {
  fitter->FitFCN();
  auto result = fitter->Result();
  
  HelixParams resultParams;
  resultParams.R0 = result.GetParams()[0];
  resultParams.a  = result.GetParams()[1];
  resultParams.s0 = result.GetParams()[2];
  resultParams.b  = result.GetParams()[3];
  
  double L  = result.GetParams()[4];
  Point vertex = pointsProcessor.GetPointOnTrack(L, track, eventVertex);
  
  double x0 = result.GetParams()[5];
  double y0 = result.GetParams()[6];
  double z0 = result.GetParams()[7];
  
  Point origin(x0, y0, z0);
  
  helix.helixParams = resultParams;
  helix.SetVertex(vertex);
  helix.UpdateOrigin(origin);
  helix.chi2 = result.MinFcnValue();
  
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
  FixParameter(fitter, 1, "a" , 0);
  SetParameter(fitter, 2, "s0", 0, -10000, 10000);
  FixParameter(fitter, 3, "b" , 0);
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
