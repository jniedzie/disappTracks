//  Fitter.cpp
//
//  Created by Jeremi Niedziela on 14/12/2018.

#include "Fitter.hpp"
#include "Logger.hpp"

Fitter::Fitter(double _maxExecTime) :
eventVertex(Point(0, 0, 0)),
maxExecTime(_maxExecTime)
{
  chi2Function = [&](const double *par) {
    double L  = par[4];
    double x0 = par[5];
    double y0 = par[6];
    double z0 = par[7];

    // Use current L value to calculate decay vertex
    fittingPoints[0] = make_shared<Point>(pointsProcessor.GetPointOnTrack(L, track, eventVertex));
    
    double tMin = pointsProcessor.GetTforPoint(*fittingPoints.front(), Point(x0,y0,z0), charge);
    
    double f=0;
    int nPointFitted=0;
    // Then add distances between this helix and all fitting points
    for(auto &p : fittingPoints){
      if(p->GetSubDetName()=="missing") continue;
      
      double totalDistance = pointsProcessor.GetMinHelixToPointDistance(par, tMin, *p, 0, charge);
      double totalDistEdge = inf;
      
      // In case of endcaps, check if one of the strip edges doesn't give better results
      if(p->IsEndcapHit()){
        double alpha = -atan2(p->GetX(), p->GetY());

        double minImprovement = 0.01;
        Point p_m = *p;
        double prevDist = totalDistance;
        double sizeY = p->GetYerr();
        int i=1;
        do{
          Point p_l(p_m.GetX() + (sizeY * sin(alpha))/pow(2, i),
                    p_m.GetY() + (sizeY * cos(alpha))/pow(2, i),
                    p_m.GetZ());
          
          Point p_r(p_m.GetX() - (sizeY * sin(alpha))/pow(2, i),
                    p_m.GetY() - (sizeY * cos(alpha))/pow(2, i),
                    p_m.GetZ());
        
          double dist_l = pointsProcessor.GetMinHelixToPointDistance(par, tMin, p_l, alpha, charge);
          double dist_r = pointsProcessor.GetMinHelixToPointDistance(par, tMin, p_r, alpha, charge);
          
          if(dist_l < dist_r) p_m = p_l;
          else                p_m = p_r;
          
          if(min(dist_l, dist_r) > totalDistEdge) break;
          
          prevDist = totalDistEdge;
          totalDistEdge = min(dist_l, dist_r);
          i++;
        }
        while( (prevDist-totalDistEdge)/prevDist > minImprovement );
      }
      if(totalDistEdge < totalDistance) totalDistance = totalDistEdge;
      
      f += totalDistance;
      nPointFitted++;
    }
    
    return f/(nPointFitted+nDegreesOfFreedom);
  };
}

Fitter::~Fitter()
{
  
}

Helices Fitter::FitHelices(const Points &_points,
                           const Track &_track,
                           const Point &_eventVertex,
                           int _nTrackerLayers)
{
  startTime = now();
  
  points        = _points;
  pointsByLayer = pointsProcessor.SortByLayer(points);
  pointsByDisk  = pointsProcessor.SortByDisk(points);
  track         = _track;
  eventVertex   = _eventVertex;
  nTrackLayers  = _nTrackerLayers > 0 ? _nTrackerLayers : track.GetNtrackerLayers();
  
  InitLparams();
  charge        = track.GetCharge();

  Helices fittedHelices = PerformFittingCycle();
  if(ShouldStop()) return vector<Helix>();
  
  if(nTrackLayers < config.params["check_opposite_charge_below_Nlayers"]){
    Log(1)<<"Checking opposite charge for default n layers\n";
    charge = -charge;
    Helices fittedHelicesOpposite = PerformFittingCycle();
    if(ShouldStop()) return vector<Helix>();
    fittedHelices.insert(fittedHelices.end(), fittedHelicesOpposite.begin(), fittedHelicesOpposite.end());
    charge = -charge;
  }
  
  if(config.params["allow_one_less_layer"]){
    Log(1)<<"Assuming one less layer\n";
    nTrackLayers-=1;
    InitLparams();
    Helices fittedHelicesOneLess = PerformFittingCycle();
    if(ShouldStop()) return vector<Helix>();
    fittedHelices.insert(fittedHelices.end(), fittedHelicesOneLess.begin(), fittedHelicesOneLess.end());
    
    if(nTrackLayers < config.params["check_opposite_charge_below_Nlayers"]){
      Log(1)<<"Checking opposite charge for one less layer\n";
      charge = -charge;
      Helices fittedHelicesOpposite = PerformFittingCycle();
      if(ShouldStop()) return vector<Helix>();
      fittedHelices.insert(fittedHelices.end(), fittedHelicesOpposite.begin(), fittedHelicesOpposite.end());
      charge = -charge;
    }
    
    nTrackLayers+=1;
  }
  if(config.params["allow_one_more_layer"]){
    Log(1)<<"Assuming one more layer\n";
    nTrackLayers+=1;
    InitLparams();
    Helices fittedHelicesOneMore = PerformFittingCycle();
    if(ShouldStop()) return vector<Helix>();
    fittedHelices.insert(fittedHelices.end(), fittedHelicesOneMore.begin(), fittedHelicesOneMore.end());
    
    if(nTrackLayers < config.params["check_opposite_charge_below_Nlayers"]){
      Log(1)<<"Checking opposite charge for one more layer\n";
      charge = -charge;
      Helices fittedHelicesOpposite = PerformFittingCycle();
      if(ShouldStop()) return vector<Helix>();
      fittedHelices.insert(fittedHelices.end(), fittedHelicesOpposite.begin(), fittedHelicesOpposite.end());
      charge = -charge;
    }
    
    nTrackLayers-=1;
  }
  
  return fittedHelices;
}

Helices Fitter::PerformFittingCycle()
{
  Helices fittedHelices = GetSeeds();
  if(ShouldStop()){
//    cout<<"Execution stopped due to exceeded time"<<endl;
    return vector<Helix>();
  }
  ExtendSeeds(fittedHelices);
  if(ShouldStop()){
//    cout<<"Execution stopped due to exceeded time"<<endl;
    return vector<Helix>();
  }
  if(config.params["merge_final_helices"]) MergeHelices(fittedHelices);
  if(ShouldStop()){
//    cout<<"Execution stopped due to exceeded time"<<endl;
    return vector<Helix>();
  }
  RemoveShortHelices(fittedHelices);
  return fittedHelices;
}

vector<set<int>> Fitter::GetLayersAndDisks()
{
  int lastHitIndex = track.GetNnotEmptyDedxHits()-1;
  bool endcapTrack = track.GetDetTypeForHit(lastHitIndex) == 2; // 2 is endcap
  
  set<int> middleHitLayers;
  set<int> middleHitDisks;
  set<int> lastHitLayers;
  set<int> lastHitDisks;
  
  if(!endcapTrack){
    int lastHitLayer = nTrackLayers-1;
    
    // regardless of the case, for barrel tracks next hit can always be just in the next barrel layer:
    if(lastHitLayer+1 < pointsByLayer.size()) middleHitLayers.insert(lastHitLayer+1);
    if(lastHitLayer+2 < pointsByLayer.size()) lastHitLayers.insert(lastHitLayer+2);
    
    if(lastHitLayer >= 0 && lastHitLayer <= 2){
      middleHitDisks.insert(0);
      lastHitDisks.insert(1);
    }
    if(lastHitLayer >=4 && lastHitLayer <=6){
      middleHitDisks.insert(3);
      lastHitDisks.insert(4);
    }
  }
  else{
    int lastHitDisk = track.GetLayerForHit(lastHitIndex)-1; // iterates from 1
    
    // regardless of the case, for endcap tracks next hit can always be just in the next endcap layer:
    
    if(lastHitDisk+1 < pointsByDisk.size()) middleHitDisks.insert(lastHitDisk+1);
    if(lastHitDisk+2 < pointsByDisk.size()) lastHitDisks.insert(lastHitDisk+2);
    
    if(lastHitDisk >= 0 && lastHitDisk <= 2){
      middleHitLayers.insert(4);
      lastHitLayers.insert(5);
    }
    if(lastHitDisk >= 3 && lastHitDisk <= 5 ){
      middleHitLayers.insert(8);
      lastHitLayers.insert(9);
    }
  }
  vector<set<int>> result = { middleHitLayers, middleHitDisks, lastHitLayers, lastHitDisks };
  return result;
}

pair<Points, Points> Fitter::GetSeedPoints()
{
  Points middlePoints, lastPoints;
  
  auto layersAndDisks = GetLayersAndDisks();
  set<int> middleHitLayers = layersAndDisks[0];
  set<int> middleHitDisks  = layersAndDisks[1];
  set<int> lastHitLayers   = layersAndDisks[2];
  set<int> lastHitDisks    = layersAndDisks[3];
  
  int signZ = sgn(track.GetEta());
  
  for(int middleHitLayer : middleHitLayers){
    middlePoints.insert(middlePoints.end(),
                        pointsByLayer[middleHitLayer].begin(),
                        pointsByLayer[middleHitLayer].end());
  }
  for(int middleHitDisk : middleHitDisks){
    middlePoints.insert(middlePoints.end(),
                        pointsByDisk[GetDisksArrayIndex(middleHitDisk, signZ)].begin(),
                        pointsByDisk[GetDisksArrayIndex(middleHitDisk, signZ)].end());
  }
  
  for(int lastHitLayer : lastHitLayers){
    lastPoints.insert(lastPoints.end(),
                      pointsByLayer[lastHitLayer].begin(),
                      pointsByLayer[lastHitLayer].end());
  }
  for(int lastHitDisk : lastHitDisks){
    lastPoints.insert(lastPoints.end(),
                      pointsByDisk[GetDisksArrayIndex(lastHitDisk, signZ)].begin(),
                      pointsByDisk[GetDisksArrayIndex(lastHitDisk, signZ)].end());
  }
  return make_pair(middlePoints, lastPoints);
}

Helices Fitter::GetSeeds()
{
  Helices seeds;
  
  // Prepare points
  Point trackPointMid = pointsProcessor.GetPointOnTrack(startL, track, eventVertex);
  auto [middlePoints, lastPoints] = GetSeedPoints();
  auto middlePointsRegrouped = pointsProcessor.RegroupNerbyPoints(middlePoints);
  auto lastPointsRegrouped   = pointsProcessor.RegroupNerbyPoints(lastPoints);
  
  int nPairs=0;
  Log(1)<<"Looking for seeds...\n";
  
  for(auto &middlePoints : middlePointsRegrouped){
    if(ShouldStop()) return vector<Helix>();
    
    auto goodMiddlePoints = pointsProcessor.GetGoodMiddleSeedHits(middlePoints,
                                                                  trackPointMid,
                                                                  eventVertex,
                                                                  charge);
    if(goodMiddlePoints.size()==0) continue;
    
    for(auto &lastPoints : lastPointsRegrouped){
      if(ShouldStop()) return vector<Helix>();
      
      auto goodLastPoints = pointsProcessor.GetGoodLastSeedHits(lastPoints,
                                                                goodMiddlePoints,
                                                                trackPointMid,
                                                                charge);
      if(goodLastPoints.size()==0) continue;
      
      nPairs++;
      
      auto helix = FitSeed(goodMiddlePoints, goodLastPoints);
  
      if(!helix) continue;
      if(helix->GetChi2() > config.params["seed_max_chi2"]) continue;
      
      seeds.push_back(*helix);
    }
  }
  
  Log(1)<<"Tested pairs: "<<nPairs<<"\n";
  Log(1)<<"Number of valid seeds: "<<seeds.size()<<"\n";
  
  return seeds;
}

unique_ptr<Helix> Fitter::FitSeed(const Points &middleHits, const Points &lastHits)
{
  auto fitter = GetSeedFitter();
  if(!fitter) return nullptr;
  
  fittingPoints = { make_shared<Point>(0,0,0) };
  fittingPoints.insert(fittingPoints.end(), middleHits.begin(), middleHits.end());
  fittingPoints.insert(fittingPoints.end(), lastHits.begin(), lastHits.end());
  
  fitter->FitFCN();
  auto result = fitter->Result();
  auto resultHelix = make_unique<Helix>(result, track, eventVertex, charge);
  resultHelix->SetLastPoints(middleHits);
  resultHelix->SetLastPoints(lastHits);

  return resultHelix;
}

vector<set<int>> Fitter::GetNextPointLayersAndDisks(const Helix &helix)
{
  auto lastPoints = helix.GetLastPoints();
  int lastPointLayer = lastPoints.front()->GetLayer();
  int lastPointDisk = abs(lastPoints.front()->GetDisk())-1;
  
  set<int> nextPointLayers;
  set<int> nextPointDisks;
  
  if(lastPointLayer >=0 ){
    if(helix.IsIncreasing()){
      if(lastPointLayer+1 < pointsByLayer.size()) nextPointLayers.insert(lastPointLayer+1);
    }
    else{
      if(lastPointLayer-1 >= 0) nextPointLayers.insert(lastPointLayer-1);
    }
  }
  
  if(lastPointDisk+1 < pointsByDisk.size()) nextPointDisks.insert(lastPointDisk+1);
  
  if(lastPointDisk >= 0 && lastPointDisk <= 2)     nextPointLayers.insert(4);
  if(lastPointDisk >= 3 && lastPointDisk <= 5)    nextPointLayers.insert(8);
  
  if(lastPointLayer >= 0 && lastPointLayer <= 2)  nextPointDisks.insert(0);
  if(lastPointLayer >= 3 && lastPointLayer <= 7)  nextPointDisks.insert(3);
  if(lastPointLayer >= 7 && lastPointLayer <= 13) nextPointDisks.insert(6);
  
  vector<set<int>> result = { nextPointLayers, nextPointDisks };
  return result;
}

Points Fitter::GetPossibleNextPoints(const Helix &helix)
{
  auto nextPointLayersAndDisks = GetNextPointLayersAndDisks(helix);
  set<int> nextPointLayers = nextPointLayersAndDisks[0];
  set<int> nextPointDisks  = nextPointLayersAndDisks[1];
  
  int signZ = sgn(track.GetEta());
  
  Points possiblePointsAll, possiblePoints;
  
  for(int nextPointLayer : nextPointLayers){
    possiblePointsAll.insert(possiblePointsAll.end(),
                             pointsByLayer[nextPointLayer].begin(),
                             pointsByLayer[nextPointLayer].end());
  }
  
  for(int nextPointDisk : nextPointDisks){
    possiblePointsAll.insert(possiblePointsAll.end(),
                             pointsByDisk[GetDisksArrayIndex(nextPointDisk, signZ)].begin(),
                             pointsByDisk[GetDisksArrayIndex(nextPointDisk, signZ)].end());
  }
  
  // remove points that are already on helix
  auto helixPoints = helix.GetPoints();
  
  for(auto &pointCandidate : possiblePointsAll){
    if(find_if(helixPoints.begin(), helixPoints.end(),
               [&](const shared_ptr<Point> &p){ return *pointCandidate == *p;}) == helixPoints.end()){
      possiblePoints.push_back(pointCandidate);
    }
  }
  
  return possiblePoints;
}

Points Fitter::GetPossibleTurningBackPoints(const Helix &helix)
{
  auto lastPoints = helix.GetLastPoints();
  int lastPointLayer = lastPoints.front()->GetLayer();
  Points turningBackPointsAll;
  
  
  if(lastPointLayer<pointsByLayer.size()){
    turningBackPointsAll = pointsByLayer[lastPointLayer];
  }
  Points turningBackPoints;
  
  // remove points that are already on helix
  auto helixPoints = helix.GetPoints();
  
  for(auto &pointCandidate : turningBackPointsAll){
    if(find_if(helixPoints.begin(), helixPoints.end(),
               [&](const shared_ptr<Point> &p){ return *pointCandidate == *p;}) == helixPoints.end()){
      turningBackPoints.push_back(pointCandidate);
    }
  }
  return turningBackPoints;
}

Points Fitter::GetGoodNextPoints(const Helix &helix, const Points &points, const set<int> &nextPointLayers)
{
  vector<double> lastPointsT;
  
  for(size_t index : helix.GetLastPointsIndices()){
    lastPointsT.push_back(helix.GetPointT(index));
  }
 
  Points goodPoints;
  auto lastPoints         = helix.GetLastPoints();
  auto secondToLastPoints = helix.GetSecontToLastPoints();
  
  for(shared_ptr<Point> point : points){
    if(helix.GetNlayers() >= config.params["min_layers_for_delta_xy"] &&
       !point->IsEndcapHit()){
      if(!helixProcessor.IsPointCloseToHelixInLayer(helix, *point, *nextPointLayers.begin(), true)){
        Log(3)<<"Δxy too high\n";
        continue;
      }
    }
    else{
      if(!pointsProcessor.IsPhiGood(lastPoints, secondToLastPoints, point, charge)){
        Log(3)<<"Δφ too high\n";
        continue;
      }
    }
    if(!pointsProcessor.IsZgood(lastPoints, point)){
      Log(3)<<"Δz too high\n";
      continue;
    }
    
    double previousT = lastPointsT.front();
    double t = pointsProcessor.GetTforPoint(*point, helix.GetOrigin(), helix.GetCharge());
    while(fabs(t+2*TMath::Pi()-previousT) < fabs(t-previousT)) t += 2*TMath::Pi();
    while(fabs(t-2*TMath::Pi()-previousT) < fabs(t-previousT)) t -= 2*TMath::Pi();
    
    if(!pointsProcessor.IsTgood(lastPointsT, t)) continue;
    goodPoints.push_back(point);
  }
  return goodPoints;
}

Points Fitter::GetGoodTurningBackPoints(const Helix &helix,
                                        const Points &points)
{
  Points goodTurningBackPoints;
  auto lastPoints = helix.GetLastPoints();
  int lastPointLayer = lastPoints.front()->GetLayer();
  
  for(auto &point : points){
    if(!helixProcessor.IsPointCloseToHelixInLayer(helix, *point, lastPointLayer, false)){
      Log(3)<<"Δxy too high\n";
      continue;
    }
    if(!pointsProcessor.IsZgood(lastPoints, point)){
      Log(3)<<"Δz too high\n";
      continue;
    }
    goodTurningBackPoints.push_back(point);
  }
  return goodTurningBackPoints;
}

unique_ptr<Helix> Fitter::TryToExtendHelix(const Helix &helix, const Points &points, bool turnsBack)
{
  if(points.size()==0) return nullptr;
  
  /// Extend helix by the new point and refit its params
  Helix helixCopy(helix);
  helixCopy.SetLastPoints(points);
  if(turnsBack) helixCopy.SetFirstTurningPointIndex((int)helixCopy.GetNpoints()-1);
  RefitHelix(helixCopy);
  
  // check if chi2 is small enough
  if(helixCopy.GetChi2() > config.params["track_max_chi2"]){
    Log(3)<<"Chi2 too high\n";
    return nullptr;
  }
  helixCopy.SetIsPreviousHitMissing(false);
  
  return make_unique<Helix>(helixCopy);
}

void Fitter::ExtendSeeds(Helices &helices)
{
  bool finished;
  int nSteps=0;
  Log(1)<<"Extending seeds...\n";
  
  do{
    if(ShouldStop()) return;
    
    Log(1)<<"Helices before "<<nSteps<<" step: "<<helices.size()<<"\n";
    
    finished = true;
    Helices nextStepHelices;
    
    // for all helices from previous step
    for(Helix &helix : helices){
      if(ShouldStop()) return;
      
      if(helix.GetIsFinished()){
        nextStepHelices.push_back(helix);
      }
      else{
        auto lastPoints    = helix.GetLastPoints();
        int lastPointLayer = lastPoints.front()->GetLayer();
        int lastPointDisk  = abs(lastPoints.front()->GetDisk())-1;
        if(lastPointLayer < 0 && lastPointDisk < 0) continue;
        
        auto nextPointLayersAndDisks = GetNextPointLayersAndDisks(helix);
        set<int> nextPointLayers = nextPointLayersAndDisks[0];
        
        // fist, check if helix crosses next layer
        Point pA, pB;
        bool crossesNextLayer = helixProcessor.GetIntersectionWithLayer(helix, *nextPointLayers.begin(), pA, pB);
        
        bool lastPointIsEndcap = lastPoints.front()->IsEndcapHit();
        bool turnsBack = !(crossesNextLayer || lastPointIsEndcap);
        
        Points possiblePoints = turnsBack ? GetPossibleTurningBackPoints(helix) : GetPossibleNextPoints(helix);
        vector<Points> possiblePointsRegrouped = pointsProcessor.RegroupNerbyPoints(possiblePoints);
        Helices extendedHelices;
        
        // Find points that could extend this helix
        for(Points &points : possiblePointsRegrouped){
          if(ShouldStop()) return;
          
          Points goodPoints = turnsBack ? GetGoodTurningBackPoints(helix, points) : GetGoodNextPoints(helix, points, nextPointLayers);
        
          unique_ptr<Helix> extendedHelix = TryToExtendHelix(helix, goodPoints, turnsBack);
          if(extendedHelix) extendedHelices.push_back(*extendedHelix);
        }

        // if it was possible to extend the helix
        if(extendedHelices.size() != 0){
          nextStepHelices.insert(nextStepHelices.end(),
                                 extendedHelices.begin(),
                                 extendedHelices.end());
          finished = false;
        }
        else{ // if helix could not be extended
          
          // if not missing hits are allowed
          if(helix.GetNmissingHits()       >= config.params["max_n_missing_hits"] ||
             helix.GetNmissingHitsInRow()  >= config.params["max_n_missing_hits_in_raw"] ||
             
             // TODO: For the moment no missing hits if last one was in endcaps!!
             lastPoints.front()->IsEndcapHit()){
            
            // if last hit was a missing hit, remove it
            if(helix.GetLastPoints().back()->GetSubDetName()=="missing"){
              helix.RemoveLastPoints();
              RefitHelix(helix);
            }
            
            helix.SetIsFinished(true);
            nextStepHelices.push_back(helix);
          }
          // add a missing hit
          else{
            int missingHitLayer = helix.GetFirstTurningPointIndex() > 0 ? lastPointLayer-1 : lastPointLayer+1;
            if(!crossesNextLayer) missingHitLayer = lastPointLayer;
            shared_ptr<Point> missingHit = helixProcessor.GetPointCloseToHelixInLayer(helix, missingHitLayer);

            helix.SetLastPoints({missingHit});
            missingHit = helix.GetLastPoints().back();
            missingHit->SetZ(lastPoints.front()->GetZ());
            missingHit->SetSubDetName("missing");
            helix.IncreaseMissingHits();
            
            nextStepHelices.push_back(helix);
            finished = false;
          }
        }
      }
    }
    
    helices.clear();
    for(auto &h : nextStepHelices) helices.push_back(h);
    nSteps++;
  }
  while(!finished);
  Log(1)<<"Candidates found: "<<helices.size()<<"\n";
}

void Fitter::RefitHelix(Helix &helix)
{
  helix.SetChi2(inf);
  fittingPoints = helix.GetPoints();
  charge = helix.GetCharge();
  
  ROOT::Fit::Fitter *fitter = new ROOT::Fit::Fitter();
  int nPar = 8;
  ROOT::Math::Functor fitFunction = ROOT::Math::Functor(chi2Function, nPar);
  double pStart[nPar];
  fitter->SetFCN(fitFunction, pStart);
  
  double startR0     = helix.GetRadius(helix.GetTmin());
  double startRslope = helix.GetRadiusFactor();
  double startS0     = helix.GetSlope(helix.GetTmin());
  double startSslope = helix.GetSlopeFactor();
  double startX0     = helix.GetOrigin().GetX();
  double startY0     = helix.GetOrigin().GetY();
  double startZ0     = helix.GetOrigin().GetZ();
  
  if(startR0      < config.params["min_R0"]     || startR0     > config.params["max_R0"]      ||
     startRslope  < config.params["min_Rslope"]  || startRslope > config.params["max_Rslope"] ||
     startS0      < config.params["min_S0"]      || startS0     > config.params["max_S0"]     ||
     startSslope  < config.params["min_Sslope"]  || startSslope > config.params["max_Sslope"] ||
     startX0      < config.params["min_X0"]      || startX0     > config.params["max_X0"]     ||
     startY0      < config.params["min_Y0"]      || startY0     > config.params["max_Y0"]     ||
     startZ0      < config.params["min_Z0"]      || startZ0     > config.params["max_Z0"]){
    Log(1)<<"ERROR -- wrong params in RefitHelix, which should never happen...\n";
    return;
  }
  
  SetParameter(fitter, 0, "R0", startR0     , config.params["min_R0"]     , config.params["max_R0"]);
  SetParameter(fitter, 1, "a" , startRslope , config.params["min_Rslope"] , config.params["max_Rslope"]);
  SetParameter(fitter, 2, "s0", startS0     , config.params["min_S0"]     , config.params["max_S0"]);
  SetParameter(fitter, 3, "b" , startSslope , config.params["min_Sslope"] , config.params["max_Sslope"]);
  SetParameter(fitter, 4, "L" , startL      , minL      , maxL);
  SetParameter(fitter, 5, "x0", startX0     , config.params["min_X0"]     , config.params["max_X0"]);
  SetParameter(fitter, 6, "y0", startY0     , config.params["min_Y0"]     , config.params["max_Y0"]);
  SetParameter(fitter, 7, "z0", startZ0     , config.params["min_Z0"]     , config.params["max_Z0"]);

  for(int i=0; i<nPar; i++) pStart[i] = fitter->Config().ParSettings(i).Value();
  
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
  
  helix.SetParams(resultParams);
  helix.SetVertex(vertex);
  helix.UpdateOrigin(origin);
  helix.SetChi2(result.MinFcnValue());
}

void Fitter::MergeHelices(Helices &helices)
{
  Log(1)<<"Merging overlapping helices...";
  
  Helices helicesToMerge;
  Helices tooShortToMerge;
  // Remove very short candidates which should not even be merged with others
  for(auto helix : helices){
    if(ShouldStop()) return;
    
    if(helix.GetNpoints() >= config.params["candidate_min_n_points"]) helicesToMerge.push_back(helix);
    else                                                              tooShortToMerge.push_back(helix);
  }
  
  // Merge similar candidates
  while(LinkAndMergeHelices(helicesToMerge));
  
  // Insert back those that were too short to merge
  helicesToMerge.insert(helicesToMerge.end(), tooShortToMerge.begin(), tooShortToMerge.end());
  
  Log(1)<<" merged down to: "<<helicesToMerge.size()<<"\n";
  helices = helicesToMerge;
}

ROOT::Fit::Fitter* Fitter::GetSeedFitter()
{
  auto fitter = new ROOT::Fit::Fitter();
  int nPar = 8;
  auto fitFunction = ROOT::Math::Functor(chi2Function, nPar);
  double pStart[nPar];
  fitter->SetFCN(fitFunction, pStart);
  
  // Calculate initial parameters as good as we can at this point.
  Point trackPoint = pointsProcessor.GetPointOnTrack(startL, track, eventVertex);
  
  double startX0, minX0, maxX0, startY0, minY0, maxY0;
  GetXYranges(trackPoint, startX0, minX0, maxX0, startY0, minY0, maxY0);
  
  if(startX0 < minX0 || startX0 > maxX0){
    Log(1)<<"ERROR -- x0:"<<startX0<<"\tmin:"<<minX0<<"\tmax:"<<maxX0<<"\n";
    if(config.params["require_good_starting_values"]) return nullptr;
  }
  if(startY0 < minY0 || startY0 > maxY0){
    Log(1)<<"ERROR -- y0:"<<startY0<<"\tmin:"<<minY0<<"\tmax:"<<maxY0<<"\n";
    if(config.params["require_good_starting_values"]) return nullptr;
  }
  
  // -- calculate slope from the track momentum direction (pion usually follows this direction)
  double startS0 = config.params["start_R0"] * trackPoint.GetVectorSlopeC();
  
  if(startS0 < config.params["min_S0"] || startS0 > config.params["max_S0"]){
    Log(1)<<"ERROR -- S0:"<<startS0<<"\tmin:"<<config.params["min_S0"]<<"\tmax:"<<config.params["max_S0"]<<"\n";
    if(config.params["require_good_starting_values"]) return nullptr;
  }
  
  // -- get t param of the track point and calculate Z position of the vertex
  Point origin(startX0, startY0, 0);
  double tTrack = pointsProcessor.GetTforPoint(trackPoint, origin, charge);
  double startZ0 = -charge * (trackPoint.GetZ() - startS0 * tTrack);

  if(startZ0 < config.params["min_Z0"] || startZ0 > config.params["max_Z0"]){
    Log(1)<<"ERROR -- z0:"<<startZ0<<"\tmin:"<<config.params["min_Z0"]<<"\tmax:"<<config.params["max_Z0"]<<"\n";
    if(config.params["require_good_starting_values"]) return nullptr;
  }
  
  // Set calculated initial param values
  SetParameter(fitter, 0, "R0", config.params["start_R0"]  , config.params["min_R0"]   , config.params["max_R0"] ); // limits from MC
  SetParameter(fitter, 2, "s0", startS0 , config.params["min_S0"]   , config.params["max_S0"] );
  SetParameter(fitter, 4, "L" , startL  , minL    , maxL  );
  SetParameter(fitter, 5, "x0", startX0 , config.params["min_X0"]   , config.params["max_X0"] );
  SetParameter(fitter, 6, "y0", startY0 , config.params["min_Y0"]   , config.params["max_Y0"] );
  SetParameter(fitter, 7, "z0", startZ0 , config.params["min_Z0"]   , config.params["max_Z0"] );
  
  // With 3 points we don't know how fast will radius and slope decrease:
  FixParameter(fitter, 1, "a" ,  0.0000001);
  FixParameter(fitter, 3, "b" , -0.0000001);
  
//  for(int i=0; i<nPar; i++) pStart[i] = fitter->Config().ParSettings(i).Value();
//  fitter->SetFCN(fitFunction, pStart);
  
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

void Fitter::GetXYranges(const Point &trackPoint,
                         double &startX0, double &minX0, double &maxX0,
                         double &startY0, double &minY0, double &maxY0)
{
  double x = trackPoint.GetX();
  double y = trackPoint.GetY();
  
  minX0 = maxX0 = x;
  minY0 = maxY0 = y;
  
  startX0 = x + sgn(y)*charge*config.params["start_R0"]/sqrt(pow(x/y, 2)+1);
  startY0 = y - sgn(x)*charge*config.params["start_R0"]/sqrt(pow(y/x, 2)+1);
  
  if(charge*y < 0)  minX0 -= config.params["max_R0"]/sqrt(pow(x/y, 2)+1);
  else              maxX0 += config.params["max_R0"]/sqrt(pow(x/y, 2)+1);
  
  if(charge*x > 0)  minY0 -= config.params["max_R0"]/sqrt(pow(y/x, 2)+1);
  else              maxY0 += config.params["max_R0"]/sqrt(pow(y/x, 2)+1);
}

void Fitter::RemoveShortHelices(Helices &helices)
{
  Log(1)<<"Removing very short merged helices...";
  Helices longHelices;
  for(auto helix : helices){
    if(ShouldStop()) return;
    
    if(helix.GetNpoints() >= config.params["track_min_n_points"] &&
       helix.GetNlayers() >= config.params["track_min_n_layers"]){
      helix.SetNrecLayers((int)helix.GetNlayers());
      helix.SetNrecHits((int)helix.GetNpoints());
      longHelices.push_back(helix);
    }
  }
  Log(1)<<" long merged helices: "<<longHelices.size()<<"\n";
  helices = longHelices;
}

bool Fitter::LinkAndMergeHelices(Helices &helices)
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
    Points points1 = helix1.GetPoints();
    
    sort(points1.begin(), points1.end());
    
    for(int iHelix2=iHelix1+1; iHelix2<helixLinks.size(); iHelix2++){
      Helix &helix2 = helixLinks[iHelix2].first;
      Points points2 = helix2.GetPoints();
      sort(points2.begin(), points2.end());
      
      Points samePoints;
      set_intersection(points1.begin(), points1.end(),
                       points2.begin(), points2.end(),
                       back_inserter(samePoints));

      size_t nDifferentPoints = max(points1.size()-samePoints.size(),
                                    points2.size()-samePoints.size());
      
      
      if(nDifferentPoints <= config.params["merging_max_different_point"]){
        helixLinks[iHelix1].second.push_back(iHelix2);
        merged = true;
      }
    }
  }
  
  // Merge all linked helices
  vector<int> alreadyUsedHelices;
  vector<int> toRemove;
  
  for(int iHelix=0; iHelix<helixLinks.size(); iHelix++){
    if(helixLinks[iHelix].second.size() == 0) continue; // if there's nothing to link, skip it
    
    Helix &helix1 = helixLinks[iHelix].first;
    Helix helix1Copy(helix1);
    
    Points allPoints = helix1.GetPoints();
    vector<double> allPointsT = helix1.GetPointsT();

    alreadyUsedHelices.push_back(iHelix);
    
    vector<int> childIndices;
    
    for(auto iChild : helixLinks[iHelix].second){
      if(find(alreadyUsedHelices.begin(), alreadyUsedHelices.end(), iChild) != alreadyUsedHelices.end()) continue;
      
      alreadyUsedHelices.push_back(iChild);
      
      Helix &helix2 = helices[iChild];
      Points points2 = helix2.GetPoints();
      vector<double> pointsT2 = helix2.GetPointsT();
      
      for(size_t index=0; index<points2.size(); index++){
        auto p = points2[index];
        if(p->GetLayer()<0) continue; // don't merge in the decay vertex
        if(find_if(allPoints.begin(), allPoints.end(), [&](const shared_ptr<Point> &p1){return *p == *p1;}) != allPoints.end()) continue; // don't add the same point twice
        
        allPoints.push_back(p);
        allPointsT.push_back(pointsT2[index]);
      }
      
      toRemove.push_back(iChild);
      childIndices.push_back(iChild);
    }
    
    helix1Copy.SetPointsAndSortByT(allPoints, allPointsT);
    RefitHelix(helix1Copy);
    
    if(helix1Copy.GetChi2() < config.params["track_max_chi2"]){
      helix1 = helix1Copy;
    }
    else{
      // if merged helix had very bad chi2, cancel the whole merging operation
      for(int i=0; i<toRemove.size();){
        if(find(childIndices.begin(), childIndices.end(), i) != childIndices.end()) toRemove.erase(toRemove.begin()+i);
        else i++;
      }
      for(int i=0; i<alreadyUsedHelices.size();){
        if(find(childIndices.begin(), childIndices.end(), i) != childIndices.end())
          alreadyUsedHelices.erase(alreadyUsedHelices.begin()+i);
        else i++;
      }
    }
  }
  
  // Store all helices that were not merged into another helix
  Helices goodHelices;
  for(int iHelix=0; iHelix<helixLinks.size(); iHelix++){
    if(find(toRemove.begin(), toRemove.end(), iHelix) != toRemove.end()) continue;
    goodHelices.push_back(helixLinks[iHelix].first);
  }
  helices = goodHelices;
  
  return merged;
}

void Fitter::InitLparams()
{
  int lastHitIndex = track.GetNnotEmptyDedxHits()-1;
  bool endcapTrack = track.GetDetTypeForHit(lastHitIndex) == 2; // 2 is endcap
  
  if(!endcapTrack){
    minL   = layerRanges[nTrackLayers-1].GetMin();
    maxL   = layerRanges[nTrackLayers].GetMax();
    
  }
  else{
    size_t diskIndex = track.GetLayerForHit(lastHitIndex);
    
    minL = fabs(diskRanges[diskIndex-1].front().GetMin() * tan(track.GetTheta()));
    maxL = fabs(diskRanges[diskIndex].back().GetMax() * tan(track.GetTheta()));
  }
  
  startL = (minL+maxL)/2.;
}

bool Fitter::ShouldStop()
{
  return maxExecTime > 0 ? duration(startTime, now()) > maxExecTime : false;
}
