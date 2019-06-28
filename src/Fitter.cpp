//  Fitter.cpp
//
//  Created by Jeremi Niedziela on 14/12/2018.

#include "Fitter.hpp"

Fitter::Fitter() :
eventVertex(Point(0, 0, 0))
{
  chi2Function = [&](const double *par) {
    double R0 = par[0];
    double a  = par[1];
    double s0 = par[2];
    double b  = par[3];
    double L  = par[4];
    double x0 = par[5];
    double y0 = par[6];
    double z0 = par[7];

    // First add distance to the vertex
    Point origin(x0, y0, z0);
    fittingPoints[0] = make_shared<Point>(pointsProcessor.GetPointOnTrack(L, track, eventVertex));
    
    double tMin = pointsProcessor.GetTforPoint(*fittingPoints.front(), origin, charge);
    
    double t, distX, distY, distZ, x, y, z;
    double f=0;
    
    // Then add distances to all other points
    for(auto &p : fittingPoints){
      if(p->GetSubDetName()=="missing") continue;
      
      t = pointsProcessor.GetTforPoint(*p, origin, charge);
      
      x = x0 + GetRofT(R0, a, tMin, t, charge)*cos(t);
      y = y0 + GetRofT(R0, a, tMin, t, charge)*sin(t);
      z = charge*z0 + GetSofT(s0, b, tMin, t, charge)*t;
      
      // calculate distance between helix and point's boundary (taking into account its errors)
      distX = distY = distZ = 0;
      
      if(fabs(x-p->GetX()) > p->GetXerr()){
        double distX_1 = x - (p->GetX() + p->GetXerr());
        double distX_2 = x - (p->GetX() - p->GetXerr());
        distX = min(pow(distX_1, 2), pow(distX_2, 2));
      }
      if(fabs(y-p->GetY()) > p->GetYerr()){
        double distY_1 = y - (p->GetY() + p->GetYerr());
        double distY_2 = y - (p->GetY() - p->GetYerr());
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
      
      f += distX + distY + distZ;
    }
    return f/(3*fittingPoints.size()+nDegreesOfFreedom);
  };
}

Fitter::~Fitter()
{
  
}

vector<Helix> Fitter::FitHelices(const vector<shared_ptr<Point>> &_points,
                                 const Track &_track,
                                 const Point &_eventVertex)
{
  points       = _points;
  track        = _track;
  eventVertex  = _eventVertex;
  nTrackLayers = track.GetNtrackerLayers();
  InitLparams();
  charge = track.GetCharge();

  vector<Helix> fittedHelices = PerformFittingCycle();
  
  if(nTrackLayers < config.checkOppositeChargeBelowNlayers){
    if(config.verbosity>0) cout<<"Checking opposite charge for default n layers"<<endl;
    charge = -charge;
    vector<Helix> fittedHelicesOpposite = PerformFittingCycle();
    fittedHelices.insert(fittedHelices.end(), fittedHelicesOpposite.begin(), fittedHelicesOpposite.end());
    charge = -charge;
  }
  
  if(config.allowOneLessLayer){
    if(config.verbosity>0) cout<<"Assuming one less layer"<<endl;
    nTrackLayers-=1;
    InitLparams();
    vector<Helix> fittedHelicesOneLess = PerformFittingCycle();
    fittedHelices.insert(fittedHelices.end(), fittedHelicesOneLess.begin(), fittedHelicesOneLess.end());
    
    if(nTrackLayers < config.checkOppositeChargeBelowNlayers){
      if(config.verbosity>0) cout<<"Checking opposite charge for one less layer"<<endl;
      charge = -charge;
      vector<Helix> fittedHelicesOpposite = PerformFittingCycle();
      fittedHelices.insert(fittedHelices.end(), fittedHelicesOpposite.begin(), fittedHelicesOpposite.end());
      charge = -charge;
    }
    
    nTrackLayers+=1;
  }
  if(config.allowOneMoreLayer){
    if(config.verbosity>0) cout<<"Assuming one more layer"<<endl;
    nTrackLayers+=1;
    InitLparams();
    vector<Helix> fittedHelicesOneMore = PerformFittingCycle();
    fittedHelices.insert(fittedHelices.end(), fittedHelicesOneMore.begin(), fittedHelicesOneMore.end());
    
    if(nTrackLayers < config.checkOppositeChargeBelowNlayers){
      if(config.verbosity>0) cout<<"Checking opposite charge for one more layer"<<endl;
      charge = -charge;
      vector<Helix> fittedHelicesOpposite = PerformFittingCycle();
      fittedHelices.insert(fittedHelices.end(), fittedHelicesOpposite.begin(), fittedHelicesOpposite.end());
      charge = -charge;
    }
    
    nTrackLayers-=1;
  }
  
  return fittedHelices;
}

vector<Helix> Fitter::PerformFittingCycle()
{
  vector<vector<shared_ptr<Point>>> pointsByLayer = pointsProcessor.SortByLayer(points);
  
  vector<Helix> fittedHelices = GetSeeds(pointsByLayer);
  ExtendSeeds(fittedHelices, pointsByLayer);
  if(config.mergeFinalHelices) MergeHelices(fittedHelices);
  RemoveShortHelices(fittedHelices);
  
  if(config.verbosity>0) cout<<"Refitting surviving helices...";
  for(auto &helix : fittedHelices){
    if(helix.GetShouldRefit()) RefitHelix(helix);
  }
  return fittedHelices;
}

vector<Helix> Fitter::GetSeeds(vector<vector<shared_ptr<Point>>> pointsByLayer)
{
  // find possible middle and last seeds' points
  auto middlePointsRegrouped = pointsProcessor.RegroupNerbyPoints(pointsByLayer[nTrackLayers]);
  auto lastPointsRegrouped   = pointsProcessor.RegroupNerbyPoints(pointsByLayer[nTrackLayers+1]);
  
  double lmin = layerR[nTrackLayers-1];
  double lmax = layerR[nTrackLayers];
  //  double lmin = layerRanges[nTrackLayers-1].GetMin();
  //  double lmax = layerRanges[nTrackLayers].GetMax();
  double lstart = (lmin+lmax)/2;
  
  Point trackPointMid = pointsProcessor.GetPointOnTrack(minL, track, eventVertex);
  
  vector<Helix> seeds;
  int nPairs=0;
  if(config.verbosity>0) cout<<"Looking for seeds..."<<endl;
  for(auto &middlePoints : middlePointsRegrouped){
    auto goodMiddlePoints = pointsProcessor.GetGoodMiddleSeedHits(middlePoints,
                                                                  trackPointMid,
                                                                  eventVertex,
                                                                  charge);
    if(goodMiddlePoints.size()==0) continue;
    
    for(auto &lastPoints : lastPointsRegrouped){
      
      auto goodLastPoints = pointsProcessor.GetGoodLastSeedHits(lastPoints,
                                                                goodMiddlePoints,
                                                                trackPointMid,
                                                                charge);
      if(goodLastPoints.size()==0) continue;
      
      nPairs++;
      vector<shared_ptr<Point>> points;
      points.insert(points.end(), goodMiddlePoints.begin(), goodMiddlePoints.end());
      points.insert(points.end(), goodLastPoints.begin()  , goodLastPoints.end());

      auto helix = FitSeed(points);
  
      if(!helix) continue;
      if(helix->GetChi2() > config.seedMaxChi2) continue;
      
      seeds.push_back(*helix);
    }
  }
  if(config.verbosity>0){
    cout<<"Tested pairs: "<<nPairs<<endl;
    cout<<"Number of valid seeds: "<<seeds.size()<<endl;
  }
  return seeds;
}

unique_ptr<Helix> Fitter::FitSeed(const vector<shared_ptr<Point>> &seedPoints)
{
  auto fitter = GetSeedFitter(seedPoints);
  if(!fitter) return nullptr;
  
  fittingPoints = { make_shared<Point>(0,0,0) };
  fittingPoints.insert(fittingPoints.end(), seedPoints.begin(), seedPoints.end());
  
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
  
  resultHelix = make_unique<Helix>(resultParams, vertex, origin, seedPoints, track);
  resultHelix->SetCharge(charge);
  resultHelix->UpdateOrigin(origin);
  resultHelix->SetChi2(result.MinFcnValue());
  
  return resultHelix;
}

void Fitter::ExtendSeeds(vector<Helix> &helices,
                         const vector<vector<shared_ptr<Point>>> &pointsByLayer)
{
  bool finished;
  int nSteps=0;
  if(config.verbosity>0) cout<<"Extending seeds..."<<endl;
  do{
    if(config.verbosity>0) cout<<"Helices before "<<nSteps<<" step: "<<helices.size()<<endl;
    
    finished = true;
    vector<Helix> nextStepHelices;
    
    // for all helices from previous step
    for(Helix &helix : helices){
      
      if(helix.GetIsFinished()){
        nextStepHelices.push_back(helix);
      }
      else{
        vector<Helix> extendedHelices;
        
        auto helixPoints        = helix.GetPoints();
        auto lastPoints         = helix.GetLastPoints();
        auto secondToLastPoints = helix.GetSecontToLastPoints();
        
        shared_ptr<Point> turningPoint = nullptr;
        if(helix.GetFirstTurningPointIndex() > 0) turningPoint = helixPoints[helix.GetFirstTurningPointIndex()];
        
        int lastPointLayer = lastPoints.front()->GetLayer();
        if(lastPointLayer < 0) continue;
        int nextPointLayer = helix.IsIncreasing() ? lastPointLayer+1 : lastPointLayer-1;
        
        // fist, check if helix crosses next layer
        Point pA, pB;
        bool crossesNextLayer = helixProcessor.GetIntersectionWithLayer(helix, nextPointLayer, pA, pB);
        
        // try to extend to next layer
        if(crossesNextLayer || !config.allowTurningBack){
  
          // Find points that could extend this helix
          vector<shared_ptr<Point>> possiblePointsAll;
          if(lastPointLayer+1 < pointsByLayer.size() && lastPointLayer-1 >= 0){
            possiblePointsAll = pointsByLayer[nextPointLayer];
          }
          vector<shared_ptr<Point>> possiblePoints;
          
          // remove points that are already on helix
          for(auto &pointCandidate : possiblePointsAll){
            if(find_if(helixPoints.begin(), helixPoints.end(), [&](const shared_ptr<Point> &p){ return *pointCandidate == *p;}) == helixPoints.end()){
              possiblePoints.push_back(pointCandidate);
            }
          }
          
          auto possiblePointsRegrouped = pointsProcessor.RegroupNerbyPoints(possiblePoints);
          
          for(auto &points : possiblePointsRegrouped){
            
            vector<shared_ptr<Point>> goodPoints;
            
            vector<double> lastPointsT;
            
            for(size_t index : helix.GetLastPointsIndices()){
              lastPointsT.push_back(helix.GetPointT(index));
            }
            
            for(auto &point : points){
              
              if(helix.GetNlayers() >= config.minLayersForDeltaXY){
                if(!helixProcessor.IsPointCloseToHelixInLayer(helix, *point, nextPointLayer, true)) continue;
              }
              else{
                if(!pointsProcessor.IsPhiGood(lastPoints, secondToLastPoints, point, charge)) continue;
              }
              if(!pointsProcessor.IsZgood(lastPoints, point)) continue;
              
              double previousT = lastPointsT.front();
              double t = pointsProcessor.GetTforPoint(*point, helix.GetOrigin(), helix.GetCharge());
              while(fabs(t+2*TMath::Pi()-previousT) < fabs(t-previousT)) t += 2*TMath::Pi();
              while(fabs(t-2*TMath::Pi()-previousT) < fabs(t-previousT)) t -= 2*TMath::Pi();
              
              if(!pointsProcessor.IsTgood(lastPointsT, t)) continue;
              goodPoints.push_back(point);
            }
            
            if(goodPoints.size()==0) continue;
            
            /// Extend helix by the new point and refit its params
            Helix helixCopy(helix);
            for(auto &point : goodPoints) helixCopy.AddPoint(point);
            RefitHelix(helixCopy);
            
            // check if chi2 is small enough
            if(helixCopy.GetChi2() > config.trackMaxChi2) continue;
            
            // if we reached this point, it means that this hit is not missing
            helixCopy.SetIsPreviousHitMissing(false);
            extendedHelices.push_back(helixCopy);
          }
        }
        else{
          // try to extend to the same layer (turning back)
          
          vector<shared_ptr<Point>> turningBackPointsAll;
          
          if(lastPointLayer<pointsByLayer.size()){
            turningBackPointsAll = pointsByLayer[lastPointLayer];
          }
          vector<shared_ptr<Point>> turningBackPoints;
          
          // remove points that are already on helix
          for(auto &pointCandidate : turningBackPointsAll){
            if(find_if(helixPoints.begin(), helixPoints.end(), [&](const shared_ptr<Point> &p){ return *pointCandidate == *p;}) == helixPoints.end()){
              turningBackPoints.push_back(pointCandidate);
            }
          }
          auto turningBackPointsRegrouped = pointsProcessor.RegroupNerbyPoints(turningBackPoints);
          
          for(auto &turningBackPoints : turningBackPointsRegrouped){

            vector<shared_ptr<Point>> goodTurningBackPoints;
            
            for(auto &point : turningBackPoints){
              if(!helixProcessor.IsPointCloseToHelixInLayer(helix, *point, lastPointLayer, false)) continue;
              if(!pointsProcessor.IsZgood(lastPoints, point)) continue;
              goodTurningBackPoints.push_back(point);
            }
            if(goodTurningBackPoints.size()==0) continue;
            
            /// Extend helix by the new point and refit its params
            Helix helixCopy(helix);
            for(auto &point : goodTurningBackPoints) helixCopy.AddPoint(point);
            helixCopy.SetFirstTurningPointIndex(helixCopy.GetNpoints()-1);
            RefitHelix(helixCopy);
            
            // check if chi2 is small enough
            if(helixCopy.GetChi2() > config.trackMaxChi2) continue;
            
            // if we reached this point, it means that this hit is not missing
            helixCopy.SetIsPreviousHitMissing(false);
            
            extendedHelices.push_back(helixCopy);
          }
        
          if(config.mergeAtTurnBack) while(LinkAndMergeHelices(extendedHelices));
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
          if(helix.GetNmissingHits()       >= config.maxNmissingHits ||
             helix.GetNmissingHitsInRow()  >= config.maxNmissingHitsInRow){
            
            // if last hit was a missing hit, remove it
            if(helix.GetLastPoints().back()->GetSubDetName()=="missing"){
              helix.RemoveLastPoint();
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

            helix.AddPoint(missingHit); // this will set missing hit's T
            missingHit = helix.GetLastPoints().back();
//            double t = missingHit->GetT();
//            missingHit->SetZ(-helix.GetCharge()*helix.GetOrigin().GetZ() + helix.GetSlope(t)*t + 10*eventVertex.GetZ());
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
  if(config.verbosity>0) cout<<"Candidates found: "<<helices.size()<<endl;
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
  
  if(startR0      < minR0      || startR0     > maxR0     ||
     startRslope  < minRslope  || startRslope > maxRslope ||
     startS0      < minS0      || startS0     > maxS0     ||
     startSslope  < minSslope  || startSslope > maxSslope ||
     startX0      < minX0      || startX0     > maxX0     ||
     startY0      < minY0      || startY0     > maxY0     ||
     startZ0      < minZ0      || startZ0     > maxZ0){
    if(config.verbosity>0) cout<<"ERROR -- wrong params in RefitHelix, which should never happen..."<<endl;
    return;
  }
  
  SetParameter(fitter, 0, "R0", startR0     , minR0     , maxR0);
  SetParameter(fitter, 1, "a" , startRslope , minRslope , maxRslope);
  SetParameter(fitter, 2, "s0", startS0     , minS0     , maxS0);
  SetParameter(fitter, 3, "b" , startSslope , minSslope , maxSslope);
  SetParameter(fitter, 4, "L" , startL      , minL      , maxL);
  SetParameter(fitter, 5, "x0", startX0     , minX0     , maxX0);
  SetParameter(fitter, 6, "y0", startY0     , minY0     , maxY0);
  SetParameter(fitter, 7, "z0", startZ0     , minZ0     , maxZ0);

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

void Fitter::MergeHelices(vector<Helix> &helices)
{
  if(config.verbosity>0) cout<<"Merging overlapping helices...";
  
  vector<Helix> helicesToMerge;
  vector<Helix> tooShortToMerge;
  // Remove very short candidates which should not even be merged with others
  for(auto helix : helices){
    if(helix.GetNpoints() >= config.candidateMinNpoints) helicesToMerge.push_back(helix);
    else                                                 tooShortToMerge.push_back(helix);
  }
  
  // Merge similar candidates
  while(LinkAndMergeHelices(helicesToMerge));
  
  // Insert back those that were too short to merge
  helicesToMerge.insert(helicesToMerge.end(), tooShortToMerge.begin(), tooShortToMerge.end());
  
  if(config.verbosity>0) cout<<" merged down to: "<<helicesToMerge.size()<<endl;
  helices = helicesToMerge;
}

ROOT::Fit::Fitter* Fitter::GetSeedFitter(const vector<shared_ptr<Point>> &points)
{
  auto fitter = new ROOT::Fit::Fitter();
  
  // This is a stupid hack to be able to set fit parameters before actually setting a fitting function
  // Params we want to set only once, but function will change in each iteration of the loop, because
  // it captures loop iterators.
  auto f = [&](const double*) {return 0;};
  int nPar = 8;
  ROOT::Math::Functor fitFunction = ROOT::Math::Functor(f, nPar);
  double pStart[nPar];
  fitter->SetFCN(fitFunction, pStart);
  
  // Calculate initial parameters as good as we can at this point.
  Point trackPoint = pointsProcessor.GetPointOnTrack(startL, track, eventVertex);
  
  double startX0, minX0, maxX0, startY0, minY0, maxY0;
  GetXYranges(trackPoint, startX0, minX0, maxX0, startY0, minY0, maxY0);
  
  if(startX0 < minX0 || startX0 > maxX0){
    if(config.verbosity>0) cout<<"ERROR -- x0:"<<startX0<<"\tmin:"<<minX0<<"\tmax:"<<maxX0<<endl;
    if(config.requireGoodStartingValues) return nullptr;
  }
  if(startY0 < minY0 || startY0 > maxY0){
    if(config.verbosity>0) cout<<"ERROR -- y0:"<<startY0<<"\tmin:"<<minY0<<"\tmax:"<<maxY0<<endl;
    if(config.requireGoodStartingValues) return nullptr;
  }
  
  // -- calculate slope from the track momentum direction (pion usually follows this direction)
  double startS0 = startR * trackPoint.GetVectorSlopeC();
  
  if(startS0 < minS0 || startS0 > maxS0){
    if(config.verbosity>0) cout<<"ERROR -- S0:"<<startS0<<"\tmin:"<<minS0<<"\tmax:"<<maxS0<<endl;
    if(config.requireGoodStartingValues) return nullptr;
  }
  
  // -- get t param of the track point and calculate Z position of the vertex
  Point origin(startX0, startY0, 0);
  double tTrack = pointsProcessor.GetTforPoint(trackPoint, origin, charge);
  double startZ0 = charge*(trackPoint.GetZ() - startS0 * tTrack);

  if(startZ0 < minZ0 || startZ0 > maxZ0){
    if(config.verbosity>0) cout<<"ERROR -- z0:"<<startZ0<<"\tmin:"<<minZ0<<"\tmax:"<<maxZ0<<endl;
    if(config.requireGoodStartingValues) return nullptr;
  }
  
  // Set calculated initial param values
  SetParameter(fitter, 0, "R0", startR  , minR0   , maxR0 ); // limits from MC
  SetParameter(fitter, 2, "s0", startS0 , minS0   , maxS0 );
  SetParameter(fitter, 4, "L" , startL  , minL    , maxL  );
  SetParameter(fitter, 5, "x0", startX0 , minX0   , maxX0 );
  SetParameter(fitter, 6, "y0", startY0 , minY0   , maxY0 );
  SetParameter(fitter, 7, "z0", startZ0 , minZ0   , maxZ0 );
  
  // With 3 points we don't know how fast will radius and slope decrease:
  FixParameter(fitter, 1, "a" ,  0.0000001);
  FixParameter(fitter, 3, "b" , -0.0000001);
  
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
  
  startX0 = x + sgn(y)*charge*startR/sqrt(pow(x/y, 2)+1);
  startY0 = y - sgn(x)*charge*startR/sqrt(pow(y/x, 2)+1);
  
  if(charge*y < 0)  minX0 -= maxR0/sqrt(pow(x/y, 2)+1);
  else              maxX0 += maxR0/sqrt(pow(x/y, 2)+1);
  
  if(charge*x > 0)  minY0 -= maxR0/sqrt(pow(y/x, 2)+1);
  else              maxY0 += maxR0/sqrt(pow(y/x, 2)+1);
}

void Fitter::RemoveShortHelices(vector<Helix> &helices)
{
  if(config.verbosity>0) cout<<"Removing very short merged helices...";
  vector<Helix> longHelices;
  for(auto helix : helices){
    if(helix.GetNpoints() >= config.trackMinNpoints &&
       helix.GetNlayers() >= config.trackMinNlayers){
      longHelices.push_back(helix);
    }
  }
  if(config.verbosity>0) cout<<" long merged helices: "<<longHelices.size()<<endl;
  helices = longHelices;
}

bool Fitter::LinkAndMergeHelices(vector<Helix> &helices)
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
      
      
      if(nDifferentPoints <= config.mergingMaxDifferentPoints){
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
    
    vector<shared_ptr<Point>> allPoints = helix1.GetPoints();
    vector<double> allPointsT = helix1.GetPointsT();

    alreadyUsedHelices.push_back(iHelix);
    
    vector<int> childIndices;
    
    for(auto iChild : helixLinks[iHelix].second){
      if(find(alreadyUsedHelices.begin(), alreadyUsedHelices.end(), iChild) != alreadyUsedHelices.end()) continue;
      
      alreadyUsedHelices.push_back(iChild);
      
      Helix &helix2 = helices[iChild];
      vector<shared_ptr<Point>> points2 = helix2.GetPoints();
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
    
    if(helix1Copy.GetChi2() < config.trackMaxChi2){
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
  vector<Helix> goodHelices;
  for(int iHelix=0; iHelix<helixLinks.size(); iHelix++){
    if(find(toRemove.begin(), toRemove.end(), iHelix) != toRemove.end()) continue;
    goodHelices.push_back(helixLinks[iHelix].first);
  }
  helices = goodHelices;
  
  return merged;
}

void Fitter::InitLparams()
{
  minL   = layerRanges[nTrackLayers-1].GetMin();
  maxL   = layerRanges[nTrackLayers].GetMax();
  startL = (minL+maxL)/2.;
}
