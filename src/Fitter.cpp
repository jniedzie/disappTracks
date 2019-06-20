//
//  Fitter.cpp
//
//  Created by Jeremi Niedziela on 14/12/2018.
//

#include "Fitter.hpp"

Fitter::Fitter() :
eventVertex(Point(0, 0, 0)),
verbose(true)
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
  if(verbose) cout<<"Looking for seeds..."<<endl;
  vector<Helix> fittedHelices = GetSeeds(pointsByLayer);
  if(verbose) cout<<"Number of valid seeds: "<<fittedHelices.size()<<endl;
  
  // Extend seeds
  if(verbose) cout<<"Extending seeds..."<<endl;
  ExtendSeeds(fittedHelices, pointsByLayer);
  if(verbose) cout<<"Candidates found: "<<fittedHelices.size()<<endl;
  
  // Remove very short candidates which should not even be merged with others
  if(verbose) cout<<"Removing short candidates...";
  vector<Helix> longHelices;
  for(int iHelix=0; iHelix<fittedHelices.size(); iHelix++){
    if(fittedHelices[iHelix].GetNpoints() >= config.candidateMinNpoints){
      longHelices.push_back(fittedHelices[iHelix]);
    }
  }
  if(verbose) cout<<" Candidates left:"<<longHelices.size()<<endl;
  
  // Merge similar candidates
  if(config.mergeFinalHelices){
    if(verbose) cout<<"Merging overlapping helices...";
    while(MergeHelices(longHelices));
    if(verbose) cout<<" merged down to: "<<longHelices.size()<<endl;
  }
  
  // Remove helices that are too short
  if(verbose) cout<<"Removing very short merged helices...";
  vector<Helix> longMergedHelices;
  for(int iHelix=0; iHelix<longHelices.size(); iHelix++){
    if(longHelices[iHelix].GetNpoints() >= config.trackMinNpoints){
      longMergedHelices.push_back(longHelices[iHelix]);
    }
  }
  if(verbose) cout<<" long merged helices: "<<longMergedHelices.size()<<endl;
  
  if(verbose) cout<<"Refitting surviving helices...";
  for(auto &helix : longMergedHelices){
    if(helix.GetShouldRefit()) RefitHelix(helix);
  }
  if(verbose) cout<<" done."<<endl;
  
  return longMergedHelices;
}

vector<Helix> Fitter::GetSeeds(vector<vector<shared_ptr<Point>>> pointsByLayer)
{
  // find possible middle and last seeds' points
  int trackLayers = track.GetNtrackerLayers();
  
  vector<shared_ptr<Point>> possibleMiddlePoints = pointsByLayer[trackLayers];
  vector<shared_ptr<Point>> possibleLastPoints   = pointsByLayer[trackLayers+1];
  Point trackPointMid = pointsProcessor.GetPointOnTrack((layerR[trackLayers-1]+layerR[trackLayers])/2., track, eventVertex);
  
  vector<vector<shared_ptr<Point>>> middlePointsRegrouped = pointsProcessor.RegroupNerbyPoints(possibleMiddlePoints, config.doubleHitsMaxDistance);
  vector<vector<shared_ptr<Point>>> lastPointsRegrouped   = pointsProcessor.RegroupNerbyPoints(possibleLastPoints, config.doubleHitsMaxDistance);
  
  vector<Helix> seeds;

  int nPairs=0;

  for(auto &middlePoints : middlePointsRegrouped){
    
    vector<shared_ptr<Point>> goodMiddlePoints;
    
    for(auto &point : middlePoints){
      
      double middleHitDeltaPhi = pointsProcessor.GetPointingAngleXY(Point(0,0,0), trackPointMid, *point);
      
      if(config.doAsymmetricConstraints)  middleHitDeltaPhi = track.GetCharge()*middleHitDeltaPhi;
      else                                middleHitDeltaPhi = fabs(middleHitDeltaPhi);
      
      if(config.seedMiddleHitDeltaPhi.IsOutside(middleHitDeltaPhi)){
        if(verbose) cout<<"Seed middle hit Δφ too large"<<endl;
        continue;
      }

      double middleHitDeltaZ = fabs(point->GetZ() - trackPointMid.GetZ());
      if(middleHitDeltaZ > config.seedMiddleHitMaxDeltaZ){
        if(verbose) cout<<"Seed middle hit Δz too large"<<endl;
        continue;
      }
        
      goodMiddlePoints.push_back(point);
    }
    if(goodMiddlePoints.size()==0) continue;
    
    for(auto &lastPoints : lastPointsRegrouped){
      
      vector<shared_ptr<Point>> goodLastPoints;
      
      for(auto &point : lastPoints){
        for(auto &middlePoint : goodMiddlePoints){
          double lastHitDeltaPhi = pointsProcessor.GetPointingAngleXY(trackPointMid, *middlePoint, *point);
          
          if(config.doAsymmetricConstraints)  lastHitDeltaPhi = track.GetCharge()*lastHitDeltaPhi;
          else                                lastHitDeltaPhi = fabs(lastHitDeltaPhi);
          
          if(config.seedLastHitDeltaPhi.IsOutside(lastHitDeltaPhi)){
            if(verbose) cout<<"Seed last hit Δφ too large"<<endl;
            continue;
          }
          
          double lastPointDeltaZ = fabs(middlePoint->GetZ() - point->GetZ());
          if(lastPointDeltaZ > config.seedLastHitMaxDeltaZ){
            if(verbose) cout<<"Seed last hit Δz too large"<<endl;
            continue;
          }
          
          goodLastPoints.push_back(point);
          break;
        }
      }
      if(goodLastPoints.size()==0) continue;
      
      nPairs++;
      vector<shared_ptr<Point>> points;
      points.insert(points.end(), goodMiddlePoints.begin(), goodMiddlePoints.end());
      points.insert(points.end(), goodLastPoints.begin(), goodLastPoints.end());

      auto helix = FitSeed(points,  track.GetCharge());
      
      if(helix){
        if(helix->GetChi2() < config.seedMaxChi2) seeds.push_back(*helix);
      }
    }
  }
  if(verbose) cout<<"Tested pairs: "<<nPairs<<endl;
  return seeds;
}

unique_ptr<Helix> Fitter::FitSeed(const vector<shared_ptr<Point>> &points, int charge)
{
  auto fitter = GetSeedFitter(points);
  
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
    Point origin(x0, y0, z0);
    auto vertex = make_shared<Point>(pointsProcessor.GetPointOnTrack(L, track, eventVertex));

    vector<shared_ptr<Point>> pointsTriplet = { vertex };
    pointsTriplet.insert(pointsTriplet.end(), points.begin(), points.end());
    pointsProcessor.SetPointsT(pointsTriplet, origin, charge);
    
    double t, distX, distY, distZ, x, y, z;
    double f=0;

    // Then add distances to all other points
    for(auto &p : pointsTriplet){
      t = p->GetT();
      x = x0 + (R0 - a*t)*cos(t);
      y = y0 + (R0 - a*t)*sin(t);
      z = -charge*z0 + (s0 - b*t)*t;
      
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
    return f/(3*pointsTriplet.size()+6); // 6 is number of free fit parameters
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
  resultHelix->SetChi2(result.MinFcnValue());
  
  return resultHelix;
}

void Fitter::ExtendSeeds(vector<Helix> &helices,
                         const vector<vector<shared_ptr<Point>>> &pointsByLayer)
{
  bool finished;
  int nSteps=0;
  do{
    if(verbose) cout<<"Helices before "<<nSteps<<" step: "<<helices.size()<<endl;
    
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
            if(find(helixPoints.begin(), helixPoints.end(), pointCandidate) == helixPoints.end()){
              possiblePoints.push_back(pointCandidate);
            }
          }
          
          vector<vector<shared_ptr<Point>>> possiblePointsRegrouped = pointsProcessor.RegroupNerbyPoints(possiblePoints, config.doubleHitsMaxDistance);
          
          for(auto &points : possiblePointsRegrouped){
            
            vector<shared_ptr<Point>> goodPoints;
            
            // Set T of possible new pionts like if they were laying on this helix
            vector<shared_ptr<Point>> tmpPoints;
            tmpPoints.insert(tmpPoints.end(), helixPoints.begin(), helixPoints.end());
            tmpPoints.insert(tmpPoints.end(), points.begin(), points.end());
            pointsProcessor.SetPointsT(tmpPoints, helix.GetOrigin(), helix.GetCharge());
            
            for(auto &point : points){
              if(!helixProcessor.IsPointCloseToHelixInLayer(helix, *point, nextPointLayer, true)) continue;
              if(!pointsProcessor.IsZgood(lastPoints, point)) continue;
              if(!pointsProcessor.IsTgood(lastPoints, point)) continue;
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
          if(fabs(lastPoints.front()->GetX()-130) < 1 &&
             fabs(lastPoints.front()->GetY()-666) < 1 &&
             fabs(lastPoints.front()->GetZ()+94) < 1){
          
          }
          
            // try to extend to the same layer (turning back)
          vector<shared_ptr<Point>> turningBackPointsAll = pointsByLayer[lastPointLayer];
          vector<shared_ptr<Point>> turningBackPoints;
          
          // remove points that are already on helix
          for(auto &pointCandidate : turningBackPointsAll){
            if(find(helixPoints.begin(), helixPoints.end(), pointCandidate) == helixPoints.end()){
              turningBackPoints.push_back(pointCandidate);
            }
          }
          vector<vector<shared_ptr<Point>>> turningBackPointsRegrouped = pointsProcessor.RegroupNerbyPoints(turningBackPoints, config.doubleHitsMaxDistance);
          
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
            helixCopy.SetFirstTurningPointIndex(helixCopy.GetNpoints());
            for(auto &point : goodTurningBackPoints) helixCopy.AddPoint(point);
            RefitHelix(helixCopy);
            
            // check if chi2 is small enough
            if(helixCopy.GetChi2() > config.trackMaxChi2) continue;
            
            // if we reached this point, it means that this hit is not missing
            helixCopy.SetIsPreviousHitMissing(false);
            
            
            extendedHelices.push_back(helixCopy);
          }
        
          if(config.mergeAtTurnBack) while(MergeHelices(extendedHelices));
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
            if(helix.GetLastPoints().size()==1 && helix.GetLastPoints().front()->GetSubDetName()=="missing"){
              helix.RemoveLastPoint();
              RefitHelix(helix);
            }
            
            helix.SetIsFinished(true);
            nextStepHelices.push_back(helix);
          }
          // add a missing hit
          else{
            int missingHitLayer = helix.GetFirstTurningPointIndex() > 0 ? lastPointLayer-1 : lastPointLayer+1;
            shared_ptr<Point> missingHit = helixProcessor.GetPointCloseToHelixInLayer(helix, missingHitLayer);

            helix.AddPoint(missingHit); // this will set missing hit's T
            double t = missingHit->GetT();
            missingHit->SetZ(-helix.GetCharge()*helix.GetOrigin().GetZ() + helix.GetSlope(t)*t + 10*eventVertex.GetZ());
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
}

bool Fitter::MergeHelices(vector<Helix> &helices)
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
    Helix &helix1 = helixLinks[iHelix].first;
    Helix helix1Copy(helix1);
    
    vector<shared_ptr<Point>> points1 = helix1.GetPoints();
    
    unordered_set<shared_ptr<Point>> uniquePoints = { make_shared<Point>(0,0,0) }; // add dummy decay vertex (layer=-1)
    for(auto &p : points1){
      if(p->GetLayer()<0) continue; // don't merge in the decay vertex
      uniquePoints.insert(p);
    }
    alreadyUsedHelices.push_back(iHelix);
    
    vector<int> childIndices;
    
    for(auto iChild : helixLinks[iHelix].second){
      if(find(alreadyUsedHelices.begin(), alreadyUsedHelices.end(), iChild) != alreadyUsedHelices.end()) continue;
      
      alreadyUsedHelices.push_back(iChild);
      
      Helix &helix2 = helices[iChild];
      vector<shared_ptr<Point>> points2 = helix2.GetPoints();
      for(auto &p : points2){
        if(p->GetLayer()<0) continue; // don't merge in the decay vertex
        uniquePoints.insert(p);
      }
      
      toRemove.push_back(iChild);
      childIndices.push_back(iChild);
    }
    vector<shared_ptr<Point>> allPoints(uniquePoints.begin(), uniquePoints.end());
    helix1Copy.SetPoints(allPoints);
    
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


void Fitter::RefitHelix(Helix &helix)
{
  vector<shared_ptr<Point>> helixPoints = helix.GetPoints();
  
  auto chi2Function = [&](const double *par) {
    double R0 = par[0];
    double a  = par[1];
    double s0 = par[2];
    double b  = par[3];
    double L  = par[4];
    double x0 = par[5];
    double y0 = par[6];
    double z0 = par[7];
    
    // Create new origin
    Point origin(x0, y0, z0);
    
    // Create new vertex
    auto vertex = make_shared<Point>(pointsProcessor.GetPointOnTrack(L, track, eventVertex));
    
    // Copy points before turning, set new vertex and update T params
    vector<shared_ptr<Point>> points = helixPoints;
    points[0] = vertex;
    pointsProcessor.SetPointsT(points, origin, helix.GetCharge());
    
    double t, distX, distY, distZ, x, y, z;
    double f=0;
    
    // Then add distances to all other points
    for(auto &p : points){
      if(p->GetSubDetName()=="missing") continue;
      
      // find helix point for this point's t
      t = p->GetT();
      x = x0 + (R0 - a*t)*cos(t);
      y = y0 + (R0 - a*t)*sin(t);
      z = -helix.GetCharge()*z0 + (s0 - b*t)*t;
      
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
    return f/(3*helixPoints.size()+8);
  };
  
  ROOT::Fit::Fitter *fitter = new ROOT::Fit::Fitter();
  int nPar = 8;
  ROOT::Math::Functor fitFunction = ROOT::Math::Functor(chi2Function, nPar);
  double pStart[nPar];
  fitter->SetFCN(fitFunction, pStart);
  
  double Lmin = layerRanges[track.GetNtrackerLayers()-1].GetMax();
  double Lmax = layerRanges[track.GetNtrackerLayers()].GetMin();
  
  SetParameter(fitter, 0, "R0", helix.GetRadius(0), 0, 1000); // from MC
  SetParameter(fitter, 1, "a" , helix.GetRadiusFactor() , 0, 10000);
  SetParameter(fitter, 2, "s0", helix.GetSlope(0), -10000, 10000);
  SetParameter(fitter, 3, "b" , helix.GetSlopeFactor() , -10000, 0);
  SetParameter(fitter, 4, "L" , (helix.GetVertex()->GetX()-10*eventVertex.GetX())/cos(track.GetPhi()),
               Lmin, Lmax);
  SetParameter(fitter, 5, "x0", helix.GetOrigin().GetX(), -1000 , 1000);
  SetParameter(fitter, 6, "y0", helix.GetOrigin().GetY(), -1000 , 1000);
  SetParameter(fitter, 7, "z0", helix.GetOrigin().GetZ(), -1000 , 1000);
  
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
  
  helix.SetParams(resultParams);
  helix.SetVertex(vertex);
  helix.UpdateOrigin(origin);
  helix.SetChi2(result.MinFcnValue());
  
}

ROOT::Fit::Fitter* Fitter::GetSeedFitter(const vector<shared_ptr<Point>> &points)
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
  
  
  // Calculate initial parameters as good as we can at this point.
  double startR = 320; // mm, from MC
  double maxR = 1000;
  
  // -- decay vertex must be after last track layer and before the next one
  double minL = layerRanges[track.GetNtrackerLayers()-1].GetMax();
  double maxL = layerRanges[track.GetNtrackerLayers()].GetMin();
  double startL = (minL+maxL)/2.; // estimate decay vertex in between of the two above
  Point trackPoint = pointsProcessor.GetPointOnTrack(startL, track, eventVertex);
  
  // -- calculate where wuold the origin X and Y be for the most probable radius and average L position
  double startX0 = 0, minX0 = trackPoint.GetX(), maxX0 = trackPoint.GetX();
  double startY0 = 0, minY0 = trackPoint.GetY(), maxY0 = trackPoint.GetY();
  
  if(trackPoint.GetX() >= 0 && trackPoint.GetY() > 0){
    startX0 = trackPoint.GetX() + track.GetCharge() * startR/sqrt(pow(trackPoint.GetX()/trackPoint.GetY(), 2)+1);
    minX0 -= track.GetCharge() > 0 ? 0 : maxR/sqrt(pow(trackPoint.GetX()/trackPoint.GetY(), 2)+1);
    maxX0 += track.GetCharge() > 0 ? maxR/sqrt(pow(trackPoint.GetX()/trackPoint.GetY(), 2)+1) : 0;
    
    startY0 = trackPoint.GetY() - track.GetCharge() * startR/sqrt(pow(trackPoint.GetY()/trackPoint.GetX(), 2)+1);
    minY0 -= track.GetCharge() > 0 ? maxR/sqrt(pow(trackPoint.GetY()/trackPoint.GetX(), 2)+1) : 0;
    maxY0 += track.GetCharge() > 0 ? 0 : maxR/sqrt(pow(trackPoint.GetY()/trackPoint.GetX(), 2)+1);
  }
  if(trackPoint.GetX() < 0 && trackPoint.GetY() > 0){
    startX0 = trackPoint.GetX() + track.GetCharge() * startR/sqrt(pow(trackPoint.GetX()/trackPoint.GetY(), 2)+1);
    minX0 -= track.GetCharge() > 0 ? 0 : maxR/sqrt(pow(trackPoint.GetX()/trackPoint.GetY(), 2)+1);
    maxX0 += track.GetCharge() > 0 ? maxR/sqrt(pow(trackPoint.GetX()/trackPoint.GetY(), 2)+1) : 0;
    
    startY0 = trackPoint.GetY() + track.GetCharge() * startR/sqrt(pow(trackPoint.GetY()/trackPoint.GetX(), 2)+1);
    minY0 -= track.GetCharge() > 0 ? 0 : maxR/sqrt(pow(trackPoint.GetY()/trackPoint.GetX(), 2)+1);
    maxY0 += track.GetCharge() > 0 ? maxR/sqrt(pow(trackPoint.GetY()/trackPoint.GetX(), 2)+1) : 0;
  }
  if(trackPoint.GetX() >= 0 && trackPoint.GetY() <= 0){
    startX0 = trackPoint.GetX() - track.GetCharge() * startR/sqrt(pow(trackPoint.GetX()/trackPoint.GetY(), 2)+1);
    minX0 -= track.GetCharge() > 0 ? maxR/sqrt(pow(trackPoint.GetX()/trackPoint.GetY(), 2)+1) : 0;
    maxX0 += track.GetCharge() > 0 ? 0 : maxR/sqrt(pow(trackPoint.GetX()/trackPoint.GetY(), 2)+1);

    startY0 = trackPoint.GetY() - track.GetCharge() * startR/sqrt(pow(trackPoint.GetY()/trackPoint.GetX(), 2)+1);
    minY0 -= track.GetCharge() > 0 ? maxR/sqrt(pow(trackPoint.GetY()/trackPoint.GetX(), 2)+1) : 0;
    maxY0 += track.GetCharge() > 0 ? 0 : maxR/sqrt(pow(trackPoint.GetY()/trackPoint.GetX(), 2)+1);
  }
  if(trackPoint.GetX() < 0 && trackPoint.GetY() <= 0){
    startX0 = trackPoint.GetX() - track.GetCharge() * startR/sqrt(pow(trackPoint.GetX()/trackPoint.GetY(), 2)+1);
    minX0 -= track.GetCharge() > 0 ? maxR/sqrt(pow(trackPoint.GetX()/trackPoint.GetY(), 2)+1) : 0;
    maxX0 += track.GetCharge() > 0 ? 0 : maxR/sqrt(pow(trackPoint.GetX()/trackPoint.GetY(), 2)+1);
    
    startY0 = trackPoint.GetY() + track.GetCharge() * startR/sqrt(pow(trackPoint.GetY()/trackPoint.GetX(), 2)+1);
    minY0 -= track.GetCharge() > 0 ? 0 : maxR/sqrt(pow(trackPoint.GetY()/trackPoint.GetX(), 2)+1);
    maxY0 += track.GetCharge() > 0 ? maxR/sqrt(pow(trackPoint.GetY()/trackPoint.GetX(), 2)+1) : 0;
  }
  
  // -- calculate slope from the track momentum direction (pion usually follows this direction)
  double startS0 = startR * trackPoint.GetVectorSlopeC();
  
  // -- get t param of the track point
  Point origin(startX0, startY0, 0);
  double tTrack = pointsProcessor.GetTforPoint(trackPoint, origin, track.GetCharge());
  
  // -- calculate Z position of the vertex
  double startZ0 = -track.GetCharge() * (trackPoint.GetZ() - startS0 * tTrack);
  double minZ0 = -1000; // to be adjusted from math of MC
  double maxZ0 = 1000;
  
  if(startX0 < minX0 || startX0 > maxX0){
    cout<<"ERROR -- x0:"<<startX0<<"\tmin:"<<minX0<<"\tmax:"<<maxX0<<endl;
  }
  if(startY0 < minY0 || startY0 > maxY0){
    cout<<"ERROR -- y0:"<<startY0<<"\tmin:"<<minY0<<"\tmax:"<<maxY0<<endl;
  }
  if(startZ0 < minZ0 || startZ0 > maxZ0){
    cout<<"ERROR -- z0:"<<startZ0<<"\tmin:"<<minZ0<<"\tmax:"<<maxZ0<<endl;
  }
  
  // Set calculated initial param values
  SetParameter(fitter, 0, "R0", startR  ,  0      , maxR  ); // limits from MC
  SetParameter(fitter, 2, "s0", startS0 , -10000  , 10000 );
  SetParameter(fitter, 4, "L" , startL  ,  minL   , maxL  );
  SetParameter(fitter, 5, "x0", startX0 , minX0   , maxX0 );
  SetParameter(fitter, 6, "y0", startY0 , minY0   , maxY0 );
  SetParameter(fitter, 7, "z0", startZ0 , minZ0   , maxZ0  );
  
  // With 3 points we don't know how fast will radius and slope decrease:
  FixParameter(fitter, 1, "a" , 0);
  FixParameter(fitter, 3, "b" , 0);
  
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
