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
  
  // Extend seeds
  cout<<"Extending seeds..."<<endl;
  ExtendSeeds(fittedHelices, pointsByLayer);
  cout<<"Candidates found: "<<fittedHelices.size()<<endl;
  
  // Remove very short candidates which should not even be merged with others
  cout<<"Removing short candidates...";
  vector<Helix> longHelices;
  for(int iHelix=0; iHelix<fittedHelices.size(); iHelix++){
    if(fittedHelices[iHelix].GetPoints().size() >= config.candidateMinNpoints){
      longHelices.push_back(fittedHelices[iHelix]);
    }
  }
  cout<<" Candidates left:"<<longHelices.size()<<endl;
  
  // Merge similar candidates
  cout<<"Merging overlapping helices...";
  while(MergeHelices(longHelices));
  cout<<" merged down to: "<<longHelices.size()<<endl;
  
  // Remove helices that are too short
  cout<<"Removing very short merged helices...";
  vector<Helix> longMergedHelices;
  for(int iHelix=0; iHelix<longHelices.size(); iHelix++){
    if(longHelices[iHelix].GetPoints().size() >= config.trackMinNpoints){
      longMergedHelices.push_back(longHelices[iHelix]);
    }
  }
  cout<<" long merged helices: "<<longMergedHelices.size()<<endl;
  
  cout<<"Refitting surviving helices...";
  for(auto &helix : longMergedHelices){
    if(helix.GetShouldRefit()) RefitHelix(helix);
  }
  cout<<" done."<<endl;
  
  return longMergedHelices;
}

vector<Helix> Fitter::GetSeeds(vector<vector<shared_ptr<Point>>> pointsByLayer)
{
  // find possible middle and last seeds' points
  int trackLayers = track.GetNtrackerLayers();
  cout<<"Track layers:"<<trackLayers<<endl;
  vector<shared_ptr<Point>> possibleMiddlePoints = pointsByLayer[trackLayers];
  vector<shared_ptr<Point>> possibleLastPoints   = pointsByLayer[trackLayers+1];
  Point trackPointMid = pointsProcessor.GetPointOnTrack((layerR[trackLayers-1]+layerR[trackLayers])/2., track, eventVertex);
  
  vector<Helix> seeds;

  int nPairs=0;
  for(auto &middlePoint : possibleMiddlePoints){
    
//    if(fabs(middlePoint->GetX() + 18) < 1.0 &&
//       fabs(middlePoint->GetY() + 275) < 1.0 &&
//       fabs(middlePoint->GetZ() - 75) < 1.0){
//      cout<<"this"<<endl;
//    }
    
    double middleHitDeltaPhi = pointsProcessor.GetPointingAngleXY(Point(0,0,0), trackPointMid, *middlePoint);
    if(middleHitDeltaPhi > config.seedMiddleHitMaxDeltaPhi){
//      cout<<"middle hit phi"<<endl;
      continue;
    }
    
    double middleHitDeltaZ = fabs(middlePoint->GetZ() - trackPointMid.GetZ());
    if(middleHitDeltaZ > config.seedMiddleHitMaxDeltaZ){
//      cout<<"middle hit z"<<endl;
      continue;
    }
    
    
    
    for(auto &lastPoint : possibleLastPoints){
      
//      if(fabs(lastPoint->GetX() -5) < 1.0 &&
//         fabs(lastPoint->GetY() + 323) < 1.0 &&
//         fabs(lastPoint->GetZ() - 29) < 1.0){
//        cout<<"this"<<endl;
//      }
      
      double lastHitDeltaPhi = pointsProcessor.GetPointingAngleXY(trackPointMid, *middlePoint, *lastPoint);
      if(lastHitDeltaPhi > config.seedLastHitMaxDeltaPhi){
//        cout<<"last hit phi"<<endl;
        continue;
      }
      
      double lastPointDeltaZ = fabs(middlePoint->GetZ() - lastPoint->GetZ());
      if(lastPointDeltaZ > config.seedLastHitMaxDeltaZ){
//        cout<<"last hit z"<<endl;
        continue;
      }
      
      nPairs++;
      auto points = { middlePoint, lastPoint };
      auto helix = FitSeed(points,  track.GetCharge());
      
      if(helix){
        helix->SetIncreasing(true); // add decreasing later
        if(helix->GetChi2() < config.seedMaxChi2) seeds.push_back(*helix);
        else{
//          cout<<"Seed chi2 out of limits"<<endl;
        }
      }
    }
  }
  cout<<"Tested pairs: "<<nPairs<<endl;
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
    vertex->SetT(pointsProcessor.GetTforPoint(*vertex, origin, charge));
    
    vector<shared_ptr<Point>> pointsTriplet = { vertex, points[0], points[1] };
    
    for(int iPoint=1;iPoint<pointsTriplet.size(); iPoint++){
      if(charge < 0){
        double t = pointsProcessor.GetTforPoint(*pointsTriplet[iPoint], origin, charge);
        while(t < pointsTriplet[iPoint-1]->GetT()) t += 2*TMath::Pi();
        pointsTriplet[iPoint]->SetT(t);
      }
      else{
        double t = pointsProcessor.GetTforPoint(*pointsTriplet[iPoint], origin, charge);
        while(t > pointsTriplet[iPoint-1]->GetT()) t -= 2*TMath::Pi();
        pointsTriplet[iPoint]->SetT(t);
      }
    }
    
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
    return f/(3*pointsTriplet.size()+6);
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
  resultHelix->SetChi2(result.MinFcnValue());
  
  return resultHelix;
}

void Fitter::ExtendSeeds(vector<Helix> &helices,
                         const vector<vector<shared_ptr<Point>>> &pointsByLayer)
{
  bool finished;
  int nSteps=0;
  do{
    cout<<"Helices before "<<nSteps<<" step: "<<helices.size()<<endl;
    
    finished = true;
    vector<Helix> nextStepHelices;
    
    // for all helices from previous step
    for(Helix &helix : helices){
      
      if(helix.GetIsFinished()){
        nextStepHelices.push_back(helix);
      }
      else{
        // Find points that could extend this helix
        int lastPointLayer = helix.GetLastPoint()->GetLayer();
        vector<shared_ptr<Point>> possiblePoints;
        
        if(lastPointLayer+1 < pointsByLayer.size() && lastPointLayer-1 >= 0){
          if(helix.GetIncreasing())  possiblePoints = pointsByLayer[lastPointLayer+1];
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
          
          if(deltaPhi > config.nextPointMaxDeltaPhi) continue;
          
          double deltaZ = fabs(helix.GetPoints()[nHelixPoints-1]->GetZ() - point->GetZ());
          
          if(deltaZ > config.nextPointMaxDeltaZ) continue;
          
          /// Extend helix by the new point and refit its params
          Helix helixCopy(helix);
          helixCopy.AddPoint(point);
          RefitHelix(helixCopy);
          
          // check if chi2 is small enough
          if(helixCopy.GetChi2() > config.trackMaxChi2) continue;

          extendedHelices.push_back(helixCopy);
        }
        // if it was possible to extend the helix
        if(extendedHelices.size() != 0){
          nextStepHelices.insert(nextStepHelices.end(),
                                 extendedHelices.begin(),
                                 extendedHelices.end());
          finished = false;
        }
        else if(!helix.GetIsFinished()){ // if helix could not be extended
          helix.SetIsFinished(true);
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
    
    unordered_set<shared_ptr<Point>> uniquePoints;
    for(auto &p : points1) uniquePoints.insert(p);
    
    alreadyUsedHelices.push_back(iHelix);
    
    vector<int> childIndices;
    
    for(auto iChild : helixLinks[iHelix].second){
      if(find(alreadyUsedHelices.begin(), alreadyUsedHelices.end(), iChild) != alreadyUsedHelices.end()) continue;
      
      alreadyUsedHelices.push_back(iChild);
      
      Helix &helix2 = helices[iChild];
      vector<shared_ptr<Point>> points2 = helix2.GetPoints();
      for(auto &p : points2) uniquePoints.insert(p);
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
    Point origin(x0, y0, z0);
    auto vertex = make_shared<Point>(pointsProcessor.GetPointOnTrack(L, track, eventVertex));
    vertex->SetT(pointsProcessor.GetTforPoint(*vertex, origin, helix.GetCharge()));
    vector<shared_ptr<Point>> points = helixPoints;
    points[0] = vertex;
    
    for(int iPoint=1;iPoint<points.size(); iPoint++){
      if(helix.GetCharge() < 0){
        double t = pointsProcessor.GetTforPoint(*points[iPoint], origin, helix.GetCharge());
        while(t < points[iPoint-1]->GetT()) t += 2*TMath::Pi();
        points[iPoint]->SetT(t);
      }
      else{
        double t = pointsProcessor.GetTforPoint(*points[iPoint], origin, helix.GetCharge());
        while(t > points[iPoint-1]->GetT()) t -= 2*TMath::Pi();
        points[iPoint]->SetT(t);
      }
    }
    
    double t, distX, distY, distZ, x, y, z;
    double f=0;
    
    // Then add distances to all other points
    for(auto &p : points){
      // find helix point for this point's t
      t = p->GetT();
      x = x0 + (R0 - a*t)*cos(t);
      y = y0 + (R0 - a*t)*sin(t);
      z = -helix.GetCharge()*z0 + (s0 - b*t)*t;
      
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
  SetParameter(fitter, 3, "b" , helix.GetSlopeFactor() , -10000, 10000);
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
  
  int zSign = 0;
  if(points.back()->GetZ() > points.front()->GetZ() > trackPoint.GetZ()) zSign =  1;
  if(points.back()->GetZ() < points.front()->GetZ() < trackPoint.GetZ()) zSign = -1;
  
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
