//
//  Fitter.cpp
//
//  Created by Jeremi Niedziela on 14/12/2018.
//

#include "Fitter.hpp"

Fitter::Fitter() :
vertex(Point(0, 0, 0))
{
//  c1 = new TCanvas("c1","c1",1500,1500);
//  c1->Divide(2,2);
}

Fitter::~Fitter()
{
  
}

vector<unique_ptr<Circle>> Fitter::FitCirclesToPoints(int pxSign, int pySign)
{
  // Prepare 2D projections in XY
  vector<Point> points2D;
  vector<vector<Point>> pointsByLine = pointsProcessor.SplitPointsIntoLines(points, config.linesToleranceForCircles);
  
  for(vector<Point> line : pointsByLine){
    if((int)line.size() >= config.minNpointsAlongZ){
      points2D.push_back(Point(line));
    }
  }
  int nPar=3;
  ROOT::Fit::Fitter *fitter = GetCirclesFitter(pxSign, pySign);
  
  // Store fitted circles for each triplet of points
  vector<unique_ptr<Circle>> circles;
  
  int nPoints = (int)points2D.size();
  double circleThickness = config.circleThickness;
  for(int i=0;i<nPoints;i++){
    for(int j=i+1;j<nPoints;j++){
      for(int k=j+1;k<nPoints;k++){
        
        auto chi2Function = [&](const double *par) {
          double f = 0;
          
          auto circle = circleProcessor.BuildCircleFromParams(par, vertex, track);
          
          f  = pow(circle->GetDistanceToPoint(points2D[i]),2);
          f += pow(circle->GetDistanceToPoint(points2D[j]),2);
          f += pow(circle->GetDistanceToPoint(points2D[k]),2);
          
          return f;
        };
        auto fitFunction = ROOT::Math::Functor(chi2Function, nPar);
        double pStart[nPar];
        fitter->SetFCN(fitFunction, pStart);
        
        if(fitter->FitFCN()) {
          auto result = fitter->Result();
          
          auto circle = circleProcessor.BuildCircleFromParams(result.GetParams(), vertex, track);
          
          int nPointsOnCircle=0;
          for(Point p : points2D){
            if(circle->GetDistanceToPoint(p) < circleThickness) nPointsOnCircle++;
          }
          if(nPointsOnCircle > 3){
            circle->SetPoints(points);
            circles.push_back(move(circle));
          }
        }
      }
    }
  }
  if(circles.size() == 0) return circles;
  circleProcessor.RemoveSimilarCircles(circles);
  return circles;
}

unique_ptr<Helix> Fitter::GetBestFittingHelix(vector<shared_ptr<Point>> _points,
                                              const Track &_track,
                                              const Point &_vertex,
                                              bool drawCircles)
{
  points = _points;
  track = _track;
  vertex = _vertex;
  
  vector<unique_ptr<Circle>> circles = GetAllCirclesForPoints();
  
  if(circles.size() == 0){
    return nullptr;
  }
  if(drawCircles){
    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    c1->cd();
    TH2D *pointsHist = new TH2D("points","points",
                            250, -250, 250,
                            250, -250, 250);
    
    for(auto &p : points){
      pointsHist->Fill(p->GetX(),p->GetY());
    }
    pointsHist->Draw("colz");
    
    for(auto &c : circles){
      auto a = c->GetArc();
      a->SetFillColorAlpha(kWhite, 0.0);
      a->SetLineWidth(1.0);
      a->SetLineColor(kRed);
      a->Draw("sameL");
    }
    c1->Update();
  }
  unique_ptr<Helix> bestHelix = nullptr;
  int maxNregularPoints = 0;
  double maxFractionRegularPoints = 0;
  
  double minPz = config.minPz;
  double maxPz = config.maxPz;
  
  for(auto &circle : circles){
    vector<shared_ptr<Point>> points = circle->GetPoints();
    
    auto testHelix = [&](double pz){
      int charge = 1;
      unique_ptr<Helix> helix = GetHelixFromCircle(circle, pz, charge);
      
      if(IsHelixBetterThanBefore(helix, maxNregularPoints, maxFractionRegularPoints)){
        bestHelix = move(helix);
      }
      
      charge = -1;
      helix = GetHelixFromCircle(circle, pz, charge);
      
      if(IsHelixBetterThanBefore(helix, maxNregularPoints, maxFractionRegularPoints)){
        bestHelix = move(helix);
      }
    };
    
    for(double pz =  maxPz; pz >=  minPz ; pz-=config.stepPz){ testHelix(pz); }
    for(double pz = -maxPz; pz <= -minPz ; pz+=config.stepPz){ testHelix(pz); }
  }
  
  return bestHelix;
}

unique_ptr<Helix> Fitter::GetHelixFromCircle(const unique_ptr<Circle> &circle, double pz, int charge)
{
  unique_ptr<Point> momentum = make_unique<Point>(charge * circle->GetMomentum()->GetX(),
                                                  charge * circle->GetMomentum()->GetY(),
                                                  pz);
  
  unique_ptr<Helix> helix = make_unique<Helix>(circle->GetDecayPoint(), momentum, charge);
  helix->SetPoints(points);
  helixProcessor.CalculateNregularPoints(helix);
  
  return helix;
}

bool Fitter::IsHelixBetterThanBefore(const unique_ptr<Helix> &helix,
                                     int &maxNregularPoints,
                                     double &maxFractionRegularPoints)
{
  int nRegularPoints = helix->GetNregularPoints();
  double fractionRegularPoints = nRegularPoints/(double)helix->GetNpoints();
  
  // Here is a condition to accept new solution as the best one
  // Accept as a new best solution if:
  // - it gives more reqular points than before or,
  // - it gives the same number of regular points, but they counstitute higher fraction of all points than before
  if(nRegularPoints > maxNregularPoints
     || (nRegularPoints == maxNregularPoints && (fractionRegularPoints - maxFractionRegularPoints > 0.001))
     ){
    maxNregularPoints = nRegularPoints;
    maxFractionRegularPoints = fractionRegularPoints;
    return true;
  }
  return false;
}

vector<unique_ptr<Circle>> Fitter::GetAllCirclesForPoints()
{
  // Collect circles for positive charge
  vector<unique_ptr<Circle>> circles, circlesTmp;
  
  circles    = FitCirclesToPoints( 1,  1);
  circlesTmp = FitCirclesToPoints(-1,  1);
  circles.insert(circles.end(), make_move_iterator(circlesTmp.begin()), make_move_iterator(circlesTmp.end()));
  circlesTmp = FitCirclesToPoints( 1, -1);
  circles.insert(circles.end(), make_move_iterator(circlesTmp.begin()), make_move_iterator(circlesTmp.end()));
  circlesTmp = FitCirclesToPoints(-1, -1);
  circles.insert(circles.end(), make_move_iterator(circlesTmp.begin()), make_move_iterator(circlesTmp.end()));
  circleProcessor.RemoveSimilarCircles(circles);
  
  return circles;
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

ROOT::Fit::Fitter* Fitter::GetCirclesFitter(int pxSign, int pySign)
{
  auto pxRange = range<double>(config.minPx, config.maxPx);
  auto pyRange = range<double>(config.minPy, config.maxPy);
  
  double maxR = layerR[track.GetNtrackerLayers()];
  double minR = layerR[track.GetNtrackerLayers()-1];
  
  // Create fitter to fit circles to 2D distribution
  ROOT::Fit::Fitter *fitter = new ROOT::Fit::Fitter();
  
  // This is a stupid hack to be able to set fit parameters before actually setting a fitting function
  // Params we want to set only once, but function will change in each iteration of the loop, because
  // it captures loop iterators.
  auto f = [&](const double*) {return 0;};
  int nPar = 3;
  ROOT::Math::Functor fitFunction = ROOT::Math::Functor(f, nPar);
  double pStart[nPar];
  fitter->SetFCN(fitFunction, pStart);
  
  double helixThickness = config.helixThickness;
  
  SetParameter(fitter, 0, "L", (maxR+minR)/2., minR-helixThickness, maxR+helixThickness);
  SetParameter(fitter, 1, "px",
               pxSign*(pxRange.GetMax()-pxRange.GetMin())/2.,
               pxSign > 0 ? pxRange.GetMin() : -pxRange.GetMax(),
               pxSign > 0 ? pxRange.GetMax() : -pxRange.GetMin());
  
  SetParameter(fitter, 2, "py",
               pySign*(pyRange.GetMax()-pyRange.GetMin())/2.,
               pySign > 0 ? pyRange.GetMin() : -pyRange.GetMax(),
               pySign > 0 ? pyRange.GetMax() : -pyRange.GetMin());
  
  return fitter;
}

vector<unique_ptr<Circle>> Fitter::FitCirclesAndAdjustFirstPoint(TripletsVector &pointTriplets,
                                                                 int pxSign, int pySign,
                                                                 double chi2threshold)
{
  int nPar=3;
  auto fitter = GetCirclesFitter(pxSign, pySign);
  
  vector<unique_ptr<Circle>> circles;
  
  for(auto &p : pointTriplets){
    auto chi2Function = [&](const double *par) {
      double f = 0;
      
      auto circle = circleProcessor.BuildCircleFromParams(par, vertex, track);
      
      f  = pow(circle->GetDistanceToPoint(circle->GetDecayPoint()),2);
      f += pow(circle->GetDistanceToPoint(*p[1]),2);
      f += pow(circle->GetDistanceToPoint(*p[2]),2);
      
      return f;
    };
    
    auto fitFunction = ROOT::Math::Functor(chi2Function, nPar);
    double pStart[nPar];
    fitter->SetFCN(fitFunction, pStart);
    
    if(fitter->FitFCN()) {
      auto result = fitter->Result();
      
      if(result.MinFcnValue() > chi2threshold) continue; // FILTER
      auto circle = circleProcessor.BuildCircleFromParams(result.GetParams(), vertex, track);
      
      p[0]->SetX(circle->GetDecayPoint().GetX());
      p[0]->SetY(circle->GetDecayPoint().GetY());
      p[0]->SetZ(circle->GetDecayPoint().GetZ());
      
      circle->SetPoints(p);
      circle->SetPz(p[1]->GetZ() - p[0]->GetZ());
      
      circles.push_back(move(circle));
    }
  }
  
  return circles;
}


vector<Helix> Fitter::FitHelix(const vector<shared_ptr<Point>> &_points,
                               const Track &_track,
                               const Point &_vertex)
{
  points = _points;
  track = _track;
  vertex = _vertex;
  
  double minPointsSeparation = 10.0; // mm
//  double chi2threshold = 1E-6;
  
  // Get only points that are not too close to each other and plot them
  vector<shared_ptr<Point>> filteredPoints = pointsProcessor.FilterNearbyPoints(points, minPointsSeparation);
  cout<<"Points after cleanup:"<<filteredPoints.size()<<endl;
  
  vector<PointsPair> pointPairs = pointsProcessor.BuildPointPairs(filteredPoints);
  
  vector<Helix> helices;
  
  int iter=0;
  for(auto pointPair : pointPairs){
    if(iter > 0) break;
    Helix helix(track, *pointPair.first, *pointPair.second, vertex);
    
    // check that the radius of a new track candidate is within allowed limits
//    if(helix.R0min > GetRadiusInMagField(config.maxPx, config.maxPy, solenoidField) ||
//       helix.R0max < GetRadiusInMagField(config.minPx, config.minPy, solenoidField)){
//      continue;
//    }
    
    // check that the range of a and b parameters allowed the pion not to gain momentum with time
//    if(helix.amax < 0 || helix.bmin > 0) continue;
    
    helices.push_back(helix);
    
    iter++;
  }

  if(helices.size()==0){
    cout<<"No seeds found..."<<endl;
    return helices;
  }
  cout<<"Seed:"<<endl; helices[0].Print();cout<<endl;
  
  auto seedID = helices[0].seedID;
  bool finished;
  long uniqueID=4872737680;
  
  int iSteps=0;
  int nSteps=1;
  
  do{
    if(iSteps==nSteps) break;
    
    finished = true;
    vector<Helix> helicesAfterExtending;
    
    // for all helices from previous step
    for(Helix &helix : helices){
      
      // that are not yet marked as finished
      if(helix.isFinished){
        helicesAfterExtending.push_back(helix);
      }
      else{
        vector<Helix> extendedHelices;
        
        // try to extend by all possible points
        for(auto &point : filteredPoints){
          
          auto helixCopy = Helix(helix);
          bool extended = helixCopy.ExtendByPoint(*point);

          if(extended){
            
            // check that the radius of a new track candidate is within allowed limits
//            if(helixCopy.R0min > GetRadiusInMagField(config.maxPx, config.maxPy, solenoidField) ||
//               helixCopy.R0max < GetRadiusInMagField(config.minPx, config.minPy, solenoidField)){
//              cout<<"Rejected because of radius"<<endl;
//              continue;
//            }
//            
//            // check that the range of a and b parameters allowed the pion not to gain momentum with time
//            if(helixCopy.amax < 0 || helixCopy.bmin > 0){
//              cout<<"Rejected because of gaining momentum"<<endl;
//              continue;
//            }
            
            if(helix.seedID == seedID){
              cout<<"After extension:"<<endl; helixCopy.Print(); cout<<endl;
            }
            
            helixCopy.uniqueID = uniqueID;
            uniqueID++;
            extendedHelices.push_back(helixCopy);
          }
        }
        // if it was possible to extend the helix
        if(extendedHelices.size() != 0){
          helicesAfterExtending.insert(helicesAfterExtending.end(),
                                       extendedHelices.begin(),
                                       extendedHelices.end());
          finished = false;
        }
        else{ // if helix could not be extended
          if(helix.seedID == seedID){
            cout<<"Could not further extend the helix"<<endl;
          }
          helix.isFinished = true;
          helicesAfterExtending.push_back(helix); // is this needed?
        }
      }
      
    }
    helices.clear();
    
    for(auto &h : helicesAfterExtending){
      helices.push_back(move(h));
    }
    iSteps++;
  }
  while(!finished);
  
  return helices;
  
  // Create all possible point triplets
//  auto pointTriplets = pointsProcessor->BuildPointTriplets(filteredPoints);
//  cout<<"N valid point triplets:"<<pointTriplets.size()<<endl;
//
//
//
//  // Get circles matching those point triplets
//  vector<unique_ptr<Circle>> circles;
//  vector<unique_ptr<Circle>> circlesTmp;
  
//  circlesTmp = FitCirclesAndAdjustFirstPoint(pointTriplets, 1, 1, chi2threshold);
//  circles.insert(circles.end(),
//                 make_move_iterator(circlesTmp.begin()),
//                 make_move_iterator(circlesTmp.end()));
  
//  circlesTmp = FitCirclesAndAdjustFirstPoint(pointTriplets,-1, 1, chi2threshold);
//  circles.insert(circles.end(),
//                 make_move_iterator(circlesTmp.begin()),
//                 make_move_iterator(circlesTmp.end()));
//
//  circlesTmp = FitCirclesAndAdjustFirstPoint(pointTriplets, 1,-1, chi2threshold);
//  circles.insert(circles.end(),
//                 make_move_iterator(circlesTmp.begin()),
//                 make_move_iterator(circlesTmp.end()));
  
//  circlesTmp = FitCirclesAndAdjustFirstPoint(pointTriplets,-1,-1, chi2threshold);
//  circles.insert(circles.end(),
//                 make_move_iterator(circlesTmp.begin()),
//                 make_move_iterator(circlesTmp.end()));
  
//  double minR = layerR[track.GetNtrackerLayers()-1];
//  double maxR = layerR[track.GetNtrackerLayers()];
//
//  double decayXmin = minR*cos(track.GetPhi());
//  double decayYmin = minR*sin(track.GetPhi());
//  double decayZmin = minR/sin(track.GetTheta())*cos(track.GetTheta());
//
//  double decayXmax = maxR*cos(track.GetPhi());
//  double decayYmax = maxR*sin(track.GetPhi());
//  double decayZmax = maxR/sin(track.GetTheta())*cos(track.GetTheta());
//
//  auto originMin = make_shared<Point>(decayXmin, decayYmin, decayZmin);
//  auto originMax = make_shared<Point>(decayXmax, decayYmax, decayZmax);
//
//  double minRadius = GetRadiusInMagField(config.minPx, config.minPy, solenoidField);
//  double maxRadius = GetRadiusInMagField(config.maxPx, config.maxPy, solenoidField);
  
//  TripletPairsVector tripletPairs = pointsProcessor->BuildPointTripletPairs(filteredPoints, originMin, originMax);
//  circles = circleProcessor->BuildCirclesFromTripletPairs(tripletPairs, range<double>(minRadius, maxRadius));
  

  /*
  // Build seeds from the circles (this will check that they make sense)
  vector<unique_ptr<ArcSet2D>> potentialPionTracks = arcSetProcessor->BuildArcSetsFromCircles(circles);
  cout<<"N track seeds:"<<potentialPionTracks.size()<<endl;
  
//  int iTrack = 0;
//  potentialPionTracks.erase(potentialPionTracks.begin(), potentialPionTracks.begin()+iTrack);
//  potentialPionTracks.erase(potentialPionTracks.begin()+1, potentialPionTracks.end());
  
  // fit more segments staring from seeds
  vector<double> alphaVector;
  int iter=0;
  for(auto &pionTrack : potentialPionTracks){
    while(1){
      iter++;
      // Get potential triplets of new points
//      auto newPointTriplets = arcSetProcessor->BuildTripletsCompatibleWithArcSet(pionTrack, filteredPoints);
      auto newPoints = arcSetProcessor->FindPossibleNextPoints(pionTrack, filteredPoints);
      
      vector<unique_ptr<Circle>> newCircles;
      unique_ptr<Circle> bestCircle = nullptr;
      
      double maxRadiiDecrease = 0.2; // max allowed relative decrease in radius
      double maxRadiiIncrease = 0.2; // max allowed relative increase in radius
      
      for(auto point : newPoints){
        auto circle = circleProcessor->GetParallelCircle(pionTrack->GetLastCircle(), point);
        
//        if(circle->GetRadius() > 1.1*pionTrack->GetLastCircle()->GetRadius()) continue;
        
        double diff = (pionTrack->GetLastCircle()->GetRadius() - circle->GetRadius())/pionTrack->GetLastCircle()->GetRadius();
        
        if((diff > 0 &&  diff < maxRadiiDecrease) ||  // radius decreased
           (diff < 0 && -diff < maxRadiiIncrease)     // radius increased
           ){
          newCircles.push_back(move(circle));
        }
      }
      
      double shortestArc = inf;
      
      for(auto &circle : newCircles){
        double arcLength = circle->GetRange().GetMax() - circle->GetRange().GetMin();
        
        if(arcLength < shortestArc){
          shortestArc = arcLength;
          bestCircle = move(circle);
        }
      }
      
      // Build a circle for each possible combination of points
//      vector<unique_ptr<Circle>> newCircles = circleProcessor->BuildCirclesFromPoints(newPointTriplets);
      
//      unique_ptr<Circle> bestCircle = circleProcessor->GetMostCompatibleCircle(newCircles, pionTrack->GetLastCircle(), alphaVector);
      
      if(!bestCircle) break;
      
      // Add best circle's arc and the last point (the new one) to the pion's track candidate
      pionTrack->AddCircle(bestCircle);
    }
  }
  
//  PlotRadiiAngles(alphaVector);
//  PlotSeeds(potentialPionTracks);
//  PlotTracks(potentialPionTracks);
//  PlotGoodTracks(potentialPionTracks);
  PlotRadiiChi2(2, potentialPionTracks);
  
  unique_ptr<ArcSet2D> bestPionTrack = arcSetProcessor->GetBestArcSet(potentialPionTracks);
  
  if(!bestPionTrack) return nullptr;
  
  PlotClusters(1, filteredPoints);
  PlotBestTrack(1, bestPionTrack);
  
  cout<<"------------------------------------------------"<<endl;
  cout<<"The best track is:"<<endl;
  bestPionTrack->Print();
  cout<<"------------------------------------------------"<<endl;
  
  auto helix = make_unique<Helix>(*bestPionTrack->GetOrigin(),
                                  make_unique<Point>(*bestPionTrack->GetCircle(0)->GetMomentum()),
                                  track.GetCharge());
  return helix;
   */
}


void Fitter::PlotSeeds(int iPad, const vector<unique_ptr<ArcSet2D>> &potentialPionTracks)
{
  c1->cd(iPad);
  
  for(auto &pionTrack : potentialPionTracks){
    // skip drawing of seed-only tracks
//    if(pionTrack->GetNarcs() == 1) continue;
    
    TGraph *seedPoints = new TGraph();
    seedPoints->SetPoint(0, pionTrack->GetPoint(0)->GetX(), pionTrack->GetPoint(0)->GetY());
    seedPoints->SetPoint(1, pionTrack->GetPoint(1)->GetX(), pionTrack->GetPoint(1)->GetY());
    seedPoints->SetPoint(2, pionTrack->GetPoint(2)->GetX(), pionTrack->GetPoint(2)->GetY());
    seedPoints->SetMarkerStyle(26);
    seedPoints->SetMarkerSize(1.8);
    seedPoints->SetMarkerColor(kMagenta);
    seedPoints->Draw("P");
    
    TArc *seedArc = pionTrack->GetArcs()[0];
    seedArc->SetFillColorAlpha(kWhite, 0.0);
    seedArc->SetLineWidth(2.0);
    seedArc->SetLineColor(kMagenta);
    seedArc->Draw("sameLonly");
  }
  c1->Update();
}

void Fitter::PlotTracks(int iPad, const vector<unique_ptr<ArcSet2D>> &potentialPionTracks)
{
  c1->cd(iPad);
  auto graphDecay = GetDecayGraph();
  graphDecay->Draw("P");
  
  for(auto &pionTrack : potentialPionTracks){
    pionTrack->Print();
    bool first = true;
    
    for(auto singleArc : pionTrack->GetArcs()){
      if(first){
        first = false;
        continue;
      }
      
      singleArc->SetFillColorAlpha(kWhite, 0.0);
      singleArc->SetLineWidth(1.0);
      singleArc->SetLineColor(kRed);
      singleArc->Draw("sameLonly");
    }
  }
  c1->Update();
}

void Fitter::PlotGoodTracks(int iPad, const vector<unique_ptr<ArcSet2D>> &potentialPionTracks)
{
  c1->cd(iPad);
  auto graphDecay = GetDecayGraph();
  graphDecay->Draw("P");
  
  for(auto &pionTrack : potentialPionTracks){
    if(pionTrack->GetCycle() < 1) continue;
    
    for(auto singleArc : pionTrack->GetArcs()){
      singleArc->SetFillColorAlpha(kWhite, 0.0);
      singleArc->SetLineWidth(1.0);
      singleArc->SetLineColor(kBlue);
      singleArc->Draw("sameLonly");
    }
  }
  c1->Update();
}

void Fitter::PlotRadiiChi2(int iPad, const vector<unique_ptr<ArcSet2D>> &potentialPionTracks)
{
  c1->cd(iPad);
  auto radiiChi2 = new TH1D("radiiChi2","radiiChi2",100,0,1E-12);
  
  for(auto &track : potentialPionTracks){
    radiiChi2->Fill(track->GetRadiiSlopeChi2());
  }
  radiiChi2->Draw();
  c1->Update();
}

void Fitter::PlotRadiiVsIter(int iPad, const vector<unique_ptr<ArcSet2D>> &potentialPionTracks)
{
  c1->cd(iPad);

  vector<int> colors = { kRed, kBlue, kGreen, kMagenta, kCyan, kYellow, kOrange, kAzure };
  int colorIter=0;
  
  auto dummyGraph = new TGraph();
  dummyGraph->SetPoint(0,-1,0);
  dummyGraph->SetPoint(1,10,0);
  dummyGraph->SetPoint(2,0,0);
  dummyGraph->SetPoint(3,0,800);
  dummyGraph->Draw("AP");
  
  for(auto &track : potentialPionTracks){
    auto graph = new TGraph();
    int iter=0;
    for(auto &circle : track->GetCircles()){
      graph->SetPoint(iter, iter, circle->GetRadius());
      iter++;
    }
    graph->SetMarkerSize(1.0);
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(colors[colorIter]);
    graph->SetLineColor(colors[colorIter++]);
    graph->Draw("PLsame");
    
    if(colorIter == colors.size()) colorIter=0;
  }
 
  c1->Update();
}

void Fitter::PlotBestTrack(int iPad, const unique_ptr<ArcSet2D> &pionTrack)
{
  c1->cd(iPad);
  auto graphDecay = GetDecayGraph();
  graphDecay->Draw("P");
  
  for(auto singleArc : pionTrack->GetArcs()){
    singleArc->SetFillColorAlpha(kWhite, 0.0);
    singleArc->SetLineWidth(2.0);
    singleArc->SetLineColor(kGreen);
    singleArc->Draw("sameLonly");
  }
  c1->Update();
}

void Fitter::PlotRadiiAngles(int iPad, const vector<double> &alphaVector)
{
  c1->cd(iPad);
  
  auto radiiAnglesHist = make_unique<TH1D>("radiiAnglesHist","radiiAnglesHist",100,-1,1);
  for(double alpha : alphaVector) radiiAnglesHist->Fill(alpha);
  radiiAnglesHist->Draw();
  
  c1->Update();
}

void Fitter::PlotClusters(int iPad, const vector<shared_ptr<Point>> &filteredPoints)
{
  c1->cd(iPad);
  
  TH2D *pointsHist = new TH2D("points","points",
                              300, -layerR[nLayers-1], layerR[nLayers-1],
                              300, -layerR[nLayers-1], layerR[nLayers-1]);
  
  pointsHist->Fill(0.0, 0.0, 5.0);
  
  for(auto &p : filteredPoints){
    pointsHist->Fill(p->GetX(),p->GetY());
  }
  pointsHist->DrawCopy("colz");
  
  c1->Update();
}

TGraph* Fitter::GetDecayGraph()
{
  auto graphDecay = new TGraph();
  
  double minR = layerR[track.GetNtrackerLayers()-1];
  double maxR = layerR[track.GetNtrackerLayers()];
  
  graphDecay->SetPoint(0,
                       minR*cos(track.GetPhi()) + 10*vertex.GetX(),
                       minR*sin(track.GetPhi()) + 10*vertex.GetY());
  
  graphDecay->SetPoint(1,
                       maxR*cos(track.GetPhi()) + 10*vertex.GetX(),
                       maxR*sin(track.GetPhi()) + 10*vertex.GetY());
  graphDecay->SetMarkerStyle(20);
  graphDecay->SetMarkerSize(1.0);
  graphDecay->SetMarkerColor(kRed);
  
  return graphDecay;
}

ROOT::Fit::Fitter* Fitter::GetHelixParamsFitter(range<double> rangeL)
{
  ROOT::Fit::Fitter *fitter = new ROOT::Fit::Fitter();
  
  // This is a stupid hack to be able to set fit parameters before actually setting a fitting function
  // Params we want to set only once, but function will change in each iteration of the loop, because
  // it captures loop iterators.
  auto f = [&](const double*) {return 0;};
  int nPar = 5;
  ROOT::Math::Functor fitFunction = ROOT::Math::Functor(f, nPar);
  double pStart[nPar];
  fitter->SetFCN(fitFunction, pStart);
  
  double minR = GetRadiusInMagField(config.minPx, config.minPy, solenoidField);
  double maxR = GetRadiusInMagField(config.maxPx, config.maxPy, solenoidField);
  
  SetParameter(fitter, 0, "R0", (minR+maxR)/2., 0, 10000);
  SetParameter(fitter, 1, "a" , (minR+maxR)/2., 0, 10000);
  SetParameter(fitter, 2, "s0", -50, -1000, 1000);
  SetParameter(fitter, 3, "b" , 5000, -10000, 10000);
  SetParameter(fitter, 4, "L" , (rangeL.GetMin()+rangeL.GetMax())/2., rangeL.GetMin(), rangeL.GetMax());
  
  return fitter;
}

HelixParams Fitter::FitHelixParams(const vector<shared_ptr<Point>> &points, const Point &nextPoint,
                                   const Point &origin, const Track &track,
                                   const Point &eventVertex, EHelixParams iParam)
{
  int nLayers = track.GetNtrackerLayers();
  double Lmin = layerR[nLayers-1];
  double Lmax = layerR[nLayers];
  
  int nPar=5;
  auto fitter = GetHelixParamsFitter(range<double>(Lmin, Lmax));
  
  double distPenalty = 10.;
  double testT = 100*TMath::Pi();
  
  double ht = config.helixThickness;
  
  auto chi2Function = [&](const double *par) {
    double R0 = par[0];
    double a  = par[1];
    double s0 = par[2];
    double b  = par[3];
    double L  = par[4];
    
    // First add distance to the vertex
    double vertexX = L * cos(track.GetPhi())    + 10*eventVertex.GetX();
    double vertexY = L * sin(track.GetPhi())    + 10*eventVertex.GetY();
    double vertexZ = L / tan(track.GetTheta())  + 10*eventVertex.GetZ();
    double t = atan2(vertexY - origin.GetY(), vertexX - origin.GetX());
    
    double x = origin.GetX() + (R0 - a*t)*cos(t);
    double y = origin.GetY() + (R0 - a*t)*sin(t);
    double z = origin.GetZ() + (s0 - b*t)*t;
    
    double distX = pow(x-vertexX, 2);
    double distY = pow(y-vertexY, 2);
    double distZ = pow(z-vertexZ, 2);
    
    double f = distX + distY + distZ;
    
    // Then add distances to all other points
    for(auto &p : points){
      t = p->GetT();
      
      // find helix point for this point's t
      x = origin.GetX() + (R0 - a*t)*cos(t);
      y = origin.GetY() + (R0 - a*t)*sin(t);
      z = origin.GetZ() + (s0 - b*t)*t;
      
      // calculate distance between helix and point's boundary (taking into account its errors)
      distX = fabs(x-p->GetX()) > p->GetXerr()+ht ? pow(x-p->GetX(), 2) : 0;
      distY = fabs(y-p->GetY()) > p->GetYerr()+ht ? pow(y-p->GetY(), 2) : 0;
      distZ = fabs(z-p->GetZ()) > p->GetZerr()+ht ? pow(z-p->GetZ(), 2) : 0;
      
      distX /= pow(p->GetXerr(), 2);
      distY /= pow(p->GetYerr(), 2);
      distZ /= pow(p->GetZerr(), 2);
      
      f += distX + distY + distZ;
    }
    
    // Add distance to the testing point
    t = nextPoint.GetT();
    
    // find helix point for this point's t
    x = origin.GetX() + (R0 - a*t)*cos(t);
    y = origin.GetY() + (R0 - a*t)*sin(t);
    z = origin.GetZ() + (s0 - b*t)*t;
    
    // calculate distance between helix and point's boundary (taking into account its errors)
    distX = fabs(x-nextPoint.GetX()) > nextPoint.GetXerr()+ht ? pow(x-nextPoint.GetX(), 2) : 0;
    distY = fabs(y-nextPoint.GetY()) > nextPoint.GetYerr()+ht ? pow(y-nextPoint.GetY(), 2) : 0;
    distZ = fabs(z-nextPoint.GetZ()) > nextPoint.GetZerr()+ht ? pow(z-nextPoint.GetZ(), 2) : 0;
    
    distX /= pow(nextPoint.GetXerr(), 2);
    distY /= pow(nextPoint.GetYerr(), 2);
    distZ /= pow(nextPoint.GetZerr(), 2);
    
    f += distX + distY + distZ;
    
    // apply additional penalty for being outside of the point errors
    f *= distPenalty;
    
    // Finally, add factor minimizing/maximizing slope or radius
    double valR = fabs(R0 - a*testT);
    double valS = fabs(s0 - b*testT);
    
    if(iParam == kMinR) f += valR;
    if(iParam == kMaxR) f += 1/valR;
    if(iParam == kMinS) f += valS;
    if(iParam == kMaxS) f += 1/valS;

    return f;
  };
  
  HelixParams resultParams;
  
  auto fitFunction = ROOT::Math::Functor(chi2Function, nPar);
  double pStart[nPar];
  for(int i=0; i<nPar; i++) pStart[i] = fitter->Config().ParSettings(i).Value();
  fitter->SetFCN(fitFunction, pStart);
  
  
//  if(fitter->FitFCN()) {
  fitter->FitFCN();
  auto result = fitter->Result();
    
  resultParams.R0 = result.GetParams()[0];
  resultParams.a  = result.GetParams()[1];
  resultParams.s0 = result.GetParams()[2];
  resultParams.b  = result.GetParams()[3];
  
  double L = result.GetParams()[4];
  
  // First add distance to the vertex
  double vertexX = L * cos(track.GetPhi())    + 10*eventVertex.GetX();
  double vertexY = L * sin(track.GetPhi())    + 10*eventVertex.GetY();
  double vertexZ = L / tan(track.GetTheta())  + 10*eventVertex.GetZ();
  double t = atan2(vertexY - origin.GetY(), vertexX - origin.GetX());
  
  resultParams.tShift     = t;
  resultParams.tMax       = nextPoint.GetT();
  resultParams.zShift = 0;
  resultParams.zShift = vertexZ - (resultParams.s0 - resultParams.b * t) * t;
  
//  }
  
  return resultParams;
}
