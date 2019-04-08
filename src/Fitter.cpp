//
//  Fitter.cpp
//
//  Created by Jeremi Niedziela on 14/12/2018.
//

#include "Fitter.hpp"

Fitter::Fitter() :
pointsProcessor(make_unique<PointsProcessor>()),
helixProcessor(make_unique<HelixProcessor>()),
circleProcessor(make_unique<CircleProcessor>()),
arcSetProcessor(make_unique<ArcSetProcessor>()),
vertex(Point(0, 0, 0))
{
  c1 = new TCanvas("c1","c1",1500,1500);
  c1->Divide(2,2);
}

Fitter::~Fitter()
{
  
}

vector<unique_ptr<Circle>> Fitter::FitCirclesToPoints(int pxSign, int pySign)
{
  // Prepare 2D projections in XY
  vector<Point> points2D;
  vector<vector<Point>> pointsByLine = pointsProcessor->SplitPointsIntoLines(points, config->linesToleranceForCircles);
  
  for(vector<Point> line : pointsByLine){
    if((int)line.size() >= config->minNpointsAlongZ){
      points2D.push_back(Point(line));
    }
  }
  int nPar=3;
  ROOT::Fit::Fitter *fitter = GetCirclesFitter(pxSign, pySign);
  
  // Store fitted circles for each triplet of points
  vector<unique_ptr<Circle>> circles;
  
  int nPoints = (int)points2D.size();
  double circleThickness = config->circleThickness;
  for(int i=0;i<nPoints;i++){
    for(int j=i+1;j<nPoints;j++){
      for(int k=j+1;k<nPoints;k++){
        
        auto chi2Function = [&](const double *par) {
          double f = 0;
          
          auto circle = circleProcessor->BuildCircleFromParams(par, vertex, track);
          
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
          
          auto circle = circleProcessor->BuildCircleFromParams(result.GetParams(), vertex, track);
          
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
  circleProcessor->RemoveSimilarCircles(circles);
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
  
  double minPz = config->minPz;
  double maxPz = config->maxPz;
  
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
    
    for(double pz =  maxPz; pz >=  minPz ; pz-=config->stepPz){ testHelix(pz); }
    for(double pz = -maxPz; pz <= -minPz ; pz+=config->stepPz){ testHelix(pz); }
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
  helixProcessor->CalculateNregularPoints(helix);
  
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
  circleProcessor->RemoveSimilarCircles(circles);
  
  return circles;
}

void Fitter::SetParameter(ROOT::Fit::Fitter *fitter, int i, string name, double start, double min, double max, bool fix)
{
  fitter->Config().ParSettings(i).SetName(name);
  fitter->Config().ParSettings(i).SetValue(start);
  fitter->Config().ParSettings(i).SetLimits((min < max) ? min : max,(min < max) ? max : min);
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
  auto pxRange = range<double>(config->minPx, config->maxPx);
  auto pyRange = range<double>(config->minPy, config->maxPy);
  
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
  
  double helixThickness = config->helixThickness;
  
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
      
      auto circle = circleProcessor->BuildCircleFromParams(par, vertex, track);
      
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
      auto circle = circleProcessor->BuildCircleFromParams(result.GetParams(), vertex, track);
      
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


unique_ptr<Helix> Fitter::FitHelix(const vector<shared_ptr<Point>> &_points,
                                   const Track &_track,
                                   const Point &_vertex)
{
  points = _points;
  track = _track;
  vertex = _vertex;
  
  double minPointsSeparation = 10.0; // mm
//  double chi2threshold = 1E-6;
  
  // Get only points that are not too close to each other and plot them
  vector<shared_ptr<Point>> filteredPoints = pointsProcessor->FilterNearbyPoints(points, minPointsSeparation);
  cout<<"Points after cleanup:"<<filteredPoints.size()<<endl;
  
  // Create all possible point triplets
  auto pointTriplets = pointsProcessor->BuildPointTriplets(filteredPoints);
  cout<<"N valid point triplets:"<<pointTriplets.size()<<endl;
  
  // Get circles matching those point triplets
  vector<unique_ptr<Circle>> circles;
  vector<unique_ptr<Circle>> circlesTmp;
  
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
  
  double minR = layerR[track.GetNtrackerLayers()-1];
  double maxR = layerR[track.GetNtrackerLayers()];
  
  double decayXmin = minR*cos(track.GetPhi());
  double decayYmin = minR*sin(track.GetPhi());
  double decayZmin = minR/sin(track.GetTheta())*cos(track.GetTheta());
  
  double decayXmax = maxR*cos(track.GetPhi());
  double decayYmax = maxR*sin(track.GetPhi());
  double decayZmax = maxR/sin(track.GetTheta())*cos(track.GetTheta());
  
  auto originMin = make_shared<Point>(decayXmin, decayYmin, decayZmin);
  auto originMax = make_shared<Point>(decayXmax, decayYmax, decayZmax);
  
  double minRadius = GetRadiusInMagField(config->minPx, config->minPy, solenoidField);
  double maxRadius = GetRadiusInMagField(config->maxPx, config->maxPy, solenoidField);
  
//  TripletPairsVector tripletPairs = pointsProcessor->BuildPointTripletPairs(filteredPoints, originMin, originMax);
//  circles = circleProcessor->BuildCirclesFromTripletPairs(tripletPairs, range<double>(minRadius, maxRadius));
  
  vector<PointsPair> pointPairs = pointsProcessor->BuildPointPairs(filteredPoints);
  
  vector<shared_ptr<Helix>> helices;
  
  for(auto pointPair : pointPairs){
    auto helix = make_shared<Helix>(track, *pointPair.first, *pointPair.second);
    helices.push_back(helix);
  }
  
  for(auto &helix : helices){
    cout<<endl;
    helix->Print();
  }
  
  return nullptr;
  
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
