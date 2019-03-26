//
//  Fitter.cpp
//
//  Created by Jeremi Niedziela on 14/12/2018.
//

#include "Fitter.hpp"

Fitter::Fitter() :
pointsProcessor(make_unique<PointsProcessor>()),
helixProcessor(make_unique<HelixProcessor>()),
circleProcessor(make_unique<CircleProcessor>())
{
  c1 = new TCanvas("c1","c1",2000,1000);
  c1->Divide(2,1);
  c1->cd(1);
  
  radiiAnglesHist = new TH1D("radiiAnglesHist","radiiAnglesHist",100,-1,1);
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
          
          auto circle = GetCircleFromFitterParams(par);
          
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
          
          auto circle = GetCircleFromFitterParams(result.GetParams());
          
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
                                              const shared_ptr<Track> _track,
                                              const unique_ptr<Point> &_vertex,
                                              bool drawCircles)
{
  points = _points;
  track = _track;
  vertex = make_unique<Point>(_vertex);
  
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
  
  double minL = layerR[track->GetLastBarrelLayer()];
  double maxL = layerR[track->GetLastBarrelLayer()+1];
  
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
  
  SetParameter(fitter, 0, "L", (maxL+minL)/2., minL-helixThickness, maxL+helixThickness);
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

vector<vector<shared_ptr<Point>>> Fitter::BuildPointTriplets(const vector<shared_ptr<Point>> &inputPoints)
{
  int nPoints = (int)inputPoints.size();
  vector<vector<shared_ptr<Point>>> pointTriplets;
  
  for(int i=0;i<nPoints;i++){
    for(int j=i+1;j<nPoints;j++){
      vector<shared_ptr<Point>> points = {make_shared<Point>(inf,inf,inf), inputPoints[i], inputPoints[j] };
      pointTriplets.push_back(points);
    }
  }
  return pointTriplets;
}

vector<vector<shared_ptr<Point>>> Fitter::BuildPointTriplets(const unique_ptr<ArcSet2D> &pionTrack,
                                                             const vector<shared_ptr<Point>> &inputPoints)
{
  vector<vector<shared_ptr<Point>>> pointTriplets;
  
  double stripSensorHalfLength = stripModuleZlength/2.;
  
  unique_ptr<Circle> circle = pionTrack->GetLastCircle();
  shared_ptr<Point> point1  = pionTrack->GetSecondToLastPoint();
  shared_ptr<Point> point2  = pionTrack->GetLastPoint();
  vector<shared_ptr<Point>> pionPoints = pionTrack->GetPoints();
  
  for(auto point : inputPoints){
    bool isValidPoint = true;
    
    // make sure that it's not the same point as already in the pion track candidate
    if(find(pionPoints.begin(), pionPoints.end(), point) != pionPoints.end()) continue; // FILTER
    
    // new point must be on the correct side of the previous arc in Z direction
    if(fabs(point1->GetZ() - point2->GetZ()) > stripSensorHalfLength){
      // we can check it only if two prevous hits are in different Z locations
      
      if(point1->GetZ() + stripSensorHalfLength < point2->GetZ() - stripSensorHalfLength  &&
         point->GetZ() < (point2->GetZ() - stripSensorHalfLength)){
        isValidPoint = false; // FILTER
      }
      
      if(point1->GetZ() - stripSensorHalfLength > point2->GetZ() + stripSensorHalfLength  &&
         point->GetZ() > (point2->GetZ() + stripSensorHalfLength)){
        isValidPoint = false; // FILTER
      }
    }
    
    // it also has to be within the radius of the helix
    double pointR = pointsProcessor->distanceXY(point, circle->GetCenter());
    if(pointR > 1.1*circle->GetRadius()) isValidPoint = false; // FILTER
    
    if(isValidPoint){
      vector<shared_ptr<Point>> points = {point1, point2, point};
      pointTriplets.push_back(points);
    }
  }
  return pointTriplets;
}

unique_ptr<Circle> Fitter::GetCircleFromFitterParams(const double *par)
{
  double L  = par[0];
  double px = par[1];
  double py = par[2];
  
  double x0 = L*cos(track->GetPhi())                          + 10*vertex->GetX();
  double y0 = L*sin(track->GetPhi())                          + 10*vertex->GetY();
  double z0 = L/sin(track->GetTheta())*cos(track->GetTheta()) + 10*vertex->GetZ();
  
  auto decayPoint  = make_unique<Point>(x0,y0,z0);
  auto momentum    = make_unique<Point>(px,py,0);
  
  return make_unique<Circle>(decayPoint, momentum);
}

vector<unique_ptr<Circle>> Fitter::FitCirclesAndAdjustFirstPoint(vector<vector<shared_ptr<Point>>> &pointTriplets,
                                                                 int pxSign, int pySign,
                                                                 double chi2threshold)
{
  int nPar=3;
  auto fitter = GetCirclesFitter(pxSign, pySign);
  
  vector<unique_ptr<Circle>> circles;
  
  for(auto &p : pointTriplets){
    auto chi2Function = [&](const double *par) {
      double f = 0;
      
      auto circle = GetCircleFromFitterParams(par);
      
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
      auto circle = GetCircleFromFitterParams(result.GetParams());
      
      p[0]->SetX(circle->GetDecayPoint()->GetX());
      p[0]->SetY(circle->GetDecayPoint()->GetY());
      p[0]->SetZ(circle->GetDecayPoint()->GetZ());
      
      circle->SetPoints(p);
      circle->SetPz(p[1]->GetZ() - p[0]->GetZ());
      
      circles.push_back(move(circle));
    }
  }
  
  return circles;
}

bool Fitter::IsValidSeed(const unique_ptr<Circle> &circle, vector<shared_ptr<Point>> pointTriplet)
{
  double phiVertex = circle->GetPointAngle(pointTriplet[0]->GetX(), pointTriplet[0]->GetY());
  double phi1      = circle->GetPointAngle(pointTriplet[1]->GetX(), pointTriplet[1]->GetY());
  double phi2      = circle->GetPointAngle(pointTriplet[2]->GetX(), pointTriplet[2]->GetY());
  
  // Reject cases where hits are on the both sides of the vertex point instead of forming a tracklet
  if((phi1 < phiVertex && phi2 > phiVertex) ||
     (phi2 < phiVertex && phi1 > phiVertex)){
    return false; // FILTER
  }
  
  double z0 = pointTriplet[0]->GetZ();
  double z1 = pointTriplet[1]->GetZ();
  double z2 = pointTriplet[2]->GetZ();
  
  // this is the range of where point 0 and point 1 can be located along Z axis
  double z0min = z0 - stripModuleZlength/2.;
  double z0max = z0 + stripModuleZlength/2.;
  double z1min = z1 - stripModuleZlength/2.;
  double z1max = z1 + stripModuleZlength/2.;
  
  // which determines range in the slope:
  double slopeMin = (z1min - z0max)/phi1;
  double slopeMax = (z1max - z0min)/phi1;
  
  // point 2 can be located somewhere between:
  double z2min = z0max + slopeMin * phi2;
  double z2max = z0min + slopeMax * phi2;
  
  if((z2 < z2min) ||
     (z2 > z2max)){
    // this condition kills correct seed... has to be re-thought
//    return false; // FILTER
  }
  
  return true;
}

range<double> Fitter::GetPhiRange(const unique_ptr<Circle> &circle,
                                  const unique_ptr<ArcSet2D> &pionTrack)
{
  vector<shared_ptr<Point>> pointTriplet = circle->GetPoints();
  
  double phiVertex = circle->GetPointAngle(pointTriplet[0]->GetX(), pointTriplet[0]->GetY());
  double phi1      = circle->GetPointAngle(pointTriplet[1]->GetX(), pointTriplet[1]->GetY());
  double phi2      = circle->GetPointAngle(pointTriplet[2]->GetX(), pointTriplet[2]->GetY());
  
  bool clockwise = pionTrack->IsClockwise();
  if(phi2 > phiVertex &&
     phi2 > phi1){
    pionTrack->IncreaseCycle();
    
    phi2 = phi2 + (clockwise ? -1 : 1) * pionTrack->GetCycle() * 2*TMath::Pi();
  }
  else if(phi2 > phiVertex &&
          phi1 > phiVertex){
    phi1 = phi1 + (clockwise ? -1 : 1) * pionTrack->GetCycle() * 2*TMath::Pi();
    phi2 = phi2 + (clockwise ? -1 : 1) * pionTrack->GetCycle() * 2*TMath::Pi();
  }
  else{
    phiVertex = phiVertex + (clockwise ? -1 : 1) * pionTrack->GetCycle() * 2*TMath::Pi();
    phi1      = phi1 + (clockwise ? -1 : 1) * pionTrack->GetCycle() * 2*TMath::Pi();
    phi2      = phi2 + (clockwise ? -1 : 1) * pionTrack->GetCycle() * 2*TMath::Pi();
  }
  
  double phiMin = min(min(phi1, phi2), phiVertex)/TMath::Pi() * 180;
  double phiMax = max(max(phi1, phi2), phiVertex)/TMath::Pi() * 180;
  
  return range<double>(phiMin, phiMax);
}

range<double> Fitter::GetSeedPhiRange(const unique_ptr<Circle> &circle,
                                      vector<shared_ptr<Point>> pointTriplet,
                                      bool &clockwise)
{
  double phiVertex = circle->GetPointAngle(pointTriplet[0]->GetX(), pointTriplet[0]->GetY());
  double phi1      = circle->GetPointAngle(pointTriplet[1]->GetX(), pointTriplet[1]->GetY());
  double phi2      = circle->GetPointAngle(pointTriplet[2]->GetX(), pointTriplet[2]->GetY());
  
  clockwise = true;
//  if(phi1 < phiVertex)  clockwise = true;
//  else                  clockwise = false;
  
  double phiMin = min(min(phi1, phi2), phiVertex)/TMath::Pi() * 180;
  double phiMax = max(max(phi1, phi2), phiVertex)/TMath::Pi() * 180;
  
  return range<double>(phiMin, phiMax);
}

vector<unique_ptr<ArcSet2D>> Fitter::BuildArcSetsFromCircles(const vector<unique_ptr<Circle>> &circles,
                                                             vector<vector<shared_ptr<Point>>> pointTriplets)
{
  vector<unique_ptr<ArcSet2D>> arcs;
  int iter=0;
  
  for(auto &circle : circles){
    
    if(!IsValidSeed(circle, pointTriplets[iter])) continue;
    
    bool clockwise;
    auto circleRange = GetSeedPhiRange(circle, pointTriplets[iter], clockwise);
    
    auto arcSet2D = make_unique<ArcSet2D>(clockwise);
    arcSet2D->AddCircle(circle, circleRange);
    arcSet2D->AddPoints(pointTriplets[iter]);
    arcs.push_back(move(arcSet2D));
    
//    circle->Print();
    
    iter++;
  }
  
  return arcs;
}

unique_ptr<Circle> Fitter::GetBestCircle(const vector<unique_ptr<Circle>> &newCircles,
                                         const unique_ptr<Circle> &previousCircle)
{
  double bestRadiiDifference = inf;
  unique_ptr<Circle> bestCircle = nullptr;
  
  double p_x = previousCircle->GetPoints()[2]->GetX();
  double p_y = previousCircle->GetPoints()[2]->GetY();
  double c1_x = previousCircle->GetCenter()->GetX();
  double c1_y = previousCircle->GetCenter()->GetY();
  
  double r1_x = c1_x - p_x;
  double r1_y = c1_y - p_y;
  double r1_mod = sqrt(r1_x*r1_x + r1_y*r1_y);
  
  for(auto &testingCircle : newCircles){
    // New track segment cannot have greater radius (within some tolerance)
    if(testingCircle->GetRadius() > 1.1*previousCircle->GetRadius()) continue; // FILTER
    
    // The center of the new circle should be withing the previous circle
    double centerDifference = pointsProcessor->distanceXY(previousCircle->GetCenter(), testingCircle->GetCenter());
    if(centerDifference > 1.1*previousCircle->GetRadius()) continue; // FILTER
  
    // Here we calculate an angle between radius of the previous circle and radius of the testing circle
    // looking from the last point of previous circle. For perfectly matching circles that would be zero
    double r2_x = testingCircle->GetCenter()->GetX() - p_x;
    double r2_y = testingCircle->GetCenter()->GetY() - p_y;
    double r2_mod = sqrt(r2_x*r2_x + r2_y*r2_y);
    double alpha = acos((r1_x*r2_x + r1_y*r2_y) / (r1_mod*r2_mod));
    radiiAnglesHist->Fill(alpha);
    
    if(alpha > 0.05) continue; // FILTER
    
    // This circle is better than previous if it's radius is closer to the desired one
    // TODO: This condition may not be the best one...
    double radiiDifference = previousCircle->GetRadius() - testingCircle->GetRadius();
    
    if(radiiDifference < bestRadiiDifference){
      bestRadiiDifference = radiiDifference;
      bestCircle = make_unique<Circle>(testingCircle);
    }
  }
  
  return bestCircle;
}

unique_ptr<Helix> Fitter::FitHelix(const vector<shared_ptr<Point>> &_points,
                                   const shared_ptr<Track> &_track,
                                   const unique_ptr<Point> &_vertex)
{
  points = _points;
  track = _track;
  vertex = make_unique<Point>(_vertex);
  
  cout<<"================================================================"<<endl;
  cout<<"Fitting starts\n\n"<<endl;
  cout<<"Initial number of points:"<<points.size()<<endl;
  
  double minPointsSeparation = 3.0;
  double chi2threshold = 1E-5;
  
  cout<<"min points separation:"<<minPointsSeparation<<endl;
  cout<<"chi2 threshold:"<<chi2threshold<<endl;
  
  //----------------------------------------------------------------------------------------
  
  vector<shared_ptr<Point>> pointsInCycle = pointsProcessor->FilterNearbyPoints(points, minPointsSeparation);
  cout<<"Points after cleanup:"<<pointsInCycle.size()<<endl;
  
  TH2D *pointsHist = new TH2D("points","points",
                              300, -layerR[nLayers-1], layerR[nLayers-1],
                              300, -layerR[nLayers-1], layerR[nLayers-1]);
//                              300, -2000, 2000,
//                              300, -2000, 2000);
  
  pointsHist->Fill(0.0, 0.0, 5.0);
  
  for(auto &p : pointsInCycle){
    pointsHist->Fill(p->GetX(),p->GetY());
  }
  pointsHist->Draw("colz");
  

  auto pointTriplets = BuildPointTriplets(pointsInCycle);
  cout<<"N valid point triplets:"<<pointTriplets.size()<<endl;
  
  // add circles for all momentum directions
  vector<unique_ptr<Circle>> circles;
  vector<unique_ptr<Circle>> circlesTmp;
  
  circlesTmp = FitCirclesAndAdjustFirstPoint(pointTriplets, 1, 1, chi2threshold);
  circles.insert(circles.end(),
                 make_move_iterator(circlesTmp.begin()),
                 make_move_iterator(circlesTmp.end()));
  
  circlesTmp = FitCirclesAndAdjustFirstPoint(pointTriplets,-1, 1, chi2threshold);
  circles.insert(circles.end(),
                 make_move_iterator(circlesTmp.begin()),
                 make_move_iterator(circlesTmp.end()));
  
  circlesTmp = FitCirclesAndAdjustFirstPoint(pointTriplets, 1,-1, chi2threshold);
  circles.insert(circles.end(),
                 make_move_iterator(circlesTmp.begin()),
                 make_move_iterator(circlesTmp.end()));
  
  circlesTmp = FitCirclesAndAdjustFirstPoint(pointTriplets,-1,-1, chi2threshold);
  circles.insert(circles.end(),
                 make_move_iterator(circlesTmp.begin()),
                 make_move_iterator(circlesTmp.end()));
  
  
  cout<<"N valid circles:"<<circles.size()<<endl;
  
  vector<unique_ptr<ArcSet2D>> potentialPionTracks = BuildArcSetsFromCircles(circles, pointTriplets);
  cout<<"N track seeds:"<<potentialPionTracks.size()<<endl;
  
  // fit more segments staring from seeds
  /*
  for(auto &pionTrack : potentialPionTracks){
    int iter=0;
    TGraph *pionTrackPoints = new TGraph();
    pionTrackPoints->SetMarkerStyle(25);
    pionTrackPoints->SetMarkerSize(1.0);
    pionTrackPoints->SetMarkerColor(kRed);
    
    while(1){
      // Get potential triplets of new points
      auto newPointTriplets = BuildPointTriplets(pionTrack, pointsInCycle);
      
      // Build a circle for each possible combination of points
      vector<unique_ptr<Circle>> newCircles = GetCirclesForPoints(newPointTriplets, chi2threshold);
      
      unique_ptr<Circle> bestCircle = GetBestCircle(newCircles, pionTrack->GetLastCircle());
      
      if(!bestCircle) break;
      
      pionTrackPoints->SetPoint(iter++,
                                bestCircle->GetPoints()[2]->GetX(),
                                bestCircle->GetPoints()[2]->GetY());
      
      // Add best circle's arc and the last point (the new one) to the pion's track candidate
      pionTrack->AddCircle(bestCircle, GetPhiRange(bestCircle, pionTrack));
      
      pionTrack->AddPoint(bestCircle->GetPoints()[2]);
    }
    
    pionTrackPoints->Draw("P");
  }
  */
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
  
  TGraph *graphDecay = new TGraph();
  graphDecay->SetPoint(0, track->GetDecayPoint()->GetX(), track->GetDecayPoint()->GetY());
  graphDecay->SetMarkerStyle(20);
  graphDecay->SetMarkerSize(1.0);
  graphDecay->SetMarkerColor(kRed);
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
  
  c1->cd(2);
  radiiAnglesHist->Draw();
  
  c1->Update();
  //---------------------
  
  unique_ptr<ArcSet2D> bestPionTrack = nullptr;
  unsigned long maxNarcs = 0;
  
  for(auto &pionTrack : potentialPionTracks){
    if(pionTrack->GetNarcs() > maxNarcs){
      maxNarcs = pionTrack->GetNarcs();
      bestPionTrack = make_unique<ArcSet2D>(pionTrack);
    }
  }
  if(!bestPionTrack) return nullptr;
  
  auto helix = make_unique<Helix>(make_unique<Point>(bestPionTrack->GetOrigin()),
                                  make_unique<Point>(bestPionTrack->GetCircle(0)->GetMomentum()),
                                  track->GetCharge());
  return helix;
}
