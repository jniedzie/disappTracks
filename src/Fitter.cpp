//
//  Fitter.cpp
//
//  Created by Jeremi Niedziela on 14/12/2018.
//

#include "Fitter.hpp"

Fitter::Fitter() :
pointsProcessor(make_unique<PointsProcessor>()),
helixProcessor(make_unique<HelixProcessor>())
{
  c1 = new TCanvas("c1","c1",1000,1000);
  c1->cd();
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
  
  double minPx = config->minPx;
  double minPy = config->minPy;
  double maxPx = config->maxPx;
  double maxPy = config->maxPy;
  
  double minL = layerR[track->GetLastBarrelLayer()];
  double maxL = layerR[track->GetLastBarrelLayer()+1];
  
  double trackTheta = track->GetTheta();
  double trackPhi = track->GetPhi();
  
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
                       pxSign*(maxPx-minPx)/2.,
                       pxSign > 0 ? minPx : -maxPx,
                       pxSign > 0 ? maxPx : -minPx);
  
  SetParameter(fitter, 2, "py",
                       pySign*(maxPy-minPy)/2.,
                       pySign > 0 ? minPy : -maxPy,
                       pySign > 0 ? maxPy : -minPy);
  
  // Store fitted circles for each triplet of points
  vector<unique_ptr<Circle>> circles;
  
  int nPoints = (int)points2D.size();
  double circleThickness = config->circleThickness;
  for(int i=0;i<nPoints;i++){
    for(int j=i+1;j<nPoints;j++){
      for(int k=j+1;k<nPoints;k++){
        
        auto chi2Function = [&](const double *par) {
          double f = 0;
          
          double L = par[0];
          double px = par[1];
          double py = par[2];
          
          double x0 = L*cos(trackPhi) + 10*vertex->GetX();
          double y0 = L*sin(trackPhi) + 10*vertex->GetY();
          double z0 = L/sin(trackTheta)*cos(trackTheta) + 10*vertex->GetZ();
          
          unique_ptr<Point> decayPoint  = make_unique<Point>(x0,y0,z0);
          unique_ptr<Point> momentum    = make_unique<Point>(px,py,0);
          Circle circle(decayPoint, momentum);
          
          f  = pow(circle.GetDistanceToPoint(points2D[i]),2);
          f += pow(circle.GetDistanceToPoint(points2D[j]),2);
          f += pow(circle.GetDistanceToPoint(points2D[k]),2);
          
          return f;
        };
        fitFunction = ROOT::Math::Functor(chi2Function, nPar);
        double pStart[nPar];
        fitter->SetFCN(fitFunction, pStart);
        
        if(fitter->FitFCN()) {
          auto result = fitter->Result();
          
          double L = result.GetParams()[0];
          double px = result.GetParams()[1];
          double py = result.GetParams()[2];
          
          double x0 = L*cos(trackPhi) + 10*vertex->GetX();
          double y0 = L*sin(trackPhi) + 10*vertex->GetY();
          double z0 = L/sin(trackTheta)*cos(trackTheta) + 10*vertex->GetZ();
          
          unique_ptr<Point> decayPoint = make_unique<Point>(x0,y0,z0);
          unique_ptr<Point> momentum = make_unique<Point>(px,py,0);
          unique_ptr<Circle> circle = make_unique<Circle>(decayPoint, momentum);
          
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
  Circle::RemoveSimilarCircles(circles);
  return circles;
}

unique_ptr<Helix> Fitter::GetBestFittingHelix(shared_ptr<vector<Point>> _points,
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
    
    for(auto p : *points){
      pointsHist->Fill(p.GetX(),p.GetY());
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
    vector<Point> points = circle->GetPoints();
    
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
  Circle::RemoveSimilarCircles(circles);
  
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

vector<shared_ptr<Point>> Fitter::GetPointsInCycle(double cycleMaxZ, double minPointsSeparation)
{
  double trackZ = track->GetDecayPoint()->GetZ();
  
  vector<shared_ptr<Point>> pointsInCycle;
  
  int pointsIter=-1;
  
  // remove points that are too far in Z to be in the first helix cycle
  for(auto point : *points){
    if(fabs(point.GetZ() - trackZ) > cycleMaxZ) continue;
    bool tooClose = false;
    
    for(auto otherPoint : pointsInCycle){
      if(pointsProcessor->distance(point, *otherPoint) < minPointsSeparation) tooClose = true;
    }
    if(tooClose) continue;
    
    pointsIter++;
    
//    if(pointsIter != 0 &&
//       pointsIter != 1) continue;
    
    pointsInCycle.push_back(make_shared<Point>(point));
  }
  
  return pointsInCycle;
}

ROOT::Fit::Fitter* Fitter::GetCirclesFitter()
{
  int pxSign = -1;
  int pySign = 1;
  
  auto pxRange = range<double>(0, 1000);
  auto pyRange = range<double>(0, 1000);
  
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

vector<vector<shared_ptr<Point>>> Fitter::BuildPointTriplets(const shared_ptr<Point> point1,
                                                             const shared_ptr<Point> point2,
                                                             const vector<shared_ptr<Point>> &inputPoints,
                                                             const unique_ptr<Circle> &circle)
{
  vector<vector<shared_ptr<Point>>> pointTriplets;
  
  double stripSensorHalfLength = 100;
  
  for(auto point : inputPoints){
    bool isValidPoint = true;
    
    // make sure that it's not the same point as one of those already in the triplet
    if(point == point1) continue;
    if(point == point2) continue;
    
    // new point must be on the correct side of the previous arc in Z direction
    if(point1->GetZ()+stripSensorHalfLength < point2->GetZ()-stripSensorHalfLength  &&
       point->GetZ() < (point2->GetZ()-stripSensorHalfLength)) isValidPoint = false;
    if(point1->GetZ()-stripSensorHalfLength > point2->GetZ()+stripSensorHalfLength  &&
       point->GetZ() > (point2->GetZ()+stripSensorHalfLength)) isValidPoint = false;
    
    double pointR = sqrt(pow(point->GetX() - circle->GetCenter()->GetX(),2) +
                         pow(point->GetY() - circle->GetCenter()->GetY(),2));
    
    // it also has to be within the radius of the helix
    if(pointR > 1.1*circle->GetRadius()) isValidPoint = false;
    
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
                                                                 double chi2threshold)
{
  int nPar=3;
  auto fitter = GetCirclesFitter();
  
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
      
      if(result.MinFcnValue() > chi2threshold){
        cout<<"Trashing solution due to chi2 above threshold"<<endl;
        continue;
      }
      
      auto circle = GetCircleFromFitterParams(result.GetParams());
      
      p[0]->SetX(circle->GetDecayPoint()->GetX());
      p[0]->SetY(circle->GetDecayPoint()->GetY());
      p[0]->SetZ(circle->GetDecayPoint()->GetZ());
      
      circles.push_back(move(circle));
    }
  }
  
  return circles;
}

vector<unique_ptr<Circle>> Fitter::GetCirclesForPoints(vector<vector<shared_ptr<Point>>> &pointTriplets,
                                                       double chi2threshold)
{
  vector<unique_ptr<Circle>> circles;
  
  for(auto &p : pointTriplets){
    double x1 = p[0]->GetX();
    double y1 = p[0]->GetY();
    double x2 = p[1]->GetX();
    double y2 = p[1]->GetY();
    double x3 = p[2]->GetX();
    double y3 = p[2]->GetY();
    
    double x1_2 = x1*x1;
    double y1_2 = y1*y1;
    double x2_2 = x2*x2;
    double y2_2 = y2*y2;
    double x3_2 = x3*x3;
    double y3_2 = y3*y3;
    
    double x0 = ((x1_2+y1_2)*(y2-y3)+(x2_2+y2_2)*(y3-y1)+(x3_2+y3_2)*(y1-y2))/(2*(x1*(y2-y3)-y1*(x2-x3)+x2*y3-x3*y2));
    double y0 = ((x1_2+y1_2)*(x3-x2)+(x2_2+y2_2)*(x1-x3)+(x3_2+y3_2)*(x2-x1))/(2*(x1*(y2-y3)-y1*(x2-x3)+x2*y3-x3*y2));
    double R  = sqrt(pow(x0 - x1, 2) + pow(y0 - y1, 2));
    
    auto center     = make_unique<Point>(x0, y0, 0.0);
    auto decayPoint = make_unique<Point>(x1, y1, p[0]->GetZ());
    
    auto circle = make_unique<Circle>(decayPoint, center, R);
      
    circles.push_back(move(circle));
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
    return false;
  }
  
  double stripSensorLength = 200; // mm (to be determined more precisely later)
  
  double z0 = pointTriplet[0]->GetZ();
  double z1 = pointTriplet[1]->GetZ();
  double z2 = pointTriplet[2]->GetZ();
  
  // this is the range of where point 1 can be located along Z axis
  double z1min = z1 - stripSensorLength/2.;
  double z1max = z1 + stripSensorLength/2.;
  
  // which determines range in the slope:
  double slopeMin = (z1min - z0)/phi1;
  double slopeMax = (z1max - z0)/phi1;
  
  // point 2 can be located somewhere between:
  double z2min = z2 - stripSensorLength/2.;
  double z2max = z2 + stripSensorLength/2.;
  
  // check if point 2 lays between limits derived from positions of points 0 and 1
  double z2a = z0 + slopeMin * phi2;
  double z2b = z0 + slopeMax * phi2;
  
  if((z2a < z2min || z2a > z2max) &&
     (z2b < z2min || z2a > z2max)){
    return false;
  }
  
  return true;
}

range<double> Fitter::GetPhiRange(const unique_ptr<Circle> &circle, vector<shared_ptr<Point>> pointTriplet)
{
  double phiVertex = circle->GetPointAngle(pointTriplet[0]->GetX(), pointTriplet[0]->GetY());
  double phi1      = circle->GetPointAngle(pointTriplet[1]->GetX(), pointTriplet[1]->GetY());
  double phi2      = circle->GetPointAngle(pointTriplet[2]->GetX(), pointTriplet[2]->GetY());
  
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
    
    auto arcSet2D = make_unique<ArcSet2D>();
    arcSet2D->AddCircle(circle, GetPhiRange(circle, pointTriplets[iter]));
    arcSet2D->AddPoints(pointTriplets[iter]);
    arcs.push_back(move(arcSet2D));
    
    cout<<"Creating a seed:";
    pointTriplets[iter][0]->Print();cout<<"\t";
    pointTriplets[iter][1]->Print();cout<<"\t";
    pointTriplets[iter][2]->Print();cout<<"\n";
    circle->Print();
    
    iter++;
  }
  
  return arcs;
}

unique_ptr<Helix> Fitter::FitHelix(shared_ptr<vector<Point>> _points,
                                              const shared_ptr<Track> _track,
                                              const unique_ptr<Point> &_vertex)
{
  points = _points;
  track = _track;
  vertex = make_unique<Point>(_vertex);
  
  cout<<"================================================================"<<endl;
  cout<<"Fitting starts\n\n"<<endl;
  cout<<"Initial number of points:"<<points->size()<<endl;
  
  double cycleMaxZ = inf; // mm
  double minPointsSeparation = 3.0;
  double chi2threshold = 1;//1E-2;
  
  cout<<"max distance in Z for the first cycle:"<<cycleMaxZ<<endl;
  cout<<"min points separation:"<<minPointsSeparation<<endl;
  cout<<"chi2 threshold:"<<chi2threshold<<endl;
  
  //----------------------------------------------------------------------------------------
  
  vector<shared_ptr<Point>> pointsInCycle = GetPointsInCycle(cycleMaxZ, minPointsSeparation);
  auto pointTriplets = BuildPointTriplets(pointsInCycle);
  
  // Draw 2D histogram
  
  TH2D *pointsHist = new TH2D("points","points",
                              300, -layerR[nLayers-1], layerR[nLayers-1],
                              300, -layerR[nLayers-1], layerR[nLayers-1]);
  //                              300, -500, 300,
  //                              300, -300, 500);
  
  pointsHist->Fill(0.0,0.0,5.0);
  
  for(auto &p : pointsInCycle){
    pointsHist->Fill(p->GetX(),p->GetY());
  }
  pointsHist->Draw("colz");
  
  
  cout<<"N points in first cycle, that are not too close to each other:"<<pointsInCycle.size()<<endl;
  
  vector<unique_ptr<Circle>> circles = FitCirclesAndAdjustFirstPoint(pointTriplets, chi2threshold);
  
  vector<unique_ptr<ArcSet2D>> potentialPionTracks = BuildArcSetsFromCircles(circles, pointTriplets);
  
  potentialPionTracks.erase(potentialPionTracks.end()-1);
  
  // make another iteration
  
  for(auto &pionTrack : potentialPionTracks){
    unique_ptr<Circle> previousCircle = pionTrack->GetCircle(0);
    
    TGraph *prevPoints = new TGraph();
    prevPoints->SetPoint(0, pionTrack->GetPoint(1)->GetX(), pionTrack->GetPoint(1)->GetY());
    prevPoints->SetPoint(1, pionTrack->GetPoint(2)->GetX(), pionTrack->GetPoint(2)->GetY());
    prevPoints->SetMarkerStyle(25);
    prevPoints->SetMarkerSize(1.2);
    prevPoints->SetMarkerColor(kBlue);
    prevPoints->Draw("P");
    
    auto newPointTriplets = BuildPointTriplets(pionTrack->GetPoint(1),
                                               pionTrack->GetPoint(2),
                                               pointsInCycle, previousCircle);
    
    vector<unique_ptr<Circle>> newCircles = GetCirclesForPoints(newPointTriplets, chi2threshold);
    
    double bestRadiiDifference = inf;
    unique_ptr<Circle> bestCircle = nullptr;
    vector<shared_ptr<Point>> bestTriplet;
    
    for(uint i=0;i<newCircles.size();i++){
      auto testingCircle = make_unique<Circle>(newCircles[i]);
      auto testingPoints = newPointTriplets[i];

      // New track segment cannot have greater radius (within some tolerance)
      if(testingCircle->GetRadius() > 1.1*previousCircle->GetRadius()) continue;

      double radiiDifference = previousCircle->GetRadius() - testingCircle->GetRadius();

      if(radiiDifference < bestRadiiDifference){
        bestRadiiDifference = radiiDifference;
        bestCircle = make_unique<Circle>(testingCircle);
        bestTriplet = testingPoints;
      }
    }
    
    if(!bestCircle){
      // here I should remove this arcSet
      continue;
    }
    
    TGraph *newPoints = new TGraph();
    newPoints->SetPoint(0, bestTriplet[0]->GetX(), bestTriplet[0]->GetY());
    newPoints->SetPoint(1, bestTriplet[1]->GetX(), bestTriplet[1]->GetY());
    newPoints->SetPoint(2, bestTriplet[2]->GetX(), bestTriplet[2]->GetY());
    newPoints->SetMarkerStyle(25);
    newPoints->SetMarkerSize(1.0);
    newPoints->SetMarkerColor(kRed);
    newPoints->Draw("P");
    
    pionTrack->AddCircle(bestCircle, GetPhiRange(bestCircle, bestTriplet));
    pionTrack->AddPoint(bestTriplet[2]);
  }
  
   
  

  
  TGraph *graphDecay = new TGraph();
  graphDecay->SetPoint(1, track->GetDecayPoint()->GetX(), track->GetDecayPoint()->GetY());
  graphDecay->SetMarkerStyle(20);
  graphDecay->SetMarkerSize(1.0);
  graphDecay->SetMarkerColor(kRed);
  graphDecay->Draw("P");
  
  
  for(auto &pionTrack : potentialPionTracks){
    pionTrack->Print();
    
    TGraph *graph = new TGraph();
    graph->SetPoint(0, pionTrack->GetOrigin()->GetX(), pionTrack->GetOrigin()->GetY());
    graph->SetMarkerStyle(25);
    graph->SetMarkerSize(1.0);
    graph->SetMarkerColor(kGreen);
    graph->Draw("P");
    bool first = true;
    for(auto singleArc : pionTrack->GetArcs()){
      if(first){
        first = false;
//        continue;
      }
      
      singleArc->SetFillColorAlpha(kWhite, 0.0);
      singleArc->SetLineWidth(1.0);
      singleArc->SetLineColor(kRed);
      singleArc->Draw("sameLonly");
    }
  }
  
  c1->Update();
  //---------------------
  
  auto momentum = circles[0]->GetMomentum();
  auto helix = make_unique<Helix>(track->GetDecayPoint(), momentum, track->GetCharge());
  return helix;
}
