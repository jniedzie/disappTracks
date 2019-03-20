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



unique_ptr<Helix> Fitter::FitHelix(shared_ptr<vector<Point>> _points,
                                              const shared_ptr<Track> _track,
                                              const unique_ptr<Point> &_vertex)
{
  points = _points;
  track = _track;
  vertex = make_unique<Point>(_vertex);
  
  double pz = -85;
  double trackZ = track->GetDecayPoint()->GetZ();
  double cycleMaxZ = 500; // mm
  double minPointsSeparation = 3.0;
  double chi2threshold = 1E-2;
  
  int pxSign = -1;
  int pySign = 1;
  
  //----------------------------------------------------------------------------------------
  
  vector<Point> pointsInCycle;
  
  int pointsIter=-1;
  
  // remove points that are too far in Z to be in the first helix cycle
  for(auto point : *points){
    pointsIter++;
    
    if(pointsIter != 10 &&
       pointsIter != 11) continue;
    
    if(fabs(point.GetZ() - trackZ) < cycleMaxZ){
      bool tooClose = false;
      
      for(auto otherPoint : pointsInCycle){
        if(pointsProcessor->distance(point, otherPoint) < minPointsSeparation) tooClose = true;
      }
      if(!tooClose) pointsInCycle.push_back(point);
    }
  }
  
//  vector<vector<Point>> pointsByLine = pointsProcessor->SplitPointsIntoLines(points, 10 /* mm */);
  
  auto pxRange = range<double>(50, 500);
  auto pyRange = range<double>(50, 500);
  
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
  
  
  
  int nPoints = (int)pointsInCycle.size();
  
  // Store fitted circles for each triplet of points
  vector<unique_ptr<Circle>> circles;
  vector<vector<Point>> pointTriplets;
  
  for(int i=0;i<nPoints;i++){
    for(int j=i+1;j<nPoints;j++){
      vector<Point> points = {Point(0,0,0), pointsInCycle[i], pointsInCycle[j] };
      pointTriplets.push_back(points);
    }
  }
  
  double trackTheta = track->GetTheta();
  double trackPhi = track->GetPhi();
  
  for(auto &points : pointTriplets){
    cout<<"Fitting points triplet:";
    points[0].Print();cout<<",\t";
    points[1].Print();cout<<",\t";
    points[2].Print();cout<<endl;
    
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
      
      f  = pow(circle.GetDistanceToPoint(decayPoint),2);
      f += pow(circle.GetDistanceToPoint(points[1]),2);
      f += pow(circle.GetDistanceToPoint(points[2]),2);
      
      return f;
    };
    
    fitFunction = ROOT::Math::Functor(chi2Function, nPar);
    double pStart[nPar];
    fitter->SetFCN(fitFunction, pStart);
    
    if(fitter->FitFCN()) {
      auto result = fitter->Result();
      
      if(result.MinFcnValue() > chi2threshold) continue;
      
      double L = result.GetParams()[0];
      double px = result.GetParams()[1];
      double py = result.GetParams()[2];
      
      double x0 = L*cos(trackPhi) + 10*vertex->GetX();
      double y0 = L*sin(trackPhi) + 10*vertex->GetY();
      double z0 = L/sin(trackTheta)*cos(trackTheta) + 10*vertex->GetZ();
      
      unique_ptr<Point> decayPoint = make_unique<Point>(x0,y0,z0);
      unique_ptr<Point> momentum = make_unique<Point>(px,py,0);
      unique_ptr<Circle> circle = make_unique<Circle>(decayPoint, momentum);
      
      points[0].SetX(decayPoint->GetX());
      points[0].SetY(decayPoint->GetY());
      points[0].SetZ(decayPoint->GetZ());
      
      circles.push_back(move(circle));
    }
    
    cout<<"Points after fitting:";
    points[0].Print();cout<<",\t";
    points[1].Print();cout<<",\t";
    points[2].Print();cout<<endl;
  }
  
  vector<unique_ptr<ArcSet2D>> arcs;
  int iter=0;
  
  for(auto &circle : circles){
    
    unique_ptr<Point> center = circle->GetCenter();
    
    double phiVertex = -TMath::Pi()/2. +atan2( (pointTriplets[iter][0].GetX()-center->GetX()),
                                              -(pointTriplets[iter][0].GetY()-center->GetY()));
    
    double phi1 = -TMath::Pi()/2. +atan2( (pointTriplets[iter][1].GetX()-center->GetX()),
                                         -(pointTriplets[iter][1].GetY()-center->GetY()));
    
    double phi2 = -TMath::Pi()/2. +atan2( (pointTriplets[iter][2].GetX()-center->GetX()),
                                         -(pointTriplets[iter][2].GetY()-center->GetY()));
    
    iter++;
    
    // Reject cases where hits are on the both sides of the vertex point instead of forming a tracklet
    if((phi1 < phiVertex && phi2 > phiVertex) ||
       (phi2 < phiVertex && phi1 > phiVertex)){
      continue;
    }
    
    double phiMin = min(min(phi1, phi2), phiVertex);
    double phiMax = max(max(phi1, phi2), phiVertex);
    
//    cout<<"Phi vertex:"<<phiVertex<<endl;
//    cout<<"Phi 1:"<<phi1<<endl;
//    cout<<"Phi 2:"<<phi2<<endl;
//    cout<<"Phi min:"<<phiMin<<endl;
//    cout<<"Phi max:"<<phiMax<<endl;
    
    auto arcSet2D = make_unique<ArcSet2D>();
    arcSet2D->AddCircle(circle, range<double>(phiMin/TMath::Pi() * 180,
                                              phiMax/TMath::Pi() * 180));
    arcs.push_back(move(arcSet2D));
  }
  
  
//  Circle::RemoveSimilarCircles(circles);
  
  // Draw 2D histogram
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->cd();
  TH2D *pointsHist = new TH2D("points","points",
                              500, -layerR[nLayers-1], layerR[nLayers-1],
                              500, -layerR[nLayers-1], layerR[nLayers-1]);
  
  for(auto p : pointsInCycle){
    pointsHist->Fill(p.GetX(),p.GetY());
  }
  for(auto points : pointTriplets){
    pointsHist->Fill(points[0].GetX(), points[0].GetY(), 5);
  }
    
  pointsHist->Draw("colz");
  
  for(auto &aa : arcs){
    for(auto a : aa->GetArcs()){
      a->SetFillColorAlpha(kWhite, 0.0);
      a->SetLineWidth(1.0);
      a->SetLineColor(kRed);
      a->Draw("sameLonly");
    }
  }
  c1->Update();
  
  auto momentum = make_unique<Point>(0,0,0);
  auto helix = make_unique<Helix>(track->GetDecayPoint(), momentum, track->GetCharge());
  return helix;
}
