//
//  Fitter.cpp
//
//  Created by Jeremi Niedziela on 14/12/2018.
//

#include "Fitter.hpp"

Fitter::Fitter(shared_ptr<FitterConfig> _config) :
config(_config),
pointsProcessor(make_unique<PointsProcessor>())
{
  
}

Fitter::~Fitter()
{
  
}

vector<unique_ptr<Circle>> Fitter::FitCirclesToPoints(int pxSign, int pySign, int charge)
{
  // Prepare 2D projections in XY
  vector<Point> points2D;
  vector<vector<Point>> pointsByLine = pointsProcessor->SplitPointsIntoLines(points, config->GetLinesToleranceForCircles());
  
  for(vector<Point> line : pointsByLine){
    if(line.size() >= config->GetMinPointsAlongZ()){
      points2D.push_back(Point(line));
    }
  }
  
  double minPx = config->GetMinPx();
  double minPy = config->GetMinPy();
  double maxPx = config->GetMaxPx();
  double maxPy = config->GetMaxPy();
  double minL = config->GetMinL();
  double maxL = config->GetMaxL();
  
  // Create fitter to fit circles to 2D distribution
  ROOT::Fit::Fitter *fitter = new ROOT::Fit::Fitter();
  auto f = [&](const double *par) {return 0;};
  int nPar = 3;
  ROOT::Math::Functor fitFunction = ROOT::Math::Functor(f, nPar);
  double pStart[nPar];
  fitter->SetFCN(fitFunction, pStart);
  
  
  double helixThickness = config->GetHelixThickness();
  
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
  double circleThickness = config->GetCircleThickness();
  for(int i=0;i<nPoints;i++){
    for(int j=i+1;j<nPoints;j++){
      for(int k=j+1;k<nPoints;k++){
        
        auto chi2Function = [&](const double *par) {
          double f = 0;
          
          double L = par[0];
          double px = par[1];
          double py = par[2];
          
          double x0 = L*sin(trackTheta)*cos(trackPhi);
          double y0 = L*sin(trackTheta)*sin(trackPhi);
          double z0 = L*cos(trackTheta);
          
          unique_ptr<Point> decayPoint  = make_unique<Point>(x0,y0,z0);
          unique_ptr<Point> momentum    = make_unique<Point>(px,py,0);
          Circle circle(decayPoint, momentum, charge, config);
          
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
          
          double x0 = L*sin(trackTheta)*cos(trackPhi);
          double y0 = L*sin(trackTheta)*sin(trackPhi);
          double z0 = L*cos(trackTheta);
          
          unique_ptr<Point> decayPoint = make_unique<Point>(x0,y0,z0);
          unique_ptr<Point> momentum = make_unique<Point>(px,py,0);
          unique_ptr<Circle> circle = make_unique<Circle>(decayPoint, momentum, charge, config);
          
          int nPoints=0;
          for(Point p : points2D){
            if(circle->GetDistanceToPoint(p) < circleThickness) nPoints++;
          }
          if(nPoints > 3){
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

unique_ptr<Helix> Fitter::GetBestFittingHelix(vector<Point> _points, double _trackTheta, double _trackPhi,
                                              bool drawCircles)
{
  points = _points;
  trackTheta = _trackTheta;
  trackPhi = _trackPhi;
  
  vector<unique_ptr<Circle>> circles = GetAllCirclesForPoints();
  
  if(circles.size() == 0){
    cout<<"No circles were found"<<endl;
    return nullptr;
  }
  if(drawCircles){
    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    c1->cd();
    TH2D *pointsHist = new TH2D("points","points",
                            250, -250, 250,
                            250, -250, 250);
    
    for(auto p : points){
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
  
  double minPz = config->GetMinPz();
  double maxPz = config->GetMaxPz();
  
  for(auto &circle : circles){
    vector<Point> points = circle->GetPoints();
    
    for(double pz = maxPz; pz >= minPz ; pz-=config->GetStepPz()){
      unique_ptr<Helix> helix = make_unique<Helix>(circle, pz);
      helix->SetPoints(points);
      helix->CalculateNregularPoints();
      
      int nRegularPoints = helix->GetNregularPoints();
      double fractionRegularPoints = nRegularPoints/(double)helix->GetNpoints();
      
      // Here is a condition to accept new solution as the best one
      // Accept as a new best solution if:
      // - it gives more reqular points than before or,
      // - it gives the same number of regular points, but they counstitute higher fraction of all points than before
      if(nRegularPoints < maxNregularPoints) continue;
      else if(nRegularPoints == maxNregularPoints){
        if(fractionRegularPoints - maxFractionRegularPoints < 0.001) continue;
      }
      
      // If we reach till this point, save this solution as the best one so far
      bestHelix = move(helix);
      maxNregularPoints = nRegularPoints;
      maxFractionRegularPoints = fractionRegularPoints;
    }
    /*
     for(double pz = -maxPz; pz <= -minPz ; pz+=config->GetStepPz()){
     double c = Point(circle->GetMomentum()->GetX(), circle->GetMomentum()->GetY(), pz).GetVectorSlopeC();
     unique_ptr<Helix> helix = make_unique<Helix>(c, circle, config);
     helix->SetPoints(points);
     helix->CalculateNregularPoints();
     
     int nRegularPoints = helix->GetNregularPoints();
     double fractionRegularPoints = nRegularPoints/(double)helix->GetNpoints();
     
     // Here is a condition to accept new solution as the best one
     // Accept as a new best solution if:
     // - it gives more reqular points than before or,
     // - it gives the same number of regular points, but they counstitute higher fraction of all points than before
     if(nRegularPoints < maxNregularPoints) continue;
     else if(nRegularPoints == maxNregularPoints){
     if(fractionRegularPoints - maxFractionRegularPoints < 0.001) continue;
     }
     
     // If we reach till this point, save this solution as the best one so far
     helix->SetPz(pz);
     bestHelix = move(helix);
     maxNregularPoints = nRegularPoints;
     maxFractionRegularPoints = fractionRegularPoints;
     }
     */
  }
  
  return bestHelix;
}

vector<unique_ptr<Circle>> Fitter::GetAllCirclesForPoints()
{
  // Collect circles for positive charge
  int charge = 1;
  vector<unique_ptr<Circle>> circles, circlesNeg, circlesTmp;
  
  circles    = FitCirclesToPoints( 1,  1, charge);
  circlesTmp = FitCirclesToPoints(-1,  1, charge);
  circles.insert(circles.end(), make_move_iterator(circlesTmp.begin()), make_move_iterator(circlesTmp.end()));
  circlesTmp = FitCirclesToPoints( 1, -1, charge);
  circles.insert(circles.end(), make_move_iterator(circlesTmp.begin()), make_move_iterator(circlesTmp.end()));
  circlesTmp = FitCirclesToPoints(-1, -1, charge);
  circles.insert(circles.end(), make_move_iterator(circlesTmp.begin()), make_move_iterator(circlesTmp.end()));
  Circle::RemoveSimilarCircles(circles);
  
  charge = -1;
  
  // Collect circles for nevgative charge
  circlesNeg = FitCirclesToPoints( 1,  1, charge);
  circlesTmp = FitCirclesToPoints(-1,  1, charge);
  circlesNeg.insert(circlesNeg.end(), make_move_iterator(circlesTmp.begin()), make_move_iterator(circlesTmp.end()));
  circlesTmp = FitCirclesToPoints( 1, -1, charge);
  circlesNeg.insert(circlesNeg.end(), make_move_iterator(circlesTmp.begin()), make_move_iterator(circlesTmp.end()));
  circlesTmp = FitCirclesToPoints(-1, -1, charge);
  circlesNeg.insert(circlesNeg.end(), make_move_iterator(circlesTmp.begin()), make_move_iterator(circlesTmp.end()));
  Circle::RemoveSimilarCircles(circlesNeg);
  
  // Merge circles for positive and negative charges
  circles.insert(circles.end(), make_move_iterator(circlesNeg.begin()), make_move_iterator(circlesNeg.end()));
  
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
