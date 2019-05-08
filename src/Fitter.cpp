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
  
  double maxChi2 = 1.0;
  double deltaPhiMax = TMath::Pi()/4.;
  
  vector<vector<shared_ptr<Point>>> pointsByLayer = pointsProcessor.SortByLayer(points);
  
  // find possible middle and last seeds' points
  int trackLayers = track.GetNtrackerLayers();
  vector<shared_ptr<Point>> possibleMiddlePoints = pointsByLayer[trackLayers];
  vector<shared_ptr<Point>> possibleLastPoints   = pointsByLayer[trackLayers+1];
  
  // build seeds for all possible combinations of points
  vector<Helix> fittedHelices;
  
  for(auto &middlePoint : possibleMiddlePoints){
    for(auto &lastPoint : possibleLastPoints){
      vector<shared_ptr<Point>> points = { middlePoint, lastPoint };
      unique_ptr<Helix> helix = FitSeed(points);
      if(!helix) continue;
      if(helix->chi2 > maxChi2) continue;
      
      helix->increasing = true; // add decreasing later
      fittedHelices.push_back(*helix);
    }
  }
  
  ExtendSeeds(fittedHelices, pointsByLayer, maxChi2, deltaPhiMax);
  MergeHelices(fittedHelices);
  
  vector<Helix> longHelices;
  
  // Remove helices that even after merging have only 3 hits
  for(int iHelix=0; iHelix<fittedHelices.size(); iHelix++){
    if(fittedHelices[iHelix].GetPoints().size() > 4){
      longHelices.push_back(fittedHelices[iHelix]);
    }
  }
  
  for(auto &helix : longHelices){
    if(helix.shouldRefit) RefitHelix(helix);
  }
  
  return longHelices;
}

unique_ptr<Helix> Fitter::FitSeed(const vector<shared_ptr<Point>> &points)
{
  double Lmin = layerR[track.GetNtrackerLayers()-1];
  double Lmax = layerR[track.GetNtrackerLayers()];
  
  auto fitter = GetSeedFitter(range<double>(Lmin, Lmax));
  
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
    Point vertex = pointsProcessor.GetPointOnTrack(L, track, eventVertex);
    double t = atan2(vertex.GetY() - y0, vertex.GetX() - x0);
    
    // Point on helix for t_vertex
    Point vertexClosest(x0 + (R0 - a*t)*cos(t),
                        y0 + (R0 - a*t)*sin(t),
                        z0 + (s0 - b*t)*t);
    
    double distX = pow(vertexClosest.GetX()-vertex.GetX(), 2);
    double distY = pow(vertexClosest.GetY()-vertex.GetY(), 2);
    double distZ = pow(vertexClosest.GetZ()-vertex.GetZ(), 2);
    
    double f = distX + distY + distZ;
    
    // Then add distances to all other points
    for(auto &p : points){
      t = atan2(p->GetY() - y0, p->GetX() - x0);
      
      // find helix point for this point's t
      double x = x0 + (R0 - a*t)*cos(t);
      double y = y0 + (R0 - a*t)*sin(t);
      double z = z0 + (s0 - b*t)*t;
      
      // calculate distance between helix and point's boundary (taking into account its errors)
      distX = pow(x-p->GetX(), 2);
      distY = pow(y-p->GetY(), 2);
      distZ = pow(z-p->GetZ(), 2);
      
      distX /= p->GetXerr() > 0 ? pow(p->GetXerr(), 2) : 1;
      distY /= p->GetYerr() > 0 ? pow(p->GetYerr(), 2) : 1;
      distZ /= p->GetZerr() > 0 ? pow(p->GetZerr(), 2) : 1;
      
      f += distX + distY + distZ;
    }
    return f;
  };
  

  int nPar=8;
  auto fitFunction = ROOT::Math::Functor(chi2Function, nPar);
  double pStart[nPar];
  for(int i=0; i<nPar; i++) pStart[i] = fitter->Config().ParSettings(i).Value();
  fitter->SetFCN(fitFunction, pStart);
  
  unique_ptr<Helix> resultHelix = nullptr;
  
  if(fitter->FitFCN()) {
    auto result = fitter->Result();
    
    // Build helix from fitters output
    HelixParams resultParams;
    resultParams.R0 = result.GetParams()[0];
    resultParams.a  = result.GetParams()[1];
    resultParams.s0 = result.GetParams()[2];
    resultParams.b  = result.GetParams()[3];
    
    double L  = result.GetParams()[4];
    double x0 = result.GetParams()[5];
    double y0 = result.GetParams()[6];
    double z0 = result.GetParams()[7];
    
    Point origin(x0, y0, z0);
    Point vertex = pointsProcessor.GetPointOnTrack(L, track, eventVertex);
    vertex.SetT(atan2(vertex.GetY() - origin.GetY(), vertex.GetX() - origin.GetX()));
    
    for(auto &p : points) p->SetT(atan2(p->GetY() - y0, p->GetX() - x0));
    
    // Check if ordering of the points is correct (third point in a cone relative to vertex+first point)
    double phi = pointsProcessor.GetPointingAngle(vertex, *points[0], *points[1]);
    
    if(phi < TMath::Pi()/2.){
      resultHelix = make_unique<Helix>(resultParams, vertex, origin, points, track);
      resultHelix->chi2 = result.MinFcnValue();
    }
  }
  
  return resultHelix;
}

void Fitter::ExtendSeeds(vector<Helix> &helices,
                         const vector<vector<shared_ptr<Point>>> &pointsByLayer,
                         double maxChi2,
                         double deltaPhiMax)
{
  bool finished;
  
  do{
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
        int lastPointLayer = helix.GetLastPoint()->GetLayer();
        vector<shared_ptr<Point>> possiblePoints;
        
        if(helix.increasing)  possiblePoints = pointsByLayer[lastPointLayer+1];
        else                  possiblePoints = pointsByLayer[lastPointLayer-1];
        
        // try to extend by all possible points
        for(auto &point : possiblePoints){
          
          // Check if new point is within some cone
          size_t nHelixPoints = helix.GetPoints().size();
          double phi = pointsProcessor.GetPointingAngle(*helix.GetPoints()[nHelixPoints-2],
                                                        *helix.GetPoints()[nHelixPoints-1],
                                                        *point);
          if(phi > deltaPhiMax) continue;
          
          /// Extend helix by the new point and refit its params
          Helix helixCopy(helix);
          helixCopy.AddPoint(point);
          RefitHelix(helixCopy);
          
          // check if chi2 is small enough
          if(helixCopy.chi2 > maxChi2) continue;
          
          extendedHelices.push_back(helixCopy);
        }
        // if it was possible to extend the helix
        if(extendedHelices.size() != 0){
          helicesAfterExtending.insert(helicesAfterExtending.end(),
                                       extendedHelices.begin(),
                                       extendedHelices.end());
          finished = false;
        }
        else{ // if helix could not be extended
          helix.isFinished = true;
          helicesAfterExtending.push_back(helix); // is this needed?
        }
      }
    }
    helices.clear();
    for(auto &h : helicesAfterExtending) helices.push_back(h);
  }
  while(!finished);
}

void Fitter::MergeHelices(vector<Helix> &helices)
{
  // Merge helices that are very similar to each other
  bool merged=true;
  
  while(merged){
    merged=false;
    
    for(int iHelix1=0; iHelix1<helices.size(); iHelix1++){
      Helix &helix1 = helices[iHelix1];
      
      for(int iHelix2=iHelix1+1; iHelix2<helices.size(); iHelix2++){
        Helix &helix2 = helices[iHelix2];
        
        vector<shared_ptr<Point>> points1 = helix1.GetPoints();
        vector<shared_ptr<Point>> points2 = helix2.GetPoints();
        
        sort(points1.begin(), points1.end());
        sort(points2.begin(), points2.end());
        
        vector<shared_ptr<Point>> samePoints;
        set_intersection(points1.begin(), points1.end(),
                         points2.begin(), points2.end(),
                         back_inserter(samePoints));
        
        double samePointsFraction = samePoints.size()/(double)points1.size();
        
//        cout<<"Helix "<<helix1.uniqueID<<"\tand "<<helix2.uniqueID<<" share "<<samePointsFraction<<" points"<<endl;
        
        // merging
        if(samePointsFraction > 0.4){
          // remove second helix
          helices.erase(helices.begin() + iHelix2);
          
          // update first helix
          unordered_set<shared_ptr<Point>> uniquePoints;
          for(auto &p : points1) uniquePoints.insert(p);
          for(auto &p : points2) uniquePoints.insert(p);
          vector<shared_ptr<Point>> allPoints(uniquePoints.begin(), uniquePoints.end());
          helix1.ReplacePoints(allPoints);
          helix1.shouldRefit = true;
          
          merged=true;
          break;
        }
      }
      
      if(merged) break;
    }
  }
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
    Point vertex = pointsProcessor.GetPointOnTrack(L, track, eventVertex);
    double t = atan2(vertex.GetY() - y0, vertex.GetX() - x0);
    
    // Point on helix for t_vertex
    Point vertexClosest(x0 + (R0 - a*t)*cos(t),
                        y0 + (R0 - a*t)*sin(t),
                        z0 + (s0 - b*t)*t);
    
    double distX = pow(vertexClosest.GetX()-vertex.GetX(), 2);
    double distY = pow(vertexClosest.GetY()-vertex.GetY(), 2);
    double distZ = pow(vertexClosest.GetZ()-vertex.GetZ(), 2);
    
    double f = distX + distY + distZ;
    
    // Then add distances to all other points
    bool first=true;
    for(auto &p : helixPoints){
      // Skip first point, as this is the old vertex that will be replaced
      if(first){ first = false; continue; }
      
      t = atan2(p->GetY() - y0, p->GetX() - x0);
      
      // find helix point for this point's t
      double x = x0 + (R0 - a*t)*cos(t);
      double y = y0 + (R0 - a*t)*sin(t);
      double z = z0 + (s0 - b*t)*t;
      
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
      if(fabs(z-p->GetZ()) > p->GetZerr()+ht){
        double distZ_1 = z - (p->GetZ() + p->GetZerr() + ht);
        double distZ_2 = z - (p->GetZ() - p->GetZerr() - ht);
        distZ = min(pow(distZ_1, 2), pow(distZ_2, 2));
      }
      
      distX /= p->GetXerr() > 0 ? pow(p->GetXerr(), 2) : fabs(p->GetX());
      distY /= p->GetYerr() > 0 ? pow(p->GetYerr(), 2) : fabs(p->GetY());
      distZ /= p->GetZerr() > 0 ? pow(p->GetZerr(), 2) : fabs(p->GetZ());
      
      f += distX + distY + distZ;
    }
    return f/(3*helixPoints.size()+8);
  };
  
  ROOT::Fit::Fitter *fitter = new ROOT::Fit::Fitter();
  int nPar = 8;
  ROOT::Math::Functor fitFunction = ROOT::Math::Functor(chi2Function, nPar);
  double pStart[nPar];
  fitter->SetFCN(fitFunction, pStart);
  
  double Lmin = layerR[track.GetNtrackerLayers()-1];
  double Lmax = layerR[track.GetNtrackerLayers()];
  
  SetParameter(fitter, 0, "R0", helix.helixParams.R0,
               GetRadiusInMagField(config.minPx, config.minPy, solenoidField),
               GetRadiusInMagField(config.maxPx, config.maxPy, solenoidField));
  SetParameter(fitter, 1, "a" , helix.helixParams.a, 0, 10000);
  SetParameter(fitter, 2, "s0", helix.helixParams.s0, -1000, 1000);
  SetParameter(fitter, 3, "b" , helix.helixParams.b, -10000, 0);
  SetParameter(fitter, 4, "L" , (helix.GetVertex()->GetX()-10*eventVertex.GetX())/cos(track.GetPhi()),
               Lmin, Lmax);
  SetParameter(fitter, 5, "x0", helix.GetOrigin().GetX(), -layerR[nLayers-1] , layerR[nLayers-1]);
  SetParameter(fitter, 6, "y0", helix.GetOrigin().GetY(), -layerR[nLayers-1] , layerR[nLayers-1]);
  SetParameter(fitter, 7, "z0", helix.GetOrigin().GetZ(), -pixelBarrelZsize  , pixelBarrelZsize);
  
  for(int i=0; i<nPar; i++) pStart[i] = fitter->Config().ParSettings(i).Value();
  
//  if(fitter->FitFCN()) {
  
  fitter->FitFCN();
    auto result = fitter->Result();
    
    double L  = result.GetParams()[4];
    Point vertex = pointsProcessor.GetPointOnTrack(L, track, eventVertex);
    
    helix.SetVertex(vertex);
    
    double x0 = result.GetParams()[5];
    double y0 = result.GetParams()[6];
    double z0 = result.GetParams()[7];
    Point origin(x0, y0, z0);
    helix.UpdateOrigin(origin);
    
    HelixParams resultParams;
    resultParams.R0 = result.GetParams()[0];
    resultParams.a  = result.GetParams()[1];
    resultParams.s0 = result.GetParams()[2];
    resultParams.b  = result.GetParams()[3];
    helix.helixParams = resultParams;
  
  helix.chi2 = result.MinFcnValue();
//  }
}

ROOT::Fit::Fitter* Fitter::GetSeedFitter(range<double> rangeL)
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
  
  double minR = GetRadiusInMagField(config.minPx, config.minPy, solenoidField);
  double maxR = GetRadiusInMagField(config.maxPx, config.maxPy, solenoidField);
  
  SetParameter(fitter, 0, "R0", (minR+maxR)/2., 0, 10000);
  SetParameter(fitter, 1, "a" , (minR+maxR)/2., 0, 10000);
  SetParameter(fitter, 2, "s0", 0, -1000, 1000);
  SetParameter(fitter, 3, "b" , 0, -10000, 0);
  SetParameter(fitter, 4, "L" , (rangeL.GetMin()+rangeL.GetMax())/2., rangeL.GetMin(), rangeL.GetMax());
  SetParameter(fitter, 5, "x0", 0, -layerR[nLayers-1] , layerR[nLayers-1]);
  SetParameter(fitter, 6, "y0", 0, -layerR[nLayers-1] , layerR[nLayers-1]);
  SetParameter(fitter, 7, "z0", 0, -pixelBarrelZsize  , pixelBarrelZsize);
  
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
