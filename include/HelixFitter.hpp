//
//  HelixFitter.hpp
//
//  Created by Jeremi Niedziela on 10/12/2018.
//

#ifndef HelixFitter_hpp
#define HelixFitter_hpp

#include "Helpers.hpp"

namespace HelixFitter {
  struct Point
  {
    Point(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
    double x,y,z;
    void Print(){cout<<"("<<x<<","<<y<<","<<z<<")"<<endl;}
    double distance(Point p){return sqrt(pow(x-p.x,2)+pow(y-p.y,2)+pow(z-p.z,2));}
    bool isPionHit = false;
  };
  
  struct Helix
  {
    Helix(double _R, double _c, double _x0, double _y0, double _z0) : R(_R), c(_c), x0(_x0), y0(_y0), z0(_z0) {}
    Helix(){}
    double R,c;
    double x0,y0,z0;
    void Print(){cout<<"R:"<<R<<"\tc:"<<c<<"\toffset:("<<x0<<","<<y0<<","<<z0<<")\tnPoints:"<<nPoints<<"\tnPionPoints:"<<nPionPoints<<endl;}
    int nPoints = 0;
    int nPionPoints = 0;
  };
  
  double helixMarkerSize = 0.1;
  double xAxisMin = -1;
  double xAxisMax =  1;
  double yAxisMin = -1;
  double yAxisMax =  1;
  double zAxisMin = 0;
  double zAxisMax = 10;
  
  double helixThickness = 2.0;
  
  vector<Point> points;
  TFitter *fitter;
  double minR, maxR, startR;
  double minC, maxC, startC;
  double minX, maxX, startX;
  double minY, maxY, startY;
  double minZ, maxZ, startZ;
  
  double minL, maxL, startL;
  
  double theta, phi;
  Helix bestSeedHelix;
  vector<Point> bestSeedPoints;
  
  void RunFitter()
  {
    double arglist[1] = {0};
    fitter->ExecuteCommand("MIGRAD", arglist, 0);
  }
  
  Helix GetFittedHelix()
  {
    double R = fitter->GetParameter(0);
    double c = fitter->GetParameter(1);
    
    double L = fitter->GetParameter(2);
    double x0 = L*sin(theta)*cos(phi);
    double y0 = L*sin(theta)*sin(phi);
    double z0 = L*cos(theta);
  
    return Helix(R,c,x0+R,y0,z0-c*TMath::Pi());
  }
  
  Point GetClosestPointOnHelix(Point p, Helix h)
  {
    double t = atan2(p.y-h.y0, p.x-h.x0);
    
    double x = h.R*cos(t) + h.x0;
    double y = h.R*sin(t) + h.y0;
    double z = h.c*t      + h.z0 + 2*h.c*TMath::Pi();
    
    double c = fabs(h.c);
    
    while(fabs(z-p.z) >= c){
      if(z < p.z) z += c;
      else        z -= c;
    }
    
    return Point(x,y,z);
  }
  
  void CountPionPointsOnHelix(Helix &helix)
  {
    for(Point p : points){
      Point q = GetClosestPointOnHelix(p,helix);
      
      double d = p.distance(q);
      
      if(d < helixThickness){
        helix.nPoints++;
        if(p.isPionHit) helix.nPionPoints++;
      }
    }
  }
  
  void fitFunction(Int_t&, Double_t *, Double_t &chi2, Double_t *par, Int_t)
  {
    double x0 = par[2]*sin(theta)*cos(phi);
    double y0 = par[2]*sin(theta)*sin(phi);
    double z0 = par[2]*cos(theta);
    
    Helix helix(par[0],par[1], x0+par[0], y0, z0-par[1]*TMath::Pi());
    
    int nPoints=0;
    
    double chi2local = 0;
    
    for(Point p : points){
      Point q = GetClosestPointOnHelix(p,helix);
      double d = p.distance(q);
      
      if(d < helixThickness){
        nPoints++;
        chi2local += d;
      }
    }
    if(nPoints > 5) chi2local /= nPoints;
    else            chi2local = 99999;
    
    chi2 = chi2local;
    
    cout<<"R:"<<par[0]<<"\tc:"<<par[1]<<"\tL:"<<par[2]<<"\tchi2/nPoints:"<<chi2<<"\n";
    
    return;
  }
  
  void fitFunctionWithSeeds(Int_t&, Double_t *, Double_t &chi2, Double_t *par, Int_t)
  {
    // par[0] is a secondary vertex position along the chargino's track
    double x0 = par[0]*sin(theta)*cos(phi);
    double y0 = par[0]*sin(theta)*sin(phi);
    double z0 = par[0]*cos(theta);
    chi2 = 0;
    
    Helix helix = Helix();
    
    double bestSeedNpoints = -1;
    double bestChi2 = 999999;
    
    vector<Point> helixPoints;
    
    for(int iSeed=0;iSeed<points.size();iSeed++){
      Point seedHit = points[iSeed];
      helixPoints.push_back(seedHit);
      // build a seed for this hit
      double seedT = atan2(seedHit.y,seedHit.x)*z0/seedHit.z;
      double R = x0/cos(seedT*z0/seedHit.z);
      double c = seedHit.z/atan2(seedHit.y,seedHit.x);
    
      // center of the helix has to be shifted:
      
      if(R < minR || R > maxR || c < minC || c > maxC){
        // then this is not a correct seed
        continue;
      }
      double chi2tmp = 9999999;
      bool first = true;
      
      // build helix from this seed (shifted approprietly)
      helix = Helix(R,c,x0+R,y0,z0-c*TMath::Pi());
      
      double nPoints=0;
      
      for(int iPoint=0;iPoint<points.size();iPoint++){ // count points along this seed helix
        if(iPoint==iSeed) continue;
        
        Point p = points[iPoint];
        Point closestOnHelix = GetClosestPointOnHelix(p, helix);
        
        if(p.distance(closestOnHelix) < helixThickness){
          nPoints++;
          helixPoints.push_back(p);
          if(first){
            chi2tmp = 0;
            first = false;
          }
          chi2tmp += pow(p.x - closestOnHelix.x, 2) +
                     pow(p.y - closestOnHelix.y, 2) +
                     pow(p.z - closestOnHelix.z, 2);
        }
      }
      if(nPoints > 0) chi2tmp /= nPoints;
//      chi2tmp += 100*pow(23-nPoints,2)/23.;
      
      if(chi2tmp < bestChi2){
        bestChi2 = chi2tmp;
        bestSeedHelix = helix;
        bestSeedHelix.nPoints = nPoints;
        
        int nPionPoints = 0;
        for(auto p : helixPoints){
          if(p.isPionHit) nPionPoints++;
        }
        bestSeedHelix.nPionPoints = nPionPoints;
        bestSeedPoints = helixPoints;
        
      }
    }
    chi2 = bestChi2;
    
    return;
  }
  
  void InitHelixFitter()
  {
    TVirtualFitter::SetDefaultFitter("Minuit");
    const int nPar = 3;
    fitter = new TFitter(nPar);
    fitter->SetFCN(fitFunction);
    
    fitter->SetParameter(0, "R", startR,  0.001, minR, maxR);
    fitter->SetParameter(1, "c", startC,  0.001, minC, maxC);
    fitter->SetParameter(2, "L", startL,  0.001, minL, maxL);
    
//    fitter->FixParameter(0);
//    fitter->FixParameter(1);
//    fitter->FixParameter(2);
    
    double args[1] = {1}; // put to 0 for results only, or to -1 for no garbage
    fitter->ExecuteCommand( "SET PRINTOUT"  , args, 1);
    //  fitter->ExecuteCommand( "SET NOWARNINGS", &args, 0);
    fitter->ExecuteCommand( "SET PRINT"     , args, 1);
    //  double fitterError[1] = {5.0};
    //  fitter->ExecuteCommand( "SET ERR", fitterError, 1);
    double strategyLevel[1] = {2};
      fitter->ExecuteCommand( "SET STR", strategyLevel, 1);
  }
  
  void InitHelixFitterWithSeeds()
  {
    TVirtualFitter::SetDefaultFitter("Minuit");
    const int nPar = 1;
    fitter = new TFitter(nPar);
    fitter->SetFCN(fitFunctionWithSeeds);
    
    fitter->SetParameter(0, "L", startL,  0.001, minL, maxL);
//    fitter->FixParameter(0);
    
    double args[1] = {1}; // put to 0 for results only, or to -1 for no garbage
    fitter->ExecuteCommand( "SET PRINTOUT"  , args, 1);
    //  fitter->ExecuteCommand( "SET NOWARNINGS", &args, 0);
    fitter->ExecuteCommand( "SET PRINT"     , args, 1);
    //  double fitterError[1] = {5.0};
    //  fitter->ExecuteCommand( "SET ERR", fitterError, 1);
//    double strategyLevel[1] = {2};
    //  fitter->ExecuteCommand( "SET STR", strategyLevel, 1);
  }
  
}



#endif /* HelixFitter_hpp */
