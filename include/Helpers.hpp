//
//  Helpers.hpp
//
//  Created by Jeremi Niedziela on 13/06/2018.
//

#ifndef Helpers_h
#define Helpers_h

#include "Config.hpp"

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TLorentzVector.h>
#include <THStack.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TFitter.h>
#include <TMath.h>
#include <TGraph2D.h>
#include <TGraph.h>
#include <Fit/Fitter.h>
#include <Math/Functor.h>
#include <TEllipse.h>
#include <TArc.h>
#include <Fit/BinData.h>
#include <TSystem.h>
#include <TEveManager.h>
#include <TEveScene.h>
#include <TEvePointSet.h>
#include <TEveJetCone.h>
#include <TEveBox.h>
#include <TApplication.h>
#include <TGeoShape.h>
#include <TGeoTube.h>
#include <TEveGeoShape.h>
#include <TF3.h>
#include <TH3F.h>

#include <vector>
#include <iostream>
#include <map>
#include <numeric>
#include <variant>
#include <any>

using namespace std;

#define inf 99999999


// Plotting style
const double fillOpacity = 1.0;
const int fillStyleBack = 1000;
const int fillStyleSignal = 1000;
const int fillStyleData = 1000;

const vector<int> signalMarkers = {
  20, // wino m=300 cτ=3
  20, // wino m=300 cτ=10
  20, // wino m=300 cτ=30
  21, // wino m=500 cτ=10
  21, // wino m=500 cτ=20
  22, // wino m=650 cτ=10
  22, // wino m=650 cτ=20
  23, // wino m=800 cτ=10
  23, // wino m=800 cτ=20
  29, // wino m=1000 cτ=10
  29, // wino m=1000 cτ=20
};

const vector<vector<int>> backColors = {
  {230, 25 , 75 },  // QCD
  {60 , 180, 75 },  // Z->μμ + jets
  {0  , 130, 200},  // tops
  {245, 130, 48 },  // VV
  {145, 30 , 180},  // W->μν + jets
  {70 , 240, 240},  // Z->νν + jets
};

const vector<vector<int>> signalColors = {
  {128, 128, 0  },  // wino m=300 cτ=3
  {170, 110, 40 },  // wino m=300 cτ=10
  {128, 0  , 0  },  // wino m=300 cτ=30
  {170, 110, 40 },  // wino m=500 cτ=10
  {128, 0  , 0  },  // wino m=500 cτ=20
  {170, 110, 40 },  // wino m=650 cτ=10
  {128, 0  , 0  },  // wino m=650 cτ=20
  {170, 110, 40 },  // wino m=800 cτ=10
  {128, 0  , 0  },  // wino m=800 cτ=20
  {170, 110, 40 },  // wino m=1000 cτ=10
  {128, 0  , 0  },  // wino m=1000 cτ=20
};

// {140, 50 , 230} {170, 200, 195} {240, 50 , 100} {255, 225, 25 } {0  , 128, 128} {230, 190, 255} {170, 100, 195} {240, 50 , 230}

const vector<vector<int>> dataColors = {
  {200 , 10, 10},  // single electron (2017B)
};

// Names of background, signal and data samples
const vector<string> backgroundTitle = {
  "QCD",
  "Z#rightarrow#mu#mu + jets",
  "tt",
  "VV",
  "W#rightarrow#mu#nu + jets",
  "Z#rightarrow#nu#nu + jets",
};

const vector<string> signalTitle = {
  "Wino m=300 c#tau=3",
  "Wino m=300 c#tau=10",
  "Wino m=300 c#tau=30",
  "Wino m=500 c#tau=10",
  "Wino m=500 c#tau=20",
  "Wino m=650 c#tau=10",
  "Wino m=650 c#tau=20",
  "Wino m=800 c#tau=10",
  "Wino m=800 c#tau=20",
  "Wino m=1000 c#tau=10",
  "Wino m=1000 c#tau=20",
};

const vector<string> dataTitle = {
  "2017B",
};

// Path to trees with background, signal and data samples (also determines which samples will be merged)
const vector<vector<string>> inFileNameBackground = {
  // QCD
  {
    "../SR_MC/QCD_HT100to200/",
    "../SR_MC/QCD_HT200to300/",
    "../SR_MC/QCD_HT300to500/",
    "../SR_MC/QCD_HT500to700/",
    "../SR_MC/QCD_HT700to1000/",
    "../SR_MC/QCD_HT1000to1500/",
    "../SR_MC/QCD_HT1500to2000/",
    "../SR_MC/QCD_HT2000toInf/",
  },
  // DY + jets
  {
//    "../SR_MC/DYJetsM50_HT100to200/",
    "../SR_MC/DYJetsM50_HT100to200e/",
//    "../SR_MC/DYJetsM50_HT200to400/",
    "../SR_MC/DYJetsM50_HT200to400e/",
//    "../SR_MC/DYJetsM50_HT400to600/",
    "../SR_MC/DYJetsM50_HT400to600e/",
    "../SR_MC/DYJetsM50_HT600to800/",
    "../SR_MC/DYJetsM50_HT800to1200/",
    "../SR_MC/DYJetsM50_HT1200to2500/",
    "../SR_MC/DYJetsM50_HT2500toInf/",
  },
  // tops
  {
    "../SR_MC/TTHad/",
    "../SR_MC/TTLep/",
    "../SR_MC/TTSemi/",
    "../SR_MC/T_tch/",
    "../SR_MC/T_tWch/",
    "../SR_MC/TBar_tch/",
    "../SR_MC/TBar_tWch/",
  },
  // VV
  {
    "../SR_MC/WW/",
    "../SR_MC/WZ/",
    "../SR_MC/ZZ/",
  },
  // W->μν + jets
  {
    "../SR_MC/WJets_HT100to200/",
    "../SR_MC/WJets_HT200to400/",
    "../SR_MC/WJets_HT400to600/",
    "../SR_MC/WJets_HT600to800/",
    "../SR_MC/WJets_HT800to1200/",
    "../SR_MC/WJets_HT1200to2500/",
    "../SR_MC/WJets_HT2500toInf/",
  },
  // Z->νν + jets
  {
    "../SR_MC/ZvvJets_HT100to200/",
    "../SR_MC/ZvvJets_HT200to400/",
    "../SR_MC/ZvvJets_HT400to600/",
    "../SR_MC/ZvvJets_HT600to800/",
    "../SR_MC/ZvvJets_HT800to1200/",
    "../SR_MC/ZvvJets_HT1200to2500/",
    "../SR_MC/ZvvJets_HT2500toInf/",
  },
};

const vector<string> inFileNameSignal = {
  "../Signal/Wino_M_300_cTau_3/",
  "../Signal/Wino_M_300_cTau_10/",
  "../Signal/Wino_M_300_cTau_30/",
  "../Signal/Wino_M_500_cTau_10/",
  "../Signal/Wino_M_500_cTau_20/",
  "../Signal/Wino_M_650_cTau_10/",
  "../Signal/Wino_M_650_cTau_20/",
  "../Signal/Wino_M_800_cTau_10/",
  "../Signal/Wino_M_800_cTau_20/",
  "../Signal/Wino_M_1000_cTau_10/",
  "../Signal/Wino_M_1000_cTau_20/",
};

const vector<vector<string>> inFileNameData = {
  {
    "../SR_DATA/MET_Run2017B_31Mar2018/",
    "../SR_DATA/MET_Run2017C_31Mar2018/",
    "../SR_DATA/MET_Run2017D_31Mar2018/"
  }
};

enum ESignal{
  kWino_M_300_cTau_3,
  kWino_M_300_cTau_10,
  kWino_M_300_cTau_30,
  kWino_M_500_cTau_10,
  kWino_M_500_cTau_20,
  kWino_M_650_cTau_10,
  kWino_M_650_cTau_20,
  kWino_M_800_cTau_10,
  kWino_M_800_cTau_20,
  kWino_M_1000_cTau_10,
  kWino_M_1000_cTau_20,
  kNsignals
};

const double signalCrossSectionTwoTracks[kNsignals] = { // (fb)
  190,  // wino m=300 cτ=3
  190,  // wino m=300 cτ=10
  190,  // wino m=300 cτ=30
  22,   // wino m=500 cτ=10
  22,   // wino m=500 cτ=20
  6.4,  // wino m=650 cτ=10
  6.4,  // wino m=650 cτ=20
  2.2,  // wino m=800 cτ=10
  2.2,  // wino m=800 cτ=20
  0.62, // wino m=1000 cτ=10
  0.62, // wino m=1000 cτ=20
};

const double signalCrossSectionOneTrack[kNsignals] = { // (fb)
  380,  // wino m=300 cτ=3
  380,  // wino m=300 cτ=10
  380,  // wino m=300 cτ=30
  45,   // wino m=500 cτ=10
  45,   // wino m=500 cτ=20
  13,  // wino m=650 cτ=10
  13,  // wino m=650 cτ=20
  4.6,  // wino m=800 cτ=10
  4.6,  // wino m=800 cτ=20
  1.3, // wino m=1000 cτ=10
  1.3, // wino m=1000 cτ=20
};

enum EBackground{
  kQCD,
  kZmumuJets,
  kTT,
  kVV,
  kWmunuJets,
  kZnunuJets,
  kNbackgrounds
};

enum EData{
  kElectron_Run2017B,
  kNdata
};

// Constants for tracker layers
const int nLayers = 14;
const double layerR[nLayers] = { 29, 68, 109, 160, 250, 340, 430, 520, 610, 696, 782, 868, 965, 1080 };

enum EVar{
  kCustom,
  
  // per event variables
  kNvertices,
  kNisoTracks,
  kNjets,
  kNjets30,
  kNjets30a,
  kMetSumEt,
  kMetPt,
  kMetMass,
  kMetEta,
  kMetPhi,
  kMetJetDphi,
  
  // per track variables
  kTrackNclusters,
  kTrackTotalDedx,
  kTrackDedxPerCluster,
  kTrackPt,
  kTrackEta,
  kTrackPhi,
  kTrackCaloEm,
  kTrackCaloHad,
  kTrackDxy,
  kTrackDz,
  kTrackCharge,
  kTrackMass,
  kTrackPid,
  kTrackMissingOuterTrackerHits,
  kTrackPixelHits,
  kTrackTrackerHits,
  kTrackRelativeIsolation,
  kTrackAbsoluteIsolation,
  kTrackMetDphi,
  kTrackDedxPerHit,
  
  // per jet variables
  kJetPt,
  kJetEta,
  kJetPhi,
  kJetTrackDr,
  kJetCHF,
  kJetNHF,
  
  // per track per layer variables
  kDedx,  ///< dE/dx per layer
  kSizeX, ///< X cluster size in each layer
  kSizeY  ///< Y cluster size in each layer
};

inline int BackColor(EBackground bck){
  return TColor::GetColor(backColors[bck][0],backColors[bck][1],backColors[bck][2]);
}

inline int SignalColor(ESignal sig){
  return TColor::GetColor(signalColors[sig][0],signalColors[sig][1],signalColors[sig][2]);
}

inline int DataColor(EData data){
  return TColor::GetColor(dataColors[data][0],dataColors[data][1],dataColors[data][2]);
}

//,,,{255, 215, 180},{0, 0, 128},{128, 128, 128},{255, 255, 255},{0, 0, 0}

//{2,63,165},{125,135,185},{190,193,212},{214,188,192},{187,119,132},{142,6,59},{74,111,227},{133,149,225},{181,187,227},{230,175,185},{224,123,145},{211,63,106},{17,198,56},{141,213,147},{198,222,199},{234,211,198},{240,185,141},{239,151,8},{15,207,192},{156,222,214},{213,234,231},{243,225,235},{246,196,225},{247,156,212}

//,{157,204,0},{194,0,136},{0,51,128},{255,164,5},{255,168,187},{66,102,0},{255,0,16},{94,241,242},{0,153,143},{224,255,102},{116,10,255},{153,0,0},{255,255,128},{255,255,0},{255,80,5}, {240,163,255}, {0,117,220}, {153,63,0}, {76,0,92}, {25,25,25},{0,92,49}, {43,206,72}, {255,204,153}, {128,128,128}, {148,255,181}, {143,124,0}


const map<EVar, tuple<string, int, double, double>> settings =
{
  {kNvertices , {"N good vertices",25,0,100}},
  {kNisoTracks , {"N iso tracks",20,0,10}},
  {kNjets , {"N jets",15,0,15}},
  {kNjets30 , {"N jets with pt > 30, |eta|<2.4",15,0.0,15}},
  {kNjets30a , {"N jets with pt > 30, |eta|<4.7",15,0.0,15}},
  {kMetSumEt , {"MET sum Et",100,-20,5000}},
  {kMetPt , {"MET pT",25,200,1000}},
  {kMetMass , {"MET mass",100,-10e-6,10e6}},
  {kMetEta , {"MET eta",100,-3.5,3.5}},
  {kMetPhi , {"MET phi",100,-3.5,3.5}},
  {kMetJetDphi , {"#Delta #phi (p_{T}^{jet},p_{T}^{MET})",25,-3.5,3.5}},
  
  {kTrackNclusters , {"N detIDs per track",20,0,22}},
  {kTrackTotalDedx , {"total dedx per track",50,0,140}},
  {kTrackDedxPerCluster , {"total dedx per track / n clusters",50,0,14}},
  {kTrackPt , {"Track p_{T} (GeV)",25,0,1000}},
  {kTrackEta , {"Track #eta",50,-3.0,3.0}},
  {kTrackPhi , {"Track #phi",50,-3.5,3.5}},
  {kTrackCaloEm , {"EM calo energy",20,0,10}},
  {kTrackCaloHad , {"Hadron calo energy",20,0,10}},
  {kTrackDxy , {"Displacement in XY",100,-0.02,0.02}},
  {kTrackDz , {"Displacement in Z",100,-0.02,0.02}},
  {kTrackCharge , {"Charge dist",100,-10,10}},
  {kTrackMass , {"Mass dist",500,0.0,0.25}},
  {kTrackPid , {"PDG PID",441,-220,220}},
  {kTrackMissingOuterTrackerHits , {"Missing outer tracker hits",20,0,20}},
  {kTrackPixelHits , {"N pixel hits",10,0,10}},
  {kTrackTrackerHits , {"N tracker hits",40,0,40}},
  {kTrackRelativeIsolation , {"Relative isolation in dR=0.3",20,0,1}},
  {kTrackAbsoluteIsolation , {"Absolute isolation in dR=0.3",20,0,1}},
  {kTrackMetDphi , {"#Delta #phi (p_{T}^{track}},p_{T}^{MET})",25,-3.5,3.5}},
  {kTrackDedxPerHit , {"dE/dx per hit",50,0,10}},
  
  {kJetPt , {"Jet pt",50,0.0,1000}},
  {kJetEta , {"Jet eta",50,-3.0,3.0}},
  {kJetPhi , { "Jet phi",50,-3.5,3.5}},
  {kJetTrackDr , {"#Delta R(jet, track)",100,0,10}},
  {kJetCHF ,{"Jet f_{CH}",100,0,1.0}},
  {kJetNHF , {"Jet f_{NH}",100,0,1.0}},
  {kDedx , {"dedx",50,0,13}},
  {kSizeX , {"sizeX",10,0,13}},
  {kSizeY , {"sizeY",10,0,13}},
};

inline bool IsPerEventVariable(EVar var)
{
  if(var == kNvertices || var == kNisoTracks || var == kNjets || var == kNjets30 ||
     var == kNjets30a || var == kMetSumEt || var == kMetPt || var == kMetMass || var == kMetEta || var == kMetPhi)
    return true;
  return false;
}

inline bool IsPerTrackVariable(EVar var)
{
  if(var == kTrackNclusters || var == kTrackTotalDedx || var == kTrackDedxPerCluster || var == kTrackPt ||
     var == kTrackEta || var == kTrackPhi || var == kTrackCaloEm || var == kTrackCaloHad ||
     var == kTrackDxy || var == kTrackDz || var == kTrackCharge || var == kTrackMass || var == kTrackPid ||
     var == kTrackMissingOuterTrackerHits || var == kTrackPixelHits || var == kTrackTrackerHits ||
     var == kTrackRelativeIsolation || var == kTrackAbsoluteIsolation || var == kTrackMetDphi ||
     var == kTrackDedxPerHit)
    return true;
  return false;
}

inline bool IsPerJetVariable(EVar var)
{
  if(var == kJetPt || var == kJetEta || var ==  kJetPhi || var ==  kJetTrackDr || var ==  kJetCHF ||
     var ==  kJetNHF  || var == kMetJetDphi)
    return true;
  return false;
}

inline bool IsPerLayerVariable(EVar var)
{
  if(var == kDedx || var == kSizeX || var == kSizeY)
    return true;
  return false;
};

template<typename ContainerType>
inline void EraseFast(ContainerType &container, size_t index)
{
  // ensure that we're not attempting to access out of the bounds of the container.
  assert(index < container.size());
  
  //Swap the element with the back element, except in the case when we're the last element.
  if (index + 1 != container.size())
    std::swap(container[index], container.back());
  
  //Pop the back of the container, deleting our old element.
  container.pop_back();
}

template <class T>
class range
{
public:
  range(T _min=-inf, T _max=inf) : min(_min), max(_max){
    if(min > max){
      throw logic_error("You try to set min grater than max!");
    }
  }
  
  inline T GetMin(){return min;}
  inline T GetMax(){return max;}
  
  inline bool IsInside(T val){
    if(val >= min && val <= max) return true;
    else return false;
  }
  
  inline bool IsOutside(T val){
    if(val < min || val > max) return true;
    else return false;
  }
    
private:
  T min;
  T max;
};

struct Point
{
  Point(double _x, double _y, double _z, double _val=0) : x(_x), y(_y), z(_z), val(_val) {}
  double x,y,z, val;
  void Print(){cout<<"("<<x<<","<<y<<","<<z<<")"<<endl;}
  double distance(Point p){return sqrt(pow(x-p.x,2)+pow(y-p.y,2)+pow(z-p.z,2));}
  bool isPionHit = false;
  
  void PerpendicularShift(double R,double c, double tShift, int charge=1){
    int xSign=1, ySign=1;
    if(x> 0 && y> 0){xSign= 1; ySign=-1;}
    if(x<=0 && y> 0){xSign= 1; ySign= 1;}
    if(x<=0 && y<=0){xSign=-1; ySign= 1;}
    if(x> 0 && y<=0){xSign=-1; ySign=-1;}
    double dx = x, dy=y;
    
    // the charge may be inverted here... to be checked later
    x +=  charge * xSign * R/sqrt(pow(dx/dy,2)+1);
    y +=  charge * ySign * R/sqrt(pow(dy/dx,2)+1);
    z += tShift*c;
  }
};

struct Helix
{
  Helix(double _R, double _c, double _x0, double _y0, double _z0, double _thickness)
  : R(_R), c(_c), x0(_x0), y0(_y0), z0(_z0), thickness(_thickness) {}
  
  Helix(double _R, double _c, Point p, double _thickness)
  : R(_R), c(_c), x0(p.x), y0(p.y), z0(p.z), thickness(_thickness) {}
  
  Helix(){}
  
  void Print(){cout<<"R:"<<R<<"\tc:"<<c<<"\toffset:("<<x0<<","<<y0<<","<<z0<<")\tnPoints:"<<nPoints<<"\tnPionPoints:"<<nPionPoints<<endl;}
  
  /// Returns vector of points along helix trajectory that hit the tracker
  /// \param tMin t parameter minimum
  /// \param tMax t parameter maximum
  /// \param tStep t parameter step
  /// \param threshold // how close to the tracker layer hits must be
  vector<Point> GetPointsHittingSilicon(double tMin=0, double tMax=5*2*TMath::Pi(), double tStep=0.01, double threshold = 1.0)
  {
    vector<Point> points;
    
    double x,y,z;
    for(double t=tMin;t<tMax;t+=tStep){
      x = R*cos(t) + x0;
      y = R*sin(t) + y0;
      z = c*t      + z0;
      
      for(int iLayer=0;iLayer<5/*nLayers*/;iLayer++){
        if(fabs(sqrt(x*x+y*y)-layerR[iLayer]) < threshold){
          points.push_back(Point(x,y,z));
        }
      }
    }
    return points;
  }
  
  Point GetClosestPoint(Point p)
  {
    double t = atan2(p.y-y0, p.x-x0);
    
    double x = R*cos(t) + x0;
    double y = R*sin(t) + y0;
    double z = c*t      + z0 + 2*c*TMath::Pi();
    
    double absC = fabs(c);
    
    while(fabs(z-p.z) >= absC){
      if(z < p.z) z += absC;
      else        z -= absC;
    }
    
    return Point(x,y,z);
  }
  
  vector<Point>* GetMatchingPoints(vector<Point> &points)
  {
    vector<Point> *result = new vector<Point>();
    
    for(Point p : points){
      Point q = GetClosestPoint(p);
      double d = p.distance(q);
      if(d < thickness) result->push_back(p);
    }
    return result;
  }
  
  void CountMatchingPoints(vector<Point> &points)
  {
    nPoints = 0;
    nPionPoints = 0;
    
    for(Point p : points){
      Point q = GetClosestPoint(p);
      double d = p.distance(q);
      
      if(d < thickness){
        nPoints++;
        if(p.isPionHit) nPionPoints++;
      }
    }
  }
  
  double thickness;
  double R,c,x0,y0,z0;
  int nPoints = 0;
  int nPionPoints = 0;
};

struct Circle
{
  Circle(double _x, double _y, double _R) : x(_x), y(_y), R(_R) {}
  Circle(){}
  
  void Print(){cout<<"Circle x:"<<x<<"\ty:"<<y<<"\tR:"<<R<<endl;}
  
  int GetNbinsOverlappingWithHist(TH2D *hist)
  {
    int nPoints = 0;
    TH2D *circle = new TH2D(*hist);
    
    for(int binX=0;binX<circle->GetNbinsX();binX++){
      for(int binY=0;binY<circle->GetNbinsY();binY++){
        circle->SetBinContent(binX,binY, 0);
      }
    }
    for(double t=0;t<2*TMath::Pi();t+=0.01){
      circle->Fill(x + R*cos(t),y + R*sin(t));
    }
    
    for(int binX=0;binX<circle->GetNbinsX();binX++){
      for(int binY=0;binY<circle->GetNbinsY();binY++){
        if(circle->GetBinContent(binX,binY) > 0 &&
           hist->GetBinContent(binX,binY) > 0){
          nPoints++;
        }
      }
    }
    return nPoints;
  }
  
  TArc* GetArc(){return new TArc(x,y,R);}
  
  double x,y,R;
};

#endif /* Helpers_h */
