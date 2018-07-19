//
//  HistSet.hpp
//  disappTracksTarget
//
//  Created by Jeremi Niedziela on 16/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#ifndef HistSet_hpp
#define HistSet_hpp

#include "Event.hpp"

#include <TH1D.h>
#include <TCanvas.h>

#include <vector>

class HistSet {
public:
  enum EVar{
    kCustom,
    
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
    
    // per jet variables
    kJetPt,
    kJetEta,
    kJetPhi,
    
    // per track per layer variables
    kDedx,  ///< dE/dx per layer
    kSizeX, ///< X cluster size in each layer
    kSizeY  ///< Y cluster size in each layer
  };
  
  HistSet(const char* title, int nBins, double min, double max);
  HistSet(EVar _var);
  ~HistSet();
  
  void FillSignal(double value){signal->Fill(value);}
  void FillBackground(double value){background->Fill(value);}
  void FillData(double value){data->Fill(value);}
  
  void FillFromEvents(Events *signalEvents, Events *backgroundEvents, Events *dataEvents);
  
  void Draw(TCanvas *c1, int pad);
  void DrawPerLayer();
  
private:
  TH1D *signal;
  TH1D *background;
  TH1D *data;
  
  EVar var;
  const char* customTitle;
  
  std::vector<TH1D*> signalPerLayer;
  std::vector<TH1D*> backgroundPerLayer;
  std::vector<TH1D*> dataPerLayer;
  
  void FillFromEventsPerLayer(Events *signalEvents, Events *backgroundEvents, Events *dataEvents);
  void FillFromEventsGlobal(Events *signalEvents, Events *backgroundEvents, Events *dataEvents);
  
  TLegend* GetLegend(double legendW = 0.15, double legendH = 0.5, double legendX = 0.75, double legendY = 0.25,const char* header="");
  
  const char* GetTitle();
  int GetNbins();
  double GetMin();
  double GetMax();
  bool ShouldNormalize();
  bool DoSumw2();
  
  void Fill(TH1D* hist, Events *events, int iLayer=-1);
};

#endif /* HistSet_hpp */
