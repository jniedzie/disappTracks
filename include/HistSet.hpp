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
    // per track variables
    
    
    // per jet variables
    kJetPt,
    kJetEta,
    kJetPhi,
    
    // per track per layer variables
    kDedx,  ///< dE/dx per layer
    kSizeX, ///< X cluster size in each layer
    kSizeY  ///< Y cluster size in each layer
  };
  
  HistSet();
  HistSet(const char* title, int nBins, double min, double max);
  ~HistSet();
  
  void FillSignal(double value){signal->Fill(value);}
  void FillBackground(double value){background->Fill(value);}
  void FillData(double value){data->Fill(value);}
  
  void FillFromEvents(Events *signalEvents, Events *backgroundEvents, Events *dataEvents, EVar var);
  
  void Draw(TCanvas *c1, int pad);
  void DrawPerLayer(EVar var);
  
private:
  TH1D *signal;
  TH1D *background;
  TH1D *data;
  
  std::vector<TH1D*> signalPerLayer;
  std::vector<TH1D*> backgroundPerLayer;
  std::vector<TH1D*> dataPerLayer;
  
  void FillFromEventsPerLayer(Events *signalEvents, Events *backgroundEvents, Events *dataEvents, EVar var);
  void FillFromEventsGlobal(Events *signalEvents, Events *backgroundEvents, Events *dataEvents, EVar var);
  
  TLegend* GetLegend(double legendW = 0.15, double legendH = 0.5, double legendX = 0.75, double legendY = 0.25,const char* header="");
  
  const char* GetTitle(EVar var);
};

#endif /* HistSet_hpp */
