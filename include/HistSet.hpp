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
 
  
  HistSet(const char* title, int nBins, double min, double max);
  HistSet(EVar _var);
  ~HistSet();
  
  inline void FillSignal(ESignal sig, double value){signal[sig]->Fill(value);}
  inline void FillBackground(EBackground bck, double value){background[bck]->Fill(value);}
  inline void FillData(double value){data->Fill(value);}
  
  void FillFromEvents(Events *signalEvents[kNsignals],
                      Events *backgroundEvents[kNbackgrounds],
                      Events *dataEvents);
  
  void Draw(TCanvas *c1, int pad);
  void DrawPerLayer();
  
  inline void SetShowNonZerBinPosX(){showNonZeroBinPosX = true;}
private:
  TH1D *signal[kNsignals];
  TH1D *data;
  
  TH1D *background[kNbackgrounds];
  
  EVar var;
  const char* customTitle;
  bool showNonZeroBinPosX;
  
  std::vector<TH1D*> signalPerLayer[kNsignals];
  std::vector<TH1D*> backgroundPerLayer[kNbackgrounds];
  std::vector<TH1D*> dataPerLayer;
  
  void FillFromEventsPerLayer(Events *signalEvents[kNsignals],
                              Events *backgroundEvents[kNbackgrounds],
                              Events *dataEvents);
  
  TLegend* GetLegend();
  
  const char* GetTitle();
  int GetNbins();
  double GetMin();
  double GetMax();
  bool ShouldNormalize();
  bool DoSumw2();
  
  void Fill(TH1D* hist, Events *events, int iLayer=-1);
  double GetNonZeroBinPosX(TH1D *hist);
};

#endif /* HistSet_hpp */
