//
//  HistSet.hpp
//
//  Created by Jeremi Niedziela on 16/07/2018.
//

#ifndef HistSet_hpp
#define HistSet_hpp

#include "EventSet.hpp"

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
  inline void FillData(EData iData, double value){data[iData]->Fill(value);}
  
  void FillFromEvents(vector<shared_ptr<EventSet>> signalEvents,
                      vector<shared_ptr<EventSet>> backgroundEvents,
                      vector<shared_ptr<EventSet>> dataEvents);
  
  void Draw(TCanvas *c1, int pad);
  void DrawPerLayer();
  
  static void DrawStandardPlots(vector<shared_ptr<EventSet>> &eventsSignal,
                                vector<shared_ptr<EventSet>> &eventsBackground,
                                vector<shared_ptr<EventSet>> &eventsData,
                                string prefix="");
  
  static void DrawPerLayerPlots(vector<shared_ptr<EventSet>> &eventsSignal,
                                vector<shared_ptr<EventSet>> &eventsBackground,
                                vector<shared_ptr<EventSet>> &eventsData);
  
  inline void SetShowNonZerBinPosX(){showNonZeroBinPosX = true;}
private:
  vector<TH1D*> signal;
  vector<TH1D*> background;
  vector<TH1D*> data;
  
  EVar var;
  const char* customTitle;
  bool showNonZeroBinPosX;
  
  vector<vector<TH1D*>> signalPerLayer;
  vector<vector<TH1D*>> backgroundPerLayer;
  vector<vector<TH1D*>> dataPerLayer;
  
  void FillFromEventsPerLayer(vector<shared_ptr<EventSet>> signalEvents,
                              vector<shared_ptr<EventSet>> backgroundEvents,
                              vector<shared_ptr<EventSet>> dataEvents);
  
  TLegend* GetLegend();
  
  const char* GetTitle();
  int GetNbins();
  double GetMin();
  double GetMax();
  bool ShouldNormalize();
  bool DoSumw2();
  
  void Fill(TH1D* hist, shared_ptr<EventSet> events, int iDetId=-1);
  double GetNonZeroBinPosX(TH1D *hist);
};

#endif /* HistSet_hpp */
