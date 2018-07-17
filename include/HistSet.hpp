//
//  HistSet.hpp
//  disappTracksTarget
//
//  Created by Jeremi Niedziela on 16/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#ifndef HistSet_hpp
#define HistSet_hpp

#include <TH1D.h>
#include <TCanvas.h>

class HistSet {
public:
  HistSet(const char* title, int nBins, double min, double max);
  ~HistSet();
  
  void FillSignal(double value){signal->Fill(value);}
  void FillBackground(double value){background->Fill(value);}
  void FillData(double value){data->Fill(value);}
  
  void Draw(TCanvas *c1, int pad);
  
private:
  TH1D *signal;
  TH1D *background;
  TH1D *data;
  
  
  TLegend* GetLegend(double legendW = 0.15, double legendH = 0.5, double legendX = 0.75, double legendY = 0.25,const char* header="");
};

#endif /* HistSet_hpp */
