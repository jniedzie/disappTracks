//
//  Helpers.hpp
//
//  Created by Jeremi Niedziela on 13/06/2018.
//

#ifndef Helpers_h
#define Helpers_h

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>

#include <vector>
#include <iostream>
#include <map>
#include <numeric>

using namespace std;

const bool analyzeData = true;

const int nLayers = 14;
const double layerR[nLayers] = { 29, 68, 109, 160, 250, 340, 430, 520, 610, 696, 782, 868, 965, 1080 };

const int shortTrackMaxNclusters = 3;

inline TLegend* GetLegend(double legendW = 0.15, double legendH = 0.5, double legendX = 0.75, double legendY = 0.25,const char* header="")
{
  
  TLegend *leg = new TLegend(legendX,legendY,legendX+legendW,legendY+legendH);
  leg->SetHeader(header);
  return leg;
}

#endif /* Helpers_h */
