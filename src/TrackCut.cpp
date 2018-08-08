//
//  TrackCut.cpp
//
//  Created by Jeremi Niedziela on 16/07/2018.
//

#include "TrackCut.hpp"

#include <iostream>

using namespace std;

TrackCut::TrackCut(ECut cutType) :
minDedxClusters(0),
maxDedxClusters(999999),
minDedxPerCluster(0.0),
minTotalDeDx(0.0),
maxTotalDeDx(999999),
minPt(0.0),
maxEmCalo(999999),
maxHadCalo(999999),
maxEta(999999)
{
  if(cutType&kEmpty) return;
  if(cutType&kPt50GeV)  minPt = 50.0;
  if(cutType&kPt200GeV) minPt = 200.0;
  if(cutType&kLowCalo){
      maxEmCalo = 0.5;
      maxHadCalo = 0.5;
  }
  if(cutType&kLowDEdx) maxTotalDeDx = 38.0;
  if(cutType&kShort){
      minDedxClusters = 3;
      maxDedxClusters = 4;
  }
  if(cutType&kMedium){
      minDedxClusters = 3;
      maxDedxClusters = 8;
  }
  if(cutType&kEta2p4) maxEta = 2.4;
}

TrackCut::~TrackCut()
{
  
}
