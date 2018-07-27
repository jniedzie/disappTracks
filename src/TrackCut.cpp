//
//  TrackCut.cpp
//
//  Created by Jeremi Niedziela on 16/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
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
maxHadCalo(999999)
{
  switch (cutType) {
    case kEmpty:
      break;
    case kHighPt:
      minPt = 200.0;
      break;
    case kLowCalo:
      maxEmCalo = 0.5;
      maxHadCalo = 0.5;
      break;
    case kLowTotal:
      maxTotalDeDx = 38.0;
      break;
    case kShort:
      minDedxClusters = 3;
      maxDedxClusters = 4;
      break;
    case kShortLowCalo:
      minDedxClusters = 3;
      maxDedxClusters = 4;
      maxEmCalo = 0.5;
      maxHadCalo = 0.5;
      break;
    case kShortHighPt:
      minDedxClusters = 3;
      maxDedxClusters = 4;
      minPt = 200.0;
      break;
    case kShortLowTotal:
      minDedxClusters = 3;
      maxDedxClusters = 4;
      maxTotalDeDx = 38.0;
      break;
    case kShortLowTotalHighPt:
      minDedxClusters = 3;
      maxDedxClusters = 4;
      maxTotalDeDx = 38.0;
      minPt = 200.0;
      break;
    case kShortAboveThreshold:
      minDedxClusters = 3;
      maxDedxClusters = 4;
      minDedxPerCluster = 2.5;
      break;
    case kMedium:
      minDedxClusters = 3;
      maxDedxClusters = 8;
      break;
    case kMediumLowCalo:
      minDedxClusters = 3;
      maxDedxClusters = 8;
      maxEmCalo = 0.5;
      maxHadCalo = 0.5;
      break;
    case kMediumHighPt:
      minDedxClusters = 3;
      maxDedxClusters = 8;
      minPt = 200.0;
      break;
    case kMediumLowTotal:
      minDedxClusters = 3;
      maxDedxClusters = 8;
      maxTotalDeDx = 38.0;
      break;
    case kMediumLowTotalHighPt:
      minDedxClusters = 3;
      maxDedxClusters = 8;
      maxTotalDeDx = 38.0;
      minPt = 200.0;
      break;
    default:
      cout<<"ERROR -- no track cut specified... in case you want a blank cut to customize, use ECut::kEmpty."<<endl;
      exit(0);
      break;
  }
}


TrackCut::~TrackCut()
{
  
}
