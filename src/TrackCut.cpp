//
//  TrackCut.cpp
//  disappTracksTarget
//
//  Created by Jeremi Niedziela on 16/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#include "TrackCut.hpp"

#include <iostream>

using namespace std;

TrackCut::TrackCut(ECut cutType) :
minDedxClusters(0),
maxDedxClusters(9999),
minDedxPerCluster(0.0),
minTotalDeDx(0.0)
{
  switch (cutType) {
    case kEmpty:
      break;
    case kShort:
      SetNdedxClusters(3, 3);
      break;
    case kShortAboveThreshold:
      SetNdedxClusters(3, 3);
      SetMinDedxPerCluster(5.0);
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
