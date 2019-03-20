//
//  ArcSet2D.cpp
//  xDisappTracks
//
//  Created by Jeremi Niedziela on 18/03/2019.
//  Copyright © 2019 Jeremi Niedziela. All rights reserved.
//

#include "ArcSet2D.hpp"

ArcSet2D::ArcSet2D()
{
  
}

ArcSet2D::~ArcSet2D()
{
  
}

void ArcSet2D::AddCircle(const unique_ptr<Circle> &circle, range<double> range)
{
  circles.push_back(make_unique<Circle>(circle));
  circlesRanges.push_back(range);
}

vector<TArc*> ArcSet2D::GetArcs()
{
  vector<TArc*> arcs;
  int iter=0;
  for(auto &circle : circles){
    auto arc = circle->GetArc();
    arc->SetPhimin(circlesRanges[iter].GetMin());
    arc->SetPhimax(circlesRanges[iter].GetMax());
    
    arcs.push_back(arc);
  }
  return arcs;
}
