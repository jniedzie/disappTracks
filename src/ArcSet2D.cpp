//
//  ArcSet2D.cpp
//  xDisappTracks
//
//  Created by Jeremi Niedziela on 18/03/2019.
//  Copyright © 2019 Jeremi Niedziela. All rights reserved.
//

#include "ArcSet2D.hpp"

ArcSet2D::ArcSet2D() :
iCycle(0)
{
  
}

ArcSet2D::ArcSet2D(const ArcSet2D &a)
{
  for(auto &c : a.circles){circles.push_back(make_unique<Circle>(c));}
  for(auto r : a.circlesRanges){circlesRanges.push_back(r);}
  for(auto &p : a.points){points.push_back(make_shared<Point>(p));}
  iCycle = a.iCycle;
}

ArcSet2D::ArcSet2D(const unique_ptr<ArcSet2D> &a)
{
  for(auto &c : a->circles){circles.push_back(make_unique<Circle>(c));}
  for(auto r : a->circlesRanges){circlesRanges.push_back(r);}
  for(auto &p : a->points){points.push_back(make_shared<Point>(p));}
  iCycle = a->iCycle;
}

ArcSet2D::~ArcSet2D()
{
  
}

void ArcSet2D::Print()
{
  cout<<"------------------------"<<endl;
  cout<<"Arc set:"<<endl;
  cout<<"\tCircles:"<<endl;
  for(auto &circle : circles){
    cout<<"\t";circle->Print();cout<<endl;
  }
  cout<<"\tCircle ranges:"<<endl;
  for(auto r : circlesRanges){
    cout<<"\t";r.Print();cout<<endl;
  }
  cout<<"\tPoints:"<<endl;
  for(auto point : points){
    cout<<"\t";point->Print();cout<<endl;
  }
}

void ArcSet2D::AddCircle(const unique_ptr<Circle> &circle)
{
  // add circle to the vector of segments
  circles.push_back(make_unique<Circle>(circle));
  
  // calculate phi range of this segment
  double phiVertex = circle->GetPointAngle(0);
  double phi1      = circle->GetPointAngle(1);
  double phi2      = circle->GetPointAngle(2);
  
  if(phi2 > phiVertex &&
     phi2 > phi1){
  
    iCycle++;
    phi2 -= iCycle * 2*TMath::Pi();
  }
  else if(phi2 > phiVertex &&
          phi1 > phiVertex){
    
    phi1 -= iCycle * 2*TMath::Pi();
    phi2 -= iCycle * 2*TMath::Pi();
  }
  else{
    phiVertex -= iCycle * 2*TMath::Pi();
    phi1      -= iCycle * 2*TMath::Pi();
    phi2      -= iCycle * 2*TMath::Pi();
  }
  
  double phiMin = min(min(phi1, phi2), phiVertex)/TMath::Pi() * 180;
  double phiMax = max(max(phi1, phi2), phiVertex)/TMath::Pi() * 180;
  
  auto r = range<double>(phiMin, phiMax);
  
  circlesRanges.push_back(r);
  
  // add last point of the circle to the collection (first two should already be there)
  points.push_back(circle->GetPoint(2));
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
    iter++;
  }
  return arcs;
}

double ArcSet2D::GetOriginPhi()
{
  return circles[0]->GetPointAngle(0);
}
