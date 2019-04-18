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
//  a.circleProcessor = make_unique<CircleProcessor>();
  for(auto &c : a.circles){circles.push_back(make_unique<Circle>(c));}
  for(auto &p : a.points){points.push_back(make_shared<Point>(*p));}
  iCycle = a.iCycle;
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
  cout<<"\tPoints:"<<endl;
  for(auto point : points){
    cout<<"\t";point->Print();cout<<endl;
  }
  cout<<"Radii chi2:"<<GetRadiiSlopeChi2()<<endl;
}

void ArcSet2D::AddCircle(const unique_ptr<Circle> &circle)
{
  // calculate phi range of this segment
  double phiMin=-inf, phiMax=inf;
  
  if(circle->GetNpoints() == 3){
    double phiStart  = circle->GetPointAngle(0);
    double phiMiddle = circle->GetPointAngle(1);
    double phiEnd    = circle->GetPointAngle(2);
    
    if(phiEnd > phiStart &&
       phiEnd > phiMiddle){
      
      iCycle++;
      phiEnd    -= iCycle * 2*TMath::Pi();
      phiMiddle -= iCycle * 2*TMath::Pi();
    }
    else if(phiEnd > phiMiddle &&
            phiEnd < phiStart){
      iCycle++;
      phiEnd -= iCycle * 2*TMath::Pi();
    }
    else if(phiEnd    > phiStart &&
            phiMiddle > phiStart){
      
      phiMiddle -= iCycle * 2*TMath::Pi();
      phiEnd    -= iCycle * 2*TMath::Pi();
    }
    else{
      phiStart  -= iCycle * 2*TMath::Pi();
      phiMiddle -= iCycle * 2*TMath::Pi();
      phiEnd    -= iCycle * 2*TMath::Pi();
    }
    
    phiMin = min(min(phiMiddle, phiEnd), phiStart)/TMath::Pi() * 180;
    phiMax = max(max(phiMiddle, phiEnd), phiStart)/TMath::Pi() * 180;
    
  }
  else if(circle->GetNpoints() == 2){
    double phiStart  = circle->GetPointAngle(0);
    double phiEnd    = circle->GetPointAngle(1);
    
    if(phiEnd > phiStart){
      iCycle++;
      phiEnd    -= iCycle * 2*TMath::Pi();
    }
    else{
      phiStart  -= iCycle * 2*TMath::Pi();
      phiEnd    -= iCycle * 2*TMath::Pi();
    }
    
    phiMin = min(phiEnd, phiStart)/TMath::Pi() * 180;
    phiMax = max(phiEnd, phiStart)/TMath::Pi() * 180;
  }
  auto r = range<double>(phiMin, phiMax);
  
  
  // add circle to the vector of segments
  circles.push_back(circleProcessor.CopyCircleAddingRange(circle, r));

  // add last point of the circle to the collection (first two should already be there)
  points.push_back(circle->GetLastPoint());
}

vector<TArc*> ArcSet2D::GetArcs()
{
  vector<TArc*> arcs;
  
  for(auto &circle : circles){
    auto arc = circle->GetArc();
    arc->SetPhimin(circle->GetRange().GetMin());
    arc->SetPhimax(circle->GetRange().GetMax());
    
    arcs.push_back(arc);
  }
  
  return arcs;
}

double ArcSet2D::GetOriginPhi()
{
  return circles[0]->GetPointAngle(0);
}

double ArcSet2D::GetRadiiSlopeChi2()
{
  auto graph = make_unique<TGraph>();
  int iter=0;
  
  for(auto &circle : circles){
    graph->SetPoint(iter, iter, circle->GetRadius());
  }
  TF1 *fun = new TF1("fun","[0]*x+[1]",0,graph->GetN());
  graph->Fit(fun);
  
  
  return fun->GetChisquare();
}
