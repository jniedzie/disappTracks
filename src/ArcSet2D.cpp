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

ArcSet2D::ArcSet2D(const ArcSet2D& a)
{
  for(auto &c : a.circles){circles.push_back(make_unique<Circle>(c));}
  for(auto r : a.circlesRanges){circlesRanges.push_back(r);}
  for(auto &p : a.points){points.push_back(make_shared<Point>(p));}
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
    iter++;
  }
  return arcs;
}

void ArcSet2D::AddPoints(vector<shared_ptr<Point>> p)
{
  points.insert(points.end(), p.begin(), p.end());
}
