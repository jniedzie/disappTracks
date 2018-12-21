//
//  Circle.cpp
//
//  Created by Jeremi Niedziela on 18/12/2018.
//

#include "Circle.hpp"

Circle::Circle(const unique_ptr<Point> &_decayPoint,
               const unique_ptr<Point> &_momentum,
               int _charge, double _thickness) :
decayPoint(make_unique<Point>(*_decayPoint)),
momentum(make_unique<Point>(*_momentum)),
charge(_charge),
thickness(_thickness)
{
  radius = GetRadiusInMagField(momentum->GetX(), momentum->GetY(), solenoidField);
  
  // take a vector perpendicular to the pion's momentum vector
  Point v = Point(charge * -momentum->GetY(),charge * momentum->GetX(), momentum->GetZ());
  
  double vTransverseLength = sqrt(pow(v.GetX(),2)+pow(v.GetY(),2));
  tShift = acos(-v.GetX()/vTransverseLength);
  double scale = radius/vTransverseLength;
  
  center = make_unique<Point>(*decayPoint);
  center->SetX(center->GetX() + scale*v.GetX());
  center->SetY(center->GetY() + scale*v.GetY());
}

void Circle::Print()
{
  cout<<"Circle center ("<<center->GetX()<<","<<center->GetY()<<")\tR:"<<radius<<endl;
}

int Circle::GetNbinsOverlappingWithHist(TH2D *hist)
{
  int nPoints = 0;
  unique_ptr<TH2D> circle = unique_ptr<TH2D>(new TH2D(*hist));
  
  for(int binX=0;binX<circle->GetNbinsX();binX++){
    for(int binY=0;binY<circle->GetNbinsY();binY++){
      circle->SetBinContent(binX,binY, 0);
    }
  }
  for(double t=0;t<2*TMath::Pi();t+=0.01){
    circle->Fill(center->GetX() + radius*cos(t),center->GetY() + radius*sin(t));
  }
  
  for(int binX=0;binX<circle->GetNbinsX();binX++){
    for(int binY=0;binY<circle->GetNbinsY();binY++){
      if(circle->GetBinContent(binX,binY) > 0 &&
         hist->GetBinContent(binX,binY) > 0){
        nPoints++;
      }
    }
  }
  return nPoints;
}

double Circle::GetDistanceToPoint(Point p)
{
  double t = atan2(p.GetY()-center->GetY(), p.GetX()-center->GetX());
  double x = radius*cos(t) + center->GetX();
  double y = radius*sin(t) + center->GetY();
  return sqrt(pow(p.GetX()-x,2)+pow(p.GetY()-y,2));
}

void Circle::SetPoints(vector<Point> _points)
{
  points.clear();
  for(auto p : _points){
    if(GetDistanceToPoint(p) < thickness) points.push_back(p);
  }
}

TArc* Circle::GetArc()
{
  return new TArc(center->GetX(),center->GetY(),radius);
}

