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
  radius = GetRadiusInMagField(momentum->x, momentum->y, solenoidField);
  
  // take a vector perpendicular to the pion's momentum vector
  Point v = Point(charge * -momentum->y,charge * momentum->x, momentum->z);
  
  double vTransverseLength = sqrt(v.x*v.x+v.y*v.y);
  tShift = acos(-v.x/vTransverseLength);
  double scale = radius/vTransverseLength;
  
  center = make_unique<Point>(*decayPoint);
  center->x += scale*v.x;
  center->y += scale*v.y;
}

void Circle::Print()
{
  cout<<"Circle center ("<<center->x<<","<<center->y<<")\tR:"<<radius<<endl;
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
    circle->Fill(center->x + radius*cos(t),center->y + radius*sin(t));
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
  double t = atan2(p.y-center->y, p.x-center->x);
  double x = radius*cos(t) + center->x;
  double y = radius*sin(t) + center->y;
  return sqrt(pow(p.x-x,2)+pow(p.y-y,2));
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
  return new TArc(center->x,center->y,radius);
}

