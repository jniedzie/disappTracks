//
//  Circle.cpp
//
//  Created by Jeremi Niedziela on 18/12/2018.
//

#include "Circle.hpp"

Circle::Circle(double _x, double _y, double _R) : x(_x), y(_y), R(_R)
{
//  tShift = acos(1/sqrt(pow(x/y,2)+1));
}

Circle::Circle()
{
  
}

Circle& Circle::operator=(const Circle &c)
{
  x = c.x;
  y = c.y;
  R = c.R;
  points = c.points;
  chi2 = c.chi2;
  tShift = c.tShift;
  return *this;
}

void Circle::Print()
{
  cout<<"Circle x:"<<x<<"\ty:"<<y<<"\tR:"<<R<<endl;
  
}

void Circle::Shift(int charge)
{
  int xSign=1, ySign=1;
  if(x> 0 && y> 0){xSign= 1; ySign=-1;}
  if(x<=0 && y> 0){xSign= 1; ySign= 1;}
  if(x<=0 && y<=0){xSign=-1; ySign= 1;}
  if(x> 0 && y<=0){xSign=-1; ySign=-1;}
  
  xInit=x;
  yInit=y;
  
  // the charge may be inverted here... to be checked later
  x += charge * xSign * R/sqrt(pow(xInit/yInit,2)+1);
  y += charge * ySign * R/sqrt(pow(yInit/xInit,2)+1);
  
//  tShift = atan2((yInit+y)/(yInit-y)-x,y);
}

void Circle::ShiftByVector(Point v, int charge)
{
  shiftVector = v;
  v = Point(charge * -v.y,charge * v.x, v.z); // take a vector perpendicular to the pion's momentum vector
  
  double vTransverseLength = sqrt(v.x*v.x+v.y*v.y);
  tShift = acos(-v.x/vTransverseLength);
  
  double scale = R/vTransverseLength;
  
  v.x *= scale;
  v.y *= scale;
  
  x += v.x;
  y += v.y;
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
    circle->Fill(x + R*cos(t),y + R*sin(t));
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

Point Circle::GetClosestPoint(Point p)
{
  double t = atan2(p.y-y, p.x-x);
  double xx = R*cos(t) + x;
  double yy = R*sin(t) + y;
  return Point(xx,yy,0);
}

TArc* Circle::GetArc()
{
  return new TArc(x,y,R);
}

