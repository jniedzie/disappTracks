//
//  Point.hpp
//
//  Created by Jeremi Niedziela on 18/12/2018.
//

#ifndef Point_hpp
#define Point_hpp

#include "Helpers.hpp"

class Point
{
public:
  Point(double _x, double _y, double _z, double _val=0);
  
  /// Prints basic info about the point
  void Print();
  
  /// Returns distance between this and another point
  double distance(Point p);
  double distanceXY(Point p);
  
  bool isPionHit = false;
  double x,y,z, val;
  
  inline double GetVectorSlopeC(){
    return tan(TMath::Pi()/2.-acos(z/sqrt(x*x+y*y+z*z)));
  }
  
private:
  double PerpendicularShift(double R,double c, int charge=1);
  
};

#endif /* Point_hpp */
