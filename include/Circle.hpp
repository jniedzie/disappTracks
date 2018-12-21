//
//  Circle.hpp
//
//  Created by Jeremi Niedziela on 18/12/2018.
//

#ifndef Circle_hpp
#define Circle_hpp

#include "Helpers.hpp"
#include "Point.hpp"

class Circle
{
public:
  
  Circle(double _x, double _y, double _R);
  
  Circle& operator=(const Circle &c);
  
  void Print();
  
  void Shift(int charge=1);
  void ShiftByVector(Point v, int charge=1);
  
  int GetNbinsOverlappingWithHist(TH2D *hist);
 
  Point GetClosestPoint(Point p);
  
  TArc* GetArc();
  
  double xInit, yInit;
  double x,y,R;
  double z;
  vector<Point> points;
  double chi2;
  double tShift;
  Point shiftVector = Point(0,0,0);
  double px,py;
  
private:
  Circle();
  
};

#endif /* Circle_hpp */
