//
//  Helix.hpp
//
//  Created by Jeremi Niedziela on 18/12/2018.
//

#ifndef Helix_hpp
#define Helix_hpp

#include "Helpers.hpp"
#include "Point.hpp"
#include "Circle.hpp"

class Helix
{
public:
  Helix();
  
  Helix(double _R, double _c, double _x0, double _y0, double _z0, int _nCycles, double _thickness);

  Helix(double _R, double _c, Point p, int _nCycles, double _thickness);
  
  Helix(double _c, Circle circle, int _nCycles, double _thickness);
  
  void Print();
  
  void Shift(int charge=1);
  
  void ShiftByVector(Point v, int charge=1);
  
  /// Returns vector of points along helix trajectory that hit the tracker
  vector<Point> GetPointsHittingSilicon();
  
  Point GetClosestPoint(Point p);
  
  vector<Point>* GetMatchingPoints(vector<Point> &points);
  
  void CountMatchingPoints(const vector<Point> &points);
  
  double thickness;
  double R,c,x0,y0,z0;
  int nPoints = 0;
  int nPionPoints = 0;
  vector<Point> points;
  double tShift;
  double tMin, tMax, tStep;
  double chi2 = inf;
  int nCycles;
  int xSign, ySign;
};

#endif /* Helix_hpp */
