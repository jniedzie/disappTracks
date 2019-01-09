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
  /// Default constructor taking XYZ coordinates and optionally a value in this point
  Point(double _x, double _y, double _z, double _value=0);
  
  /// Constructs a point that is an average of provided vector of points
  Point(vector<Point> points);
  
  /// Prints basic info about the point
  void Print();
  
  /// Returns distance between this and another point
  double distance(Point p);
  
  /// Returns distance between this and another point in the transverse plane
  double distanceXY(Point p);
  
  /// Returns slope of vector defined by this point calculated from the Y axis
  double GetVectorSlopeC();
  
  /// Tells whether or not this point belongs to a true pion's helix
  inline bool IsPionHit(){return isPionHit;}
  
  /// Separates vector of points into groups with the same XY position (within tolerance)
  static vector<vector<Point>> SplitPointsIntoLines(vector<Point> points, double tolerance);
  
  /// Returns a vector filled with random points in the pixel barrel
  /// \param nPoints Number of points that will be generated
  static vector<Point> GetRandomPoints(int nPoints);
  
  // Trivial getters
  inline double GetX(){return x;}
  inline double GetY(){return y;}
  inline double GetZ(){return z;}
  inline double GetValue(){return value;}
  
  // Trivial setters
  inline void SetX(double val){x = val;}
  inline void SetY(double val){y = val;}
  inline void SetZ(double val){z = val;}
  inline void SetIsPionHit(bool val){isPionHit = val;}
private:
  double x,y,z; ///< XYZ coordinates of the point
  double value; ///< Value at this point
  
  bool isPionHit; ///< Flag saying whether or not this point belongs to a true pion's helix
  
};

#endif /* Point_hpp */
