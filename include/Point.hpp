//  Point.hpp
//
//  Created by Jeremi Niedziela on 18/12/2018.

#ifndef Point_hpp
#define Point_hpp

#include "Helpers.hpp"

class Point
{
public:
  /// Default constructor taking XYZ coordinates and optionally a value in this point
  Point(double _x, double _y, double _z, double _value=0, string _subDetName="",
        double _errX=0, double _errY=0, double _errZ=0);
  
  /// Copy constructor
  Point(const Point &p);
  
  /// Constructs a point that is an average of provided vector of points (also errors are averaged)
  Point(vector<Point> points);
  
  /// Comparison operator (points are identical if x, y, z, value are identical within 0.000001 and
  /// sub-detector and isPionHit match
  bool operator==(const Point &p) const;
  
  /// Prints basic info about the point
  void Print() const;
  
  /// Returns slope of vector defined by this point calculated from the Y axis
  double GetVectorSlopeC() const;
  
  /// Tells whether or not this point belongs to a true pion's helix
  inline bool IsPionHit() const {return isPionHit;}
  
  // Trivial getters
  inline double GetX() const {return x;}
  inline double GetY() const {return y;}
  inline double GetZ() const {return z;}
  inline double GetXerr() const {return errX;}
  inline double GetYerr() const {return errY;}
  inline double GetZerr() const {return errZ;}
  inline double GetValue() const {return value;}
  inline string GetSubDetName() const {return subDetName;}
  
  // Trivial setters
  inline void SetX(double val){x = val;}
  inline void SetY(double val){y = val;}
  inline void SetZ(double val){z = val;}
  inline void SetIsPionHit(bool val){isPionHit = val;}
private:
  
  double x,y,z;             ///< XYZ coordinates of the point
  double errX, errY, errZ;  ///< coordinates uncertainties
  double value;             ///< Value at this point
  string subDetName;        ///< Can store name of the sub-detector
  
  bool isPionHit; ///< Flag saying whether or not this point belongs to a true pion's helix
  
  friend class PointsProcessor;
};

#endif /* Point_hpp */
