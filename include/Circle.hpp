//
//  Circle.hpp
//
//  Created by Jeremi Niedziela on 18/12/2018.
//

#ifndef Circle_hpp
#define Circle_hpp

#include "Helpers.hpp"
#include "Point.hpp"
#include "ConfigManager.hpp"

class Circle
{
public:
  /// Constructor of the circle that will shift it according to the momentum vector and pion's charge
  /// \param _decayPoint Point to which the circle must be tangent
  /// \param _momentum Pion's momentum determining shift and rotation of the circle
  Circle(const unique_ptr<Point> &_decayPoint,
         const unique_ptr<Point> &_momentum,
         const range<double> &_phiRange = range<double>(0, 2*TMath::Pi()));
  
  /// Constructor for the circle that takes its center, point where it begins and a radius
  Circle(const unique_ptr<Point> &_decayPoint,
         const unique_ptr<Point> &_center,
         double _radius,
         const range<double> &_phiRange = range<double>(0, 2*TMath::Pi()));
  
  Circle(const unique_ptr<Circle> &c);
  
  /// Prints basic information about the circle
  void Print();
  
  /// Returns closest distance between point p and this circle
  double GetDistanceToPoint(Point p);
  
  /// Sets points that belong to this circle (will automatically filter out those that are further than thickness)
  void SetPoints(const vector<shared_ptr<Point>> &_points);
  
  inline void SetPz(double pz){momentum->SetZ(pz);}
  
  /// Returns number of points laying on the circle
  int GetNpoints(){return (int)points.size();}
  
  /// Returns vector of points laying on the circle
  vector<shared_ptr<Point>> GetPoints(){return points;}
  
  /// Returns point at index i
  shared_ptr<Point> GetPoint(unsigned long i){return points[i];}
  
  /// Returns ROOT object representing a circle, that can be plotted in a canvas
  TArc* GetArc();
  
  /// Returns chargino's decay point
  unique_ptr<Point> GetDecayPoint(){return make_unique<Point>(*decayPoint);}
  
  /// Returns shifted center of the circle
  unique_ptr<Point> GetCenter(){return make_unique<Point>(*center);}
  
  /// Returns pion's momentum
  unique_ptr<Point> GetMomentum(){return make_unique<Point>(*momentum);}
  
  /// Returns circle's radius calculated from the pion's momentum
  inline double GetRadius() const {return radius;}
  
  /// Returns an angle by which circle is rotated due to the shift of its origin
  inline double GetToffset(){return tShift;}
  
  /// Returns phi angle of circle's point i
  /// \param i Index of the circle's point
  double GetPointAngle(uint i);
  
  inline range<double> GetRange(){return phiRange;}
  
private:
  unique_ptr<Point> decayPoint;   ///< Decay point of the chargino
  unique_ptr<Point> center;       ///< Center of the circle (will be automatically shifted in the constructor)
  unique_ptr<Point> momentum;     ///< Pion's momentum vector
  
  vector<shared_ptr<Point>> points; ///< Points belonging to this circle
  double radius;                    ///< Radius of the circle (calculated from the momentum)
  double tShift;                    ///< Angle by which circle is rotated due to the shift of its origin
  
  range<double> phiRange;  ///< Beginning and end of the circle in phi
  
  friend class CircleProcessor;
  friend class ArcSetProcessor;
};

#endif /* Circle_hpp */
