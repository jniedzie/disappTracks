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
  /// Constructor of the circle that will shift it according to the momentum vector and pion's charge
  /// \param _decayPoint Point to which the circle must be tangent
  /// \param _momentum Pion's momentum determining shift and rotation of the circle
  /// \param _charge Pion's charge determining direction of the shift
  Circle(const unique_ptr<Point> &_decayPoint,
         const unique_ptr<Point> &_momentum,
         int _charge,
         double _thickness);
  
  /// Prints basic information about the circle
  void Print();
  
  /// Calculates number of histogram bins that lay on the circle
  int GetNbinsOverlappingWithHist(TH2D *hist);
 
  /// Returns closest distance between point p and this circle
  double GetDistanceToPoint(Point p);
  
  /// Sets points that belong to this circle (will automatically filter out those that are further than thickness)
  void SetPoints(vector<Point> _points);
  
  /// Returns number of points laying on the circle
  int GetNpoints(){return (int)points.size();}
  
  /// Returns vector of points laying on the circle
  vector<Point> GetPoints(){return points;}
  
  /// Returns ROOT object representing a circle, that can be plotted in a canvas
  TArc* GetArc();
  
  /// Returns chargino's decay point
  unique_ptr<Point> GetDecayPoint(){return make_unique<Point>(*decayPoint);}
  
  /// Returns shifted center of the circle
  unique_ptr<Point> GetCenter(){return make_unique<Point>(*center);}
  
  /// Returns pion's momentum
  unique_ptr<Point> GetMomentum(){return make_unique<Point>(*momentum);}
  
  /// Returns circle's radius calculated from the pion's momentum
  inline double GetRadius(){return radius;}
  
  /// Returns an angle by which circle is rotated due to the shift of its origin
  inline double GetToffset(){return tShift;}
  
private:
  unique_ptr<Point> decayPoint;   ///< Decay point of the chargino
  unique_ptr<Point> center;       ///< Center of the circle (will be automatically shifted in the constructor)
  unique_ptr<Point> momentum;     ///< Pion's momentum vector
  int charge;                     ///< Charge of the particle (determines in which direction to shift)
  double thickness;
  
  vector<Point> points;           ///< Points belonging to this circle
  double radius;                  ///< Radius of the circle (calculated from the momentum)
  double tShift;                  ///< Angle by which circle is rotated due to the shift of its origin
};

#endif /* Circle_hpp */
