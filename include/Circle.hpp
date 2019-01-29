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
  /// \param _config FitterConfig object
  Circle(const unique_ptr<Point> &_decayPoint,
         const unique_ptr<Point> &_momentum,
         const shared_ptr<ConfigManager> _config);
  
  /// Prints basic information about the circle
  void Print();
  
  /// Returns closest distance between point p and this circle
  double GetDistanceToPoint(Point p);
  
  /// Sets points that belong to this circle (will automatically filter out those that are further than thickness)
  void SetPoints(const shared_ptr<vector<Point>> _points);
  
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
  inline double GetRadius() const {return radius;}
  
  /// Returns an angle by which circle is rotated due to the shift of its origin
  inline double GetToffset(){return tShift;}
  
  /// Returns FitterConfig object
  inline shared_ptr<ConfigManager> GetConfig(){return config;}
  
  static void RemoveSimilarCircles(vector<unique_ptr<Circle>> &circles);
private:
  unique_ptr<Point> decayPoint;   ///< Decay point of the chargino
  unique_ptr<Point> center;       ///< Center of the circle (will be automatically shifted in the constructor)
  unique_ptr<Point> momentum;     ///< Pion's momentum vector
  
  vector<Point> points;           ///< Points belonging to this circle
  double radius;                  ///< Radius of the circle (calculated from the momentum)
  double tShift;                  ///< Angle by which circle is rotated due to the shift of its origin
  shared_ptr<ConfigManager> config;
};

#endif /* Circle_hpp */
