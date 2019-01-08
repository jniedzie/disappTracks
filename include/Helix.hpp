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
  /// Constructor taking as an input origin and momentum vector (will be automatically shifted to begin in the origin point)
  /// \param _origin Helix origin point (e.g. chargino decay point)
  /// \param _momentum  Momentum of the particle that creates a helix
  /// \param _charge Charge of the particle (determines helix direction)
  /// \param _nCycles Number of cycles the helix should have
  /// \param _thickness Tolerance to determine whether or not a hit belongs to the helix
  /// \param _zRegularityTolerance Tolerance in Z for regularity calculation
  Helix(const unique_ptr<Point> &_origin,
        const unique_ptr<Point> &_momentum,
        int _charge, int _nCycles, double _thickness, double _zRegularityTolerance);
  
  /// Constructor taking as an input a circle and slope
  /// \param _slope Slope of the helix in Z direction
  /// \param _circle Circle that determines helix radius and center (should be already shifted by a pions vector)
  /// \param _nCycles Number of cycles the helix should have
  /// \param _thickness Tolerance to determine whether or not a hit belongs to the helix
  /// \param _zRegularityTolerance Tolerance in Z for regularity calculation
  Helix(double _slope, const unique_ptr<Circle> &_circle,
        int _nCycles, double _thickness, double _zRegularityTolerance);
  
  /// Prints basic information about the helix
  void Print();
  
  /// Returns vector of points along helix trajectory that hit the tracker
  vector<Point> GetPointsHittingSilicon();
  
  
  void SetPz(double val){momentum->SetZ(val);}
  
  /// It will pick only points that are on the helix (within its thickness)
  /// and count how many of them are pion points
  void SetPoints(const vector<Point> &_points);
  
  // Getters
  vector<Point>*  GetPoints(){return &points;}
  
  inline unique_ptr<Point>  GetOrigin(){return make_unique<Point>(*origin);}
  inline unique_ptr<Point>  GetMomentum(){return make_unique<Point>(*momentum);}
  inline double   GetRadius(){return radius;}
  inline double   GetSlope(){return slope;}
  
  inline double   GetToffset(){return tShift;}
  inline int      GetNpoints(){return (int)points.size();}
  inline int      GetNpionPoints(){return nPionPoints;}
  inline int      GetNregularPoints(){return nRegularPoints;}
  
  int nCycles;
  
  /// Calculates average of the squared distances between points (hits) and the helix
  double GetChi2();
  
private:
  vector<Point> points;   ///< Vector of points laying on the helix (withing thickness)
  double tShift;          ///< Angle by which beginning of the helix is shifted due to the shift of its origin
  double tMax;            ///< Max angle (taking into account number of cycles
  double tStep;           ///< Step determining drawing precision
  
  
  int nRegularPoints = 0; ///< Number of points that are distributed regularly along Z axis
  int nPionPoints = 0;    ///< Number of points along the helix that are true pion hits
  
  double zRegularityTolerance;  ///< Tolerance in Z for regularity calculation
  double thickness;             ///< Tolerance to determine whether or not a hit belongs to the helix
  unique_ptr<Point> origin;     ///< Center of the helix
  unique_ptr<Point> momentum;   ///< Pion's momentum vector
  double radius;                ///< Radius of the helix
  double slope;                 ///< Slope of the helix in Z direction
  int    charge;                ///< Charge of the particle (determines helix direction)
  
  Point GetClosestPoint(Point p);
  vector<vector<Point>> SplitPointsIntoLines();
  
  /// Calculates number of regular points.
  /// Splits all points into lines along Z axis. For each line, checks all possible distances between points.
  /// For each possible distance, calculates number of regular points (for all points in the collection,
  /// within zRegularityTolarance) and finds a maximum number of such points.
  void CalculateNregularPoints();
  
};

#endif /* Helix_hpp */
