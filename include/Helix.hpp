//
//  Helix.hpp
//
//  Created by Jeremi Niedziela on 18/12/2018.
//

#ifndef Helix_hpp
#define Helix_hpp

#include "Helpers.hpp"
#include "Point.hpp"
#include "PointsProcessor.hpp"
#include "CircleProcessor.hpp"
#include "ConfigManager.hpp"
#include "Track.hpp"

class Helix
{
public:
  /// Constructor taking as an input origin and momentum vector (will be automatically shifted to begin in the origin point)
  /// \param _origin Helix origin point (e.g. chargino decay point)
  /// \param _momentum  Momentum of the particle that creates a helix
  /// \param _charge Charge of the particle (determines helix direction)
  Helix(const Point &_origin,
        const unique_ptr<Point> &_momentum,
        int _charge);
  
  // Copy constructor
  Helix(const Helix &h);
  
  /// Prints basic information about the helix
  void Print();
  
  /// It will pick only points that are on the helix (within its thickness)
  /// and count how many of them are pion points
  void SetPoints(const vector<shared_ptr<Point>> &_points);
  
  // Getters
  vector<shared_ptr<Point>>  GetPoints(){return points;}
  
  inline const Point&   GetOrigin(){return origin;}
  inline unique_ptr<Point>  GetMomentum(){return make_unique<Point>(*momentum);}
  inline double   GetRadius(){return radius;}
  inline double   GetSlope(){return slope;}
  inline int      GetCharge(){return charge;}
  
  inline double   GetTmin(){return tShift;}
  inline double   GetTmax(){return tMax;}
  inline double   GetTstep(){return tStep;}
  inline uint     GetNpoints(){return (uint)points.size();}
  inline int      GetNpionPoints(){return nPionPoints;}
  inline int      GetNregularPoints(){return nRegularPoints;}
  
  inline double   GetNcycles(){return sgn(momentum->GetZ())*((sgn(momentum->GetZ())*trackerZsize) - origin.GetZ())/(fabs(slope)*2*TMath::Pi());}
  
  /// Calculates average of the squared distances between points (hits) and the helix
  double GetChi2();
  
  // limit params
  // x = L cos(φ) + (R0 - at) cos(t)
  // y = L sin(φ) + (R0 - at) sin(t)
  // z = L ctg(φ) + (s0 - bt) t
  
  Helix(const Track &_track, const Point &p1, const Point &p2, const Point &_eventVertex);
  double GetRadius(double t){
    double R0 = (R0max + R0min)/2.;
    double a = (amin + amax)/2.;
    return (R0 - a*(t));
  }
  
  double GetSlope(double t){
    double s0 = (s0min + s0max)/2.;
    double b  = (bmin + bmax)/2.;
    return (s0 - b*(t));
  }
  
  double Lmin, Lmax;
  double bmin, bmax;
  double s0min, s0max;
  double amin, amax;
  double R0min, R0max;
  int iCycles;
  bool isFinished = false;
  double slope_valmin=inf, slope_valmax=-inf;
  double radius_valmin=inf, radius_valmax=-inf;
  
  unique_ptr<Point> GetVertex(){return make_unique<Point>(*vertex);}
  bool ExtendByPoint(const Point &point);
  void CalcAndUpdateSlopeVars(double z0, double t0, double z1, double t1, double z2, double t2);
  pair<double, double> CalcSlopeVars(double z0, double t0, double z1, double t1, double z2, double t2);
  
  pair<double, double> CalcRadiiVars(double x0, double t0, double x1, double t1, double x2, double t2);
  
private:
  vector<shared_ptr<Point>> points;   ///< Vector of points laying on the helix (withing thickness)
  double tShift;          ///< Angle by which beginning of the helix is shifted due to the shift of its origin
  double tMax;            ///< Max angle (taking into account number of cycles
  double tStep;           ///< Step determining drawing precision
  
  int nRegularPoints = 0; ///< Number of points that are distributed regularly along Z axis
  int nPionPoints = 0;    ///< Number of points along the helix that are true pion hits
  
  unique_ptr<Point> vertex;     ///< Decay point (beginning) of the helix
  Point origin;                 ///< Center of the helix
  unique_ptr<Point> momentum;   ///< Pion's momentum vector
  Track track;
  double radius;                ///< Radius of the helix
  double slope;                 ///< Slope of the helix in Z direction
  double slopeAbs;              ///< Absolute value of the slope (to speed up the calculation)
  int    charge;                ///< Charge of the particle (determines helix direction)
  
  Point GetClosestPoint(Point p);
  
  unique_ptr<PointsProcessor> pointsProcessor;
  unique_ptr<CircleProcessor> circleProcessor;
  
  friend class HelixProcessor;
};

#endif /* Helix_hpp */
