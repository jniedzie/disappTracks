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
#include "Circle.hpp"
#include "ConfigManager.hpp"

class Helix
{
public:
  /// Constructor taking as an input origin and momentum vector (will be automatically shifted to begin in the origin point)
  /// \param _origin Helix origin point (e.g. chargino decay point)
  /// \param _momentum  Momentum of the particle that creates a helix
  /// \param _charge Charge of the particle (determines helix direction)
  /// \param _config Fitter config
  Helix(const unique_ptr<Point> &_origin,
        const unique_ptr<Point> &_momentum,
        int _charge,
        shared_ptr<ConfigManager> _config);
  
  /// Prints basic information about the helix
  void Print();
  
  /// It will pick only points that are on the helix (within its thickness)
  /// and count how many of them are pion points
  void SetPoints(const shared_ptr<vector<Point>> _points);
  
  // Getters
  shared_ptr<vector<Point>>  GetPoints(){return points;}
  
  inline unique_ptr<Point>  GetOrigin(){return make_unique<Point>(*origin);}
  inline unique_ptr<Point>  GetMomentum(){return make_unique<Point>(*momentum);}
  inline double   GetRadius(){return radius;}
  inline double   GetSlope(){return slope;}
  inline int      GetCharge(){return charge;}
  
  inline double   GetTmin(){return tShift;}
  inline double   GetTmax(){return tMax;}
  inline double   GetTstep(){return tStep;}
  inline int      GetNpoints(){return (int)points->size();}
  inline int      GetNpionPoints(){return nPionPoints;}
  inline int      GetNregularPoints(){return nRegularPoints;}
  
  inline double   GetNcycles(){return sgn(momentum->GetZ())*((sgn(momentum->GetZ())*pixelBarrelZsize) - origin->GetZ())/(fabs(slope)*2*TMath::Pi());}
  
  /// Calculates average of the squared distances between points (hits) and the helix
  double GetChi2();
  
private:
  shared_ptr<vector<Point>> points;   ///< Vector of points laying on the helix (withing thickness)
  double tShift;          ///< Angle by which beginning of the helix is shifted due to the shift of its origin
  double tMax;            ///< Max angle (taking into account number of cycles
  double tStep;           ///< Step determining drawing precision
  
  int nRegularPoints = 0; ///< Number of points that are distributed regularly along Z axis
  int nPionPoints = 0;    ///< Number of points along the helix that are true pion hits
  
  unique_ptr<Point> origin;     ///< Center of the helix
  unique_ptr<Point> momentum;   ///< Pion's momentum vector
  double radius;                ///< Radius of the helix
  double slope;                 ///< Slope of the helix in Z direction
  double slopeAbs;              ///< Absolute value of the slope (to speed up the calculation)
  int    charge;                ///< Charge of the particle (determines helix direction)
  
  Point GetClosestPoint(Point p);
  shared_ptr<ConfigManager> config;
  
  unique_ptr<PointsProcessor> pointsProcessor;
  
  friend class HelixProcessor;
};

#endif /* Helix_hpp */
