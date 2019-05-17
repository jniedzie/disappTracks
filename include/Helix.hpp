//  Helix.hpp
//
//  Created by Jeremi Niedziela on 18/12/2018.

#ifndef Helix_hpp
#define Helix_hpp

#include "Helpers.hpp"
#include "Point.hpp"
#include "PointsProcessor.hpp"
#include "CircleProcessor.hpp"
#include "ConfigManager.hpp"
#include "Track.hpp"

struct HelixParams
{
  HelixParams(){}
  HelixParams(double _R0, double _a, double _s0, double _b) : R0(_R0), a(_a), s0(_s0), b(_b) {}
  double R0; // initial radius
  double a;  // radius decrease rate
  double s0; // initial Z-slope
  double b;  // Z-slope decrease rate
};

class Helix
{
public:
  /// Default constructor
  Helix();
  
  /// Constructor taking as an input origin and momentum vector (will be automatically shifted to begin in the origin point)
  /// \param _origin Helix origin point (e.g. chargino decay point)
  /// \param _momentum  Momentum of the particle that creates a helix
  /// \param _charge Charge of the particle (determines helix direction)
  Helix(const Point &_origin,
        const unique_ptr<Point> &_momentum,
        int _charge);
  
  /// Copy constructor
  Helix(const Helix &h);
  
  /// Assignent operator
  Helix operator=(const Helix &h);
  
  /// Comparison operator
  bool operator==(const Helix &h);
  
  /// Prints basic information about the helix
  void Print();
  
  /// It will pick only points that are on the helix (within its thickness)
  /// and count how many of them are pion points
  void SetPoints(const vector<shared_ptr<Point>> &_points);
  
  void SetCharge(int val){charge = val;}
  
  // Getters
  vector<shared_ptr<Point>>  GetPoints() const {return points;}
  
  inline Point  GetOrigin() const {return origin;}
  inline unique_ptr<Point>  GetMomentum() const {return make_unique<Point>(*momentum);}
  inline double   GetRadius() const {return radius;}
  inline double   GetSlope() const {return slope;}
  inline int      GetCharge() const {return charge;}
  
  inline double   GetTmin() const {return tShift;}
  inline double   GetTmax() const {return tMax;}
  inline double   GetTstep() const {return tStep;}
  inline uint     GetNpoints() const {return (uint)points.size();}
  inline int      GetNpionPoints() const {return nPionPoints;}
  inline int      GetNregularPoints() const {return nRegularPoints;}
  
  inline double   GetNcycles() const {return sgn(momentum->GetZ())*((sgn(momentum->GetZ())*trackerZsize) - origin.GetZ())/(fabs(slope)*2*TMath::Pi());}
  
  // limit params
  // x = L cos(φ) + (R0 - at) cos(t)
  // y = L sin(φ) + (R0 - at) sin(t)
  // z = L ctg(φ) + (s0 - bt) t
  
  int iCycles;
  bool isFinished = false;
  uint64_t seedID;
  uint64_t uniqueID;
  
  unique_ptr<Point> GetVertex() const {return make_unique<Point>(*vertex);}
  shared_ptr<Point> GetLastPoint() const {return points[points.size()-1];}
  
  //------------------------------------------------------------------------
  // Another approach
  //
  
  Helix(const HelixParams &_params,
        const Point &_decayVertex,
        const Point &_origin,
        const vector<shared_ptr<Point>> &_points,
        const Track &_track);
  
  HelixParams helixParams;
  
  
  double GetRadius(double t) const { return (helixParams.R0 - helixParams.a*t); }
  double GetSlope(double t) const { return (helixParams.s0 - helixParams.b*t); }
  
  void AddPoint(const shared_ptr<Point> &point);
  
  void SetTmin(double _tMin){ tShift = _tMin; }
  void SetTmax(double _tMax){ tMax = _tMax; }
  void SetMomentum(const Point &_momentum){momentum = make_unique<Point>(_momentum);}
  
  void UpdateOrigin(const Point &_origin);
  
  void SetVertex(const Point &_vertex){
    points[0] = make_shared<Point>(_vertex);
    vertex = make_unique<Point>(_vertex);
  }
  double chi2;
  bool increasing = false;
  
  void ReplacePoints(vector<shared_ptr<Point>> &_points){
    points = _points;
  }
  bool shouldRefit=false;
  
private:
  unique_ptr<Point> vertex;     ///< Decay point (beginning) of the helix
  Point origin;                 ///< Center of the helix
  vector<shared_ptr<Point>> points;   ///< Vector of points laying on the helix (withing thickness)
  
  double radius;                ///< Radius of the helix
  double slope;                 ///< Slope of the helix in Z direction
  double slopeAbs;              ///< Absolute value of the slope (to speed up the calculation)
  double tShift;  ///< Angle by which beginning of the helix is shifted due to the shift of its origin
  double tMax;    ///< Max angle (taking into account number of cycles
  double tStep;   ///< Step determining drawing precision
  
  int nRegularPoints = 0; ///< Number of points that are distributed regularly along Z axis
  int nPionPoints = 0;    ///< Number of points along the helix that are true pion hits
  
  unique_ptr<Point> momentum;   ///< Pion's momentum vector
  Track track;
  int charge;                ///< Charge of the particle (determines helix direction)
  
  Point GetClosestPoint(const Point &p) const;
  Point eventVertex;
  
  friend class HelixProcessor;
};

#endif /* Helix_hpp */
