//  Helix.hpp
//
//  Created by Jeremi Niedziela on 18/12/2018.

#ifndef Helix_hpp
#define Helix_hpp

#include "Helpers.hpp"
#include "Point.hpp"
#include "PointsProcessor.hpp"
#include "ConfigManager.hpp"
#include "Track.hpp"

struct HelixParams
{
  // Helix params given parametric equation:
  // x = x0 + (R0 - at) cos(t)
  // y = y0 + (R0 - at) sin(t)
  // z = z0 + (s0 - bt) t
  
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
  /// Constructor taking as an input origin and momentum vector (will be automatically shifted to begin in the origin point)
  /// \param _origin Helix origin point (e.g. chargino decay point)
  /// \param _momentum  Momentum of the particle that creates a helix
  /// \param _charge Charge of the particle (determines helix direction)
  Helix(const Point &_origin,
        const unique_ptr<Point> &_momentum,
        int _charge);
  
  Helix(const HelixParams &_params,
        const Point &_decayVertex,
        const Point &_origin,
        const vector<shared_ptr<Point>> &_points,
        const Track &_track);
  
  /// Copy constructor
  Helix(const Helix &h);
  
  /// Assignent operator
  Helix& operator=(const Helix &h);
  
  /// Comparison operator
  bool operator==(const Helix &h);
  
  /// Prints basic information about the helix
  void Print();
  
  /// Adds point to the collection of helix's points, sets correct t for this point and updates tMax
  void AddPoint(const shared_ptr<Point> &point);
  
  /// Sets new origin, recalculates t params for all point on helix and updates tShift and tMax
  void UpdateOrigin(const Point &_origin);
  
  // Setters
  inline void   SetCharge(int val)            { charge = val; }
  inline void   SetTmin(double val)           { tShift = val; }
  inline void   SetTmax(double val)           { tMax = val; }
  inline void   SetMomentum(const Point &val) { momentum = make_unique<Point>(val); }
  inline void   SetPoints(vector<shared_ptr<Point>> &val){ points = val; }
  void          SetVertex(const Point &_vertex);
  inline void   SetShouldRefit(bool val) { shouldRefit = val; }
  inline void   SetChi2(double val) { chi2 = val; }
  inline void   SetIsFinished(bool val) { isFinished = val; }
  inline void   SetParams(HelixParams val) { helixParams = val; }
  inline void   SetIsPreviousHitMissing(bool val){ isPreviousHitMissing = val; }
  inline void   SetFirstTurningPointIndex(int val){ firstTurningPointIndex = val; }
  
  /// Increases number of missing hits and missing hits in a row (based on isPreviousHitMissing variable, which
  /// should be set by the user).
  void   IncreaseMissingHits();
  
  
  // Getters
  inline int                  GetCharge()   const {return charge; }
  inline uint64_t             GetUniqueID() const {return uniqueID; }
  
  inline double               GetTmin()     const {return tShift;}
  inline double               GetTmax()     const {return tMax;}
  inline double               GetTstep()    const {return tStep;}
  
  double                      GetRadius(double t) const;
  inline double               GetRadiusFactor()   const {return helixParams.a; }
  double                      GetSlope(double t)  const;
  inline double               GetSlopeFactor()    const {return helixParams.b; }
  
  inline Point                GetOrigin()   const {return origin;}
  inline unique_ptr<Point>    GetVertex()   const {return make_unique<Point>(*vertex);}
  inline unique_ptr<Point>    GetMomentum() const {return make_unique<Point>(*momentum);}
  
  vector<shared_ptr<Point>>   GetPoints()   const {return points;}
  vector<shared_ptr<Point>>   GetLastPoints() const;
  vector<shared_ptr<Point>>   GetSecontToLastPoints() const;
  
  inline uint                 GetNpoints()  const {return (uint)(points.size());}
  
  inline int                  GetNmissingHits() const {return nMissingHits;}
  inline int                  GetNmissingHitsInRow() const {return nMissingHitsInRow;}
  inline bool                 IsPreviousHitMissing() const {return isPreviousHitMissing;}
  
  inline int   GetFirstTurningPointIndex(){ return firstTurningPointIndex; }
  
  double                      GetNcycles()  const;
  
  inline bool GetShouldRefit() const { return shouldRefit; }
  inline double GetChi2() const { return chi2; }
  inline bool GetIsFinished() const { return isFinished; }
  inline bool IsIncreasing() const { return firstTurningPointIndex<0; }
  
private:
  unique_ptr<Point> vertex;         ///< Decay point (beginning) of the helix
  Point origin;                     ///< Center of the helix
  vector<shared_ptr<Point>> points; ///< Vector of points laying on the helix
  
  double tShift;  ///< Angle by which beginning of the helix is shifted due to the shift of its origin
  double tMax;    ///< Max angle (taking into account number of cycles
  double tStep;   ///< Step determining drawing precision
  
  unique_ptr<Point> momentum;   ///< Pion's momentum vector
  Track track;
  int charge;                ///< Charge of the particle (determines helix direction)
  
  int iCycles;
  bool isFinished;
  uint64_t seedID;
  uint64_t uniqueID;
  HelixParams helixParams;
  double chi2;
  bool shouldRefit;
  
  int firstTurningPointIndex;
  
  int nMissingHits;
  int nMissingHitsInRow;
  bool isPreviousHitMissing;
  
  friend class HelixProcessor;
};

#endif /* Helix_hpp */
