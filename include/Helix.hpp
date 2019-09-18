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
  /// \param _decayVertex Helix starting point (e.g. chargino decay vertex)
  /// \param _momentum  Momentum of the particle that creates a helix
  /// \param _charge Charge of the particle (determines helix direction)
  Helix(const Point &_decayVertex,
        const Point &_momentum,
        int _charge);
  
  Helix(const HelixParams &_params,
        const Point &_decayVertex,
        const Point &_origin,
        int _charge);
  
  Helix(const ROOT::Fit::FitResult &result,
        const Track &track,
        const Point &eventVertex,
        int _charge);
  
  /// Copy constructor
  Helix(const Helix &h);
  
  /// Assignent operator
  Helix& operator=(const Helix &h);
  
  /// Comparison operator
  bool operator==(const Helix &h);
  
  /// Prints basic information about the helix
  void Print();
  
  /// Replaces all points on helix and corresponding t params. Then sorts points by
  /// t param. Will also set collection last points (but not second to last).
  void SetPointsAndSortByT(const Points &points, const vector<double> &_pointsT);
  
  /// Will shift current last points to second to last points, save second to last in
  /// third to last. Then puts all provided points in last points collection
  /// (and all points collection too). T params will be calculated for all points
  /// and tMax will be updated.
  void SetLastPoints(Points _points);
  
  /// Increases number of missing hits and missing hits in a row (based on
  /// isPreviousHitMissing variable, which should be set by the user).
  void IncreaseMissingHits();
  
  /// Removes all last points, shifts second to last to last and third to last to
  /// second to last. If last point was missing, the missing points counters will
  /// be decreased.
  void RemoveLastPoints();
  
  /// Adds point to the collection of helix's points, sets correct t for this point
  /// and updates tMax
  void AddPoint(const shared_ptr<Point> &point, vector<size_t> *_lastPointIndices = nullptr);
  
  /// Sets new origin, recalculates t params for all point on helix and updates
  /// tShift and tMax
  void UpdateOrigin(const Point &_origin);
  
  //------------------------------------------------//
  //                    Getters                     //
  //------------------------------------------------//
  
  // Basic helix info
  inline int      GetCharge()         const {return charge;}
  inline uint64_t GetUniqueID()       const {return uniqueID;}
  inline uint64_t GetSeedID()         const {return seedID;}
  inline Point    GetMomentum()       const {return momentum;}
  inline double   GetTstep()          const {return tStep;}
  inline double   GetTmin()           const {return pointsT.front();}
         double   GetTmax(vector<size_t> *_lastPointIndices = nullptr) const;
         double   GetNcycles()        const;
         size_t   GetNlayers()        const;
         double   GetRadius(double t) const;
         double   GetSlope(double t)  const;
  inline double   GetRadiusFactor()   const {return helixParams.a;}
  inline double   GetSlopeFactor()    const {return helixParams.b;}
  inline bool     IsIncreasing()      const {return firstTurningPointIndex<0;}
  inline double   GetChi2()           const {return chi2;}
  inline bool     GetShouldRefit()    const {return shouldRefit;}
  inline bool     GetIsFinished()     const {return isFinished;}
  
  // Helix points
  inline Point              GetOrigin()                 const {return origin;}
  inline shared_ptr<Point>  GetVertex()                 const {return points.front();}
  inline Points             GetPoints()                 const {return points;}
  inline size_t             GetNpoints()                const {return points.size();}
  inline vector<double>     GetPointsT()                const {return pointsT;}
  inline double             GetPointT(size_t index)     const {return pointsT[index];}
  inline Points             GetLastPoints()             const {return lastPoints;}
  inline Points             GetSecontToLastPoints()     const {return secondToLastPoints;}
         vector<size_t>     GetLastPointsIndices()      const;
  inline int                GetFirstTurningPointIndex() const {return firstTurningPointIndex;}
  
  // Missing hits
  inline int  GetNmissingHits() 	    const {return nMissingHits;}
  inline int  GetNmissingHitsInRow()  const {return nMissingHitsInRow;}
  inline bool IsPreviousHitMissing()  const {return isPreviousHitMissing;}
  
  //------------------------------------------------//
  //                    Setters                     //
  //------------------------------------------------//
  
  inline void SetCharge(int val)                 {charge = val;}
  inline void SetMomentum(const Point &val)      {momentum = val;}
  inline void SetShouldRefit(bool val)           {shouldRefit = val;}
  inline void SetChi2(double val)                {chi2 = val;}
  inline void SetIsFinished(bool val)            {isFinished = val;}
  inline void SetParams(HelixParams val)         {helixParams = val;}
  inline void SetIsPreviousHitMissing(bool val)  {isPreviousHitMissing = val;}
  inline void SetFirstTurningPointIndex(int val) {firstTurningPointIndex = val;}
  inline void SetPoints(const Points &_points)   {points = _points;};
  inline void SetVertex(const Point &_vertex)    {points[0] = make_shared<Point>(_vertex);}
  
private:
  Point           origin;             ///< Center of the helix
  Points          points;             ///< Vector of points laying on the helix
  vector<double>  pointsT;            ///< T values of points
  Points          lastPoints;         ///< Last points on helix
  Points          secondToLastPoints; ///< Second to last points on helix
  Points          thirdToLastPoints;  ///< Second to last points on helix
  
  int firstTurningPointIndex; ///< Index of the first point after helix turns back to the same layer
  
  uint64_t  seedID;   ///< Unique ID of the helix's seed (can be the same as for other helices)
  uint64_t  uniqueID; ///< Unique ID of this helix
  Point     momentum; ///< Pion's momentum vector
  int       charge;   ///< Charge of the particle (determines helix direction)
  HelixParams helixParams; ///< Parameters defining radius, slope and their changes with time
  double tStep;       ///< Step determining drawing precision
  
  bool   isFinished;  ///< Is extending of this helix over?
  bool   shouldRefit; ///< Does this helix require refitting?
  double chi2;        ///< Chi2 of the fit to helix points
  
  int  nMissingHits;          ///< Total number of missing hits along the trajectory
  int  nMissingHitsInRow;     ///< Max number of consecutive missing hits
  bool isPreviousHitMissing;  ///< Is the point on helix a missing hit?
  int nRecLayers;             ///< Total number of reconstructed layers (should be set when fitting is done)
  int nRecHits;               ///< Total number of reconstructed hits (should be set when fitting is done)
  
  
  /// Sorts collection of all points by t param in increasing order (or decreasing
  /// if inverted==true).
  void SortPointsByT(bool inverted);
  
  friend class HelixProcessor;
};

#endif /* Helix_hpp */
