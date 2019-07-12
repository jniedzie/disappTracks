//  PointsProcessor.hpp
//
//  Created by Jeremi Niedziela on 16/01/2019.

#ifndef PointsProcessor_hpp
#define PointsProcessor_hpp

#include "Helpers.hpp"
#include "Point.hpp"
#include "ConfigManager.hpp"
#include "Track.hpp"

typedef vector<shared_ptr<Point>> Points; ///< Vector of Point objects

class PointsProcessor;
extern PointsProcessor pointsProcessor;

/// PointsProcessor provides methods performing operations on Point objects
class PointsProcessor {
public:
  /// Default constructor
  PointsProcessor();
  
  /// Default destructor
  ~PointsProcessor();
  
  /// Returns distance between this and another point
  double distance(const Point &p1,const Point &p2) const;
  double distance(shared_ptr<Point> p1, shared_ptr<Point> p2) const;
  
  /// Returns distance between this and another point in the transverse plane
  double distanceXY(Point p1, Point p2) const;
  
  /// Sorts provided points by layer 
  vector<Points> SortByLayer(const Points &points);
  
  /// Sorts points by disk (from -N to -1 and from 1 to N)
  vector<Points> SortByDisk(const Points &points);
  
  /// Returns an angle between (p2-p1) vector and (p1-p0) vector shifted to start in p1
  double GetPointingAngle(const Point &p0, const Point &p1, const Point &p2);
  
  /// Returns an angle between (p2-p1) vector and (p1-p0) vector in the transverse plane
  double GetPointingAngleXY(const Point &p0, const Point &p1, const Point &p2);
  
  /// Returns an angle between Z coordinates and the transverse plane of (p2-p1) vector
  /// and (p1-p0) vector
  double GetPointingAngleTZ(const Point &p0, const Point &p1, const Point &p2);
  
  /// Returns point laying on the track (assuming straight line) at distance L, taking
  /// into account position of the primary event vertex
  Point GetPointOnTrack(double L, const Track &track, const Point &eventVertex);
  
  
  vector<Points> RegroupNerbyPoints(const Points &points);
  
  double GetTforPoint(const Point &point, const Point &origin, int charge);
  
  void SetPointsLayers(Points &points);
  void SetPointsDisks(Points &points);
  
  bool IsPhiGood(const Points &lastPoints,
                 const Points &secondToLastPoints,
                 const shared_ptr<Point> &point,
                 int charge);
  
  bool IsZgood(const Points &lastPoints,
               const shared_ptr<Point> &point);
  
  bool IsTgood(const vector<double> &lastPointsT, double pointT);
  
  Points GetGoodMiddleSeedHits(const Points &middlePoints,
                                                  const Point &trackMidPoint,
                                                  const Point &eventVertex,
                                                  int charge);
  
  Points GetGoodLastSeedHits(const Points &lastPoints,
                                                const Points &goodMiddlePoints,
                                                const Point &trackMidPoint,
                                                int charge);
  
  /// Finds minimum distance between helix and point
  /// \param params Helix parameters: {R0, a, s0, b, L, x0, y0, z0}
  /// \param tMin t param of the decay vertex
  /// \param point Point to which the distane will be calculated
  /// \param alpha Rotation angle of the cluster in transverse plane
  /// \param charge Assumed track's charge
  double GetMinHelixToPointDistance(const double *params, double tMin,
                                    const Point &point, double alpha, int charge);
  
  /// Structs for sorting points
  struct ComparePointByZ{
    bool operator() (const shared_ptr<Point> &p1, const shared_ptr<Point> &p2){
      return (p1->GetZ() < p2->GetZ());
    }
  };
  
  struct ComparePointByX{
    bool operator() (const shared_ptr<Point> &p1, const shared_ptr<Point> &p2){
      return (p1->GetX() < p2->GetX());
    }
  };
  
  struct ComparePointByY{
    bool operator() (const shared_ptr<Point> &p1, const shared_ptr<Point> &p2){
      return (p1->GetY() < p2->GetY());
    }
  };
  
  struct ComparePointByLayer{
    bool operator() (const shared_ptr<Point> &p1, const shared_ptr<Point> &p2){
      return (p1->GetLayer() < p2->GetLayer());
    }
  };
  
  struct ComparePointByLayerInverted{
    bool operator() (const shared_ptr<Point> &p1, const shared_ptr<Point> &p2){
      return (p1->GetLayer() > p2->GetLayer());
    }
  };
  
  struct ComparePointByTime{
    bool operator() (const shared_ptr<Point> &p1, const shared_ptr<Point> &p2){
      return (p1->GetTime() < p2->GetTime());
    }
  };
  
private:
  
  bool IsGoodMiddleSeedHit(const Point &point,
                           const Point &trackMidPoint,
                           const Point &eventVertex,
                           int charge);
  
  bool IsGoodLastSeedHit(const Point &point,
                         const Point &middlePoint,
                         const Point &trackMidPoint,
                         int charge);
  
};

#endif /* PointsProcessor_hpp */
