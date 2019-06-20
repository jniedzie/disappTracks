//  PointsProcessor.hpp
//
//  Created by Jeremi Niedziela on 16/01/2019.

#ifndef PointsProcessor_hpp
#define PointsProcessor_hpp

#include "Helpers.hpp"
#include "Point.hpp"
#include "ConfigManager.hpp"
#include "Track.hpp"

typedef vector<shared_ptr<Point>> PointsTriplet; ///< Vector containing three points
typedef pair<shared_ptr<Point>, shared_ptr<Point>> PointsPair;
typedef vector<PointsTriplet> TripletsVector;    ///< Vector of triplets of points
typedef vector<pair<PointsTriplet, PointsTriplet>> TripletPairsVector; ///< Vector of pairs of point triplets

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
  
  double distanceWithUncertainZ(shared_ptr<Point> p1, shared_ptr<Point> p2, double zTolerance) const;
  
  /// Returns distance between this and another point in the transverse plane
  double distanceXY(Point p1, Point p2) const;
  
  /// Returns squared distance between this and another point in the transverse plane
  /// (should be much faster than a true distance)
  double distanceXYsquared(Point p1, Point p2) const;
  
  /// Separates vector of points into groups with the same XY position (within tolerance)
  vector<vector<Point>> SplitPointsIntoLines(const vector<shared_ptr<Point>> &points,
                                             double tolerance) const;
  
  /// Returns a vector filled with random points in the pixel barrel
  /// \param nPoints Number of points that will be generated
  vector<shared_ptr<Point>> GetRandomPoints(int nPoints) const;
  
  /// Returns a new vector containing points that are not closer to each other than given threshold
  /// \param points Input points vector
  /// \param minPointsSeparation Minimal allowed distance between two points
  vector<shared_ptr<Point>> FilterNearbyPoints(const vector<shared_ptr<Point>> &points,
                                               double minPointsSeparation);
  
  /// Creates all possible combinations of input points with a dummy first triplet point
  TripletsVector BuildPointTriplets(const vector<shared_ptr<Point>> &points);
  
  /// Creates all possible combinations of two input points
  vector<PointsPair> BuildPointPairs(const vector<shared_ptr<Point>> &points);
  
  /// Creates a vector of pairs of point triples. In each pair first one starts with provided originMin
  /// point, the second one with the originMax. Second and third element of each pair are all possible
  /// combinations of the input points provided.
  TripletPairsVector BuildPointTripletPairs(const vector<shared_ptr<Point>> &points,
                                            const shared_ptr<Point> &originMin,
                                            const shared_ptr<Point> &originMax);
  
  /// Sorts provided points by layer 
  vector<vector<shared_ptr<Point>>> SortByLayer(const vector<shared_ptr<Point>> &points);
  
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
  
  
  vector<vector<shared_ptr<Point>>> RegroupNerbyPoints(const vector<shared_ptr<Point>> &points,
                                                       double threshold);
  
  double GetTforPoint(Point &point, const Point &origin, int charge);
  
  void SetPointsLayers(vector<shared_ptr<Point>> &points);
  
  void SetPointsT(vector<shared_ptr<Point>> &points, const Point &origin, int charge);
  
  void SetPointsT(vector<shared_ptr<Point>> &points, const Point &origin, int charge,
                  shared_ptr<Point> &lastPointBeforeTurning);
  
  bool IsPhiGood(const vector<shared_ptr<Point>> &lastPoints,
                 const vector<shared_ptr<Point>> &secondToLastPoints,
                 const shared_ptr<Point> &point,
                 int charge);
  
  bool IsZgood(const vector<shared_ptr<Point>> &lastPoints,
               const shared_ptr<Point> &point);
  
  bool IsTgood(const vector<shared_ptr<Point>> &lastPoints,
               const shared_ptr<Point> &point);
  
  /// Structs for sorting point
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
  
  struct ComparePointByT{
    bool operator() (const shared_ptr<Point> &p1, const shared_ptr<Point> &p2){
      return (p1->GetT() < p2->GetT());
    }
  };
  
};

#endif /* PointsProcessor_hpp */
