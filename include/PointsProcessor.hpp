//
//  PointsProcessor.hpp
//
//  Created by Jeremi Niedziela on 16/01/2019.
//

#ifndef PointsProcessor_hpp
#define PointsProcessor_hpp

#include "Helpers.hpp"
#include "Point.hpp"
#include "ConfigManager.hpp"

/// PointsProcessor provides methods performing operations on Point objects
class PointsProcessor {
public:
  /// Default constructor
  PointsProcessor(const shared_ptr<ConfigManager> &_config);
  
  /// Default destructor
  ~PointsProcessor();
  
  /// Returns distance between this and another point
  double distance(Point p1, Point p2) const;
  
  /// Returns distance between this and another point in the transverse plane
  double distanceXY(Point p1, Point p2) const;
  
  /// Returns squared distance between this and another point in the transverse plane
  /// (should be much faster than a true distance)
  double distanceXYsquared(Point p1, Point p2) const;
  
  /// Separates vector of points into groups with the same XY position (within tolerance)
  vector<vector<Point>> SplitPointsIntoLines(const shared_ptr<vector<Point>> points, double tolerance) const;
  
  /// Returns a vector filled with random points in the pixel barrel
  /// \param nPoints Number of points that will be generated
  shared_ptr<vector<Point>> GetRandomPoints(int nPoints) const;
  
private:
  shared_ptr<ConfigManager> config;
  
};

#endif /* PointsProcessor_hpp */
