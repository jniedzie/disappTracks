//
//  AcrSetProcessor.hpp
//
//  Created by Jeremi Niedziela on 26/03/2019.
//

#ifndef AcrSetProcessor_hpp
#define AcrSetProcessor_hpp

#include "Helpers.hpp"
#include "ArcSet2D.hpp"
#include "PointsProcessor.hpp"
#include "CircleProcessor.hpp"

class ArcSetProcessor;
extern ArcSetProcessor arcSetProcessor;

class ArcSetProcessor {
public:
  /// Default constructor
  ArcSetProcessor();
  
  /// Default destructor
  ~ArcSetProcessor();
  
  /// Creates a vector of ArcSet2D, each already containing a seed.
  /// \param circles Circles from which to build the seeds
  vector<unique_ptr<ArcSet2D>> BuildArcSetsFromCircles(const vector<unique_ptr<Circle>> &circles);
  
  /// Builds a vector of point triplets that that potentially could extend an existing ArcSet2D track
  /// \param arcSet Existing track for which candidate triplets should be found
  /// \param points Vector of all points in space from which only those compatible will be selected
  TripletsVector BuildTripletsCompatibleWithArcSet(const unique_ptr<ArcSet2D> &arcSet,
                                                   const vector<shared_ptr<Point>> &points);
  
  
  /// Returns vector of points that potenially could extend an existing ArcSet2D track
  /// \param arcSet Existing track for which candidate triplets should be found
  /// \param points Vector of all points in space from which only those compatible will be selected
  vector<shared_ptr<Point>> FindPossibleNextPoints(const unique_ptr<ArcSet2D> &arcSet,
                                                   const vector<shared_ptr<Point>> &points);
  
  /// Finds the best arc set based on how close to linear the radii change from segment to segment is
  unique_ptr<ArcSet2D> GetBestArcSet(const vector<unique_ptr<ArcSet2D>> &arcSets);
  
private:
  
  bool IsValidSeed(const unique_ptr<Circle> &circle);
};

#endif /* AcrSetProcessor_hpp */
