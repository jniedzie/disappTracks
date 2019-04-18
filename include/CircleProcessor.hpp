//  CircleProcessor.hpp
//
//  Created by Jeremi Niedziela on 26/03/2019.

#ifndef CircleProcessor_hpp
#define CircleProcessor_hpp

#include "Helpers.hpp"
#include "Circle.hpp"
#include "Track.hpp"
#include "PointsProcessor.hpp"

class CircleProcessor;
extern CircleProcessor circleProcessor;

class CircleProcessor {
public:
  /// Default constructor
  CircleProcessor();
  
  
  /// Default destructor
  ~CircleProcessor();
  
  /// Removes circles which radii are within circleThickness parameter in the config
  /// If two circles are similar, the one with longer arc will be removed
  void RemoveSimilarCircles(vector<unique_ptr<Circle>> &circles);
  
  /// Takes a vector of triplets of points and returns a vector of circles matching these triplets
  vector<unique_ptr<Circle>> BuildCirclesFromPoints(const TripletsVector &points);
  
  /// For each triplet pair builds circles that go through the triplet points and stores only
  /// the circle with greater radius.
  /// \param radiusRange Circles that have radius outside of this range will be skipped
  vector<unique_ptr<Circle>> BuildCirclesFromTripletPairs(const TripletPairsVector &triplets,
                                                          range<double> radiusRange);
  
  /// Creates a circle from a parameters: par = { L, px, py }, where L is a distance from primary vertex
  /// to the decay point, px and py are momentum vector components.
  /// \param par Array of parameters { L, px, py }
  /// \param vertex Primary vertex of the event
  /// \param track Track from which the circle originates
  unique_ptr<Circle> BuildCircleFromParams(const double *par,
                                           const Point &vertex,
                                           const Track &track);
  
  /// Finds a circle that's the most compatible with the given circle
  /// \param circles Vector of circles to scan
  /// \param theCircle The circle to which the returned one will be similar
  /// \param alphaVector Will be filled with angular differences between the circle and all provided circles
  unique_ptr<Circle> GetMostCompatibleCircle(const vector<unique_ptr<Circle>> &circles,
                                             const unique_ptr<Circle> &theCircle,
                                             vector<double> &alphaVector);
  
  /// Copies privded circle, changing its phi range to the provided one
  unique_ptr<Circle> CopyCircleAddingRange(const unique_ptr<Circle> &circle,
                                           const range<double> &phiRange);
  
  /// Returns a new circle that begins at the last point of the provided circle, is tangent to it at this
  /// point and passes through the provided point
  unique_ptr<Circle> GetParallelCircle(const unique_ptr<Circle> &circle,
                                       const shared_ptr<Point> &point);
  
  /// Returns a circle passing through all three points of the specified triplet
  unique_ptr<Circle>  GetCircleFromTriplet(const PointsTriplet &triplet);
  
private:
  
  bool IsPerpendicular(const Point &p1, const Point &p2,const Point &p3);
  unique_ptr<Circle> CalcCircle(const Point &p1, const Point &p2, const Point &p3);

};

#endif /* CircleProcessor_hpp */
