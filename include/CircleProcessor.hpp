//
//  CircleProcessor.hpp
//
//  Created by Jeremi Niedziela on 26/03/2019.
//

#ifndef CircleProcessor_hpp
#define CircleProcessor_hpp

#include "Helpers.hpp"
#include "Circle.hpp"
#include "Track.hpp"
#include "PointsProcessor.hpp"


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
  
  /// Creates a circle from a parameters: par = { L, px, py }, where L is a distance from primary vertex
  /// to the decay point, px and py are momentum vector components.
  /// \param par Array of parameters { L, px, py }
  /// \param vertex Primary vertex of the event
  /// \param track Track from which the circle originates
  unique_ptr<Circle> BuildCircleFromParams(const double *par,
                                           const unique_ptr<Point> &vertex,
                                           const shared_ptr<Track> &track);
  
  /// Finds a circle that's the most compatible with the given circle
  /// \param circles Vector of circles to scan
  /// \param theCircle The circle to which the returned one will be similar
  /// \param alphaVector Will be filled with angular differences between the circle and all provided circles
  unique_ptr<Circle> GetMostCompatibleCircle(const vector<unique_ptr<Circle>> &circles,
                                             const unique_ptr<Circle> &theCircle,
                                             vector<double> &alphaVector);
  
  unique_ptr<Circle> CopyCircleAddingRange(const unique_ptr<Circle> &circle,
                                           const range<double> &phiRange);
  
private:
  unique_ptr<PointsProcessor> pointsProcessor;
  
};

#endif /* CircleProcessor_hpp */
