//
//  CircleProcessor.hpp
//
//  Created by Jeremi Niedziela on 26/03/2019.
//

#ifndef CircleProcessor_hpp
#define CircleProcessor_hpp

#include "Helpers.hpp"
#include "Circle.hpp"


class CircleProcessor {
public:
  /// Default constructor
  CircleProcessor();
  
  
  /// Default destructor
  ~CircleProcessor();
  
  /// Removes circles which radii are within circleThickness parameter in the config
  void RemoveSimilarCircles(vector<unique_ptr<Circle>> &circles);
  
  /// Takes a vector of triplets of points and returns a vector of circles matching these triplets
  vector<unique_ptr<Circle>> BuildCirclesFromPoints(const vector<vector<shared_ptr<Point>>> &points);
  
private:
  
  
};

#endif /* CircleProcessor_hpp */
