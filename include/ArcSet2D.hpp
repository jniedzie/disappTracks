//
//  ArcSet2D.hpp
//  xDisappTracks
//
//  Created by Jeremi Niedziela on 18/03/2019.
//  Copyright © 2019 Jeremi Niedziela. All rights reserved.
//

#ifndef ArcSet2D_hpp
#define ArcSet2D_hpp

#include "Helpers.hpp"

#include "Circle.hpp"

class ArcSet2D {
public:
  ArcSet2D();
  ~ArcSet2D();
  
  void Print();
  
  void AddCircle(const unique_ptr<Circle> &circle, range<double> range);
  
  inline void AddPoint(shared_ptr<Point> p){points.push_back(p);}
  
  void AddPoints(vector<shared_ptr<Point>> p);
  
  inline shared_ptr<Point> GetOrigin(){return points[0];}
  
  inline shared_ptr<Point> GetPoint(uint i){return points[i];}
  
  inline unique_ptr<Circle> GetCircle(uint i){return make_unique<Circle>(circles[i]);}
  
  vector<TArc*> GetArcs();
  
private:
  vector<unique_ptr<Circle>> circles;
  vector<range<double>> circlesRanges;
  vector<shared_ptr<Point>> points; ///< Points along the track
  
};

#endif /* ArcSet2D_hpp */
