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

#include "CircleProcessor.hpp"

class ArcSet2D {
public:
  /// Default constructor
  ArcSet2D();
  
  /// Copy constructor
  ArcSet2D(const ArcSet2D &a);
  ArcSet2D(const unique_ptr<ArcSet2D> &a);
  
  /// Default destructor
  ~ArcSet2D();
  
  void Print();
  
  double GetOriginPhi();
  
  inline unsigned long GetNarcs(){return circles.size();}
  
  inline int GetCycle(){return iCycle;}
  
  inline void IncreaseCycle(){iCycle++;}
  
  void AddCircle(const unique_ptr<Circle> &circle);
  
  inline void AddPoint(shared_ptr<Point> p){points.push_back(p);}
  
  void AddPoints(vector<shared_ptr<Point>> p);
  
  inline shared_ptr<Point> GetOrigin(){return points[0];}
  
  inline vector<shared_ptr<Point>> GetPoints(){return points;}
  
  inline shared_ptr<Point> GetPoint(uint i){return points[i];}
  
  inline shared_ptr<Point> GetLastPoint(){return points[points.size()-1];}
  
  inline shared_ptr<Point> GetSecondToLastPoint(){return points[points.size()-2];}
  
  inline unique_ptr<Circle> GetLastCircle(){return make_unique<Circle>(circles[circles.size()-1]);}
  
  inline unique_ptr<Circle> GetCircle(uint i){return make_unique<Circle>(circles[i]);}
  
  vector<TArc*> GetArcs();
  
private:
  vector<unique_ptr<Circle>> circles;     ///< Segments of the track (first one is the seed)
  vector<shared_ptr<Point>> points;       ///< Points along the track
  
  int iCycle;
  
  unique_ptr<CircleProcessor> circleProcessor;
  
  friend class ArcSetProcessor;
};

#endif /* ArcSet2D_hpp */
