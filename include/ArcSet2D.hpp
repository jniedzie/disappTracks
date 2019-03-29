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

/// Contains a vector of circle segments along a track, where each segement may have different radius
/// and center. Stores points along the track as well. Provides additional information about the quality
/// of the track, such as how linearly the radii change from segment to segment.
class ArcSet2D {
public:
  /// Default constructor
  ArcSet2D();
  
  /// Copy constructors
  ArcSet2D(const ArcSet2D &a);
  
  /// Default destructor
  ~ArcSet2D();
  
  /// Prints basic information about this track
  void Print();
  
  /// Returns phi angle of the first track point relative to the seed center
  double GetOriginPhi();
  
  /// Returns number of segments along the track
  inline unsigned long GetNarcs(){return circles.size();}
  
  /// Returns current numeber of full spiral cycles
  inline int GetCycle(){return iCycle;}
  
  /// Increases cycles counter by 1
  inline void IncreaseCycle(){iCycle++;}
  
  /// Adds new segment to the track, calculating its range and increasing cycles counter if needed.
  /// It will also add last point on added circle to the list of track's points
  void AddCircle(const unique_ptr<Circle> &circle);
  
  /// Adds point to the track
  inline void AddPoint(shared_ptr<Point> p){points.push_back(p);}
  
  /// Adds vector of points to the track
  void AddPoints(vector<shared_ptr<Point>> p);
  
  /// Returns first point of the track
  inline shared_ptr<Point> GetOrigin(){return points[0];}
  
  /// Returns total number of points along the track
  inline unsigned long GetNpoints(){return points.size();}
  
  /// Returns all points along the track
  inline vector<shared_ptr<Point>> GetPoints(){return points;}
  
  /// Returns i-th point on the track
  inline shared_ptr<Point> GetPoint(uint i){return points[i];}
  
  /// Returns last point on the track
  inline shared_ptr<Point> GetLastPoint(){return points[points.size()-1];}
  
  /// Returns second to last point on the track
  inline shared_ptr<Point> GetSecondToLastPoint(){return points[points.size()-2];}
  
  // Returns last segment of the track
  inline unique_ptr<Circle> GetLastCircle(){return make_unique<Circle>(circles[circles.size()-1]);}
  
  /// Returns i-th segemnt of the track
  inline unique_ptr<Circle> GetCircle(uint i){return make_unique<Circle>(circles[i]);}
  
  /// Returns all segments of the track
  inline vector<unique_ptr<Circle>>& GetCircles(){return circles;}
  
  /// Returns vector or TArc objects corresponding to track's segments that can be plotted in TCanvas
  vector<TArc*> GetArcs();
  
  /// Returns chi2 of a linear fit to radius vs. segment index graph
  double GetRadiiSlopeChi2();
  
private:
  vector<unique_ptr<Circle>> circles; ///< Segments of the track (first one is the seed)
  vector<shared_ptr<Point>> points;   ///< Points along the track
  
  int iCycle; ///< Current number of spiral's cycles
  
  unique_ptr<CircleProcessor> circleProcessor;
  
  friend class ArcSetProcessor;
};

#endif /* ArcSet2D_hpp */
