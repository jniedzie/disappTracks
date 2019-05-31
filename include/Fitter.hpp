//
//  Fitter.hpp
//
//  Created by Jeremi Niedziela on 14/12/2018.
//

#ifndef Fitter_hpp
#define Fitter_hpp

#include "Helpers.hpp"
#include "PointsProcessor.hpp"
#include "HelixProcessor.hpp"
#include "ConfigManager.hpp"
#include "Track.hpp"

/// Provides a method to fit a helix to a collection of points.
class Fitter {
public:
  /// Default constructor
  Fitter();
  
  /// Default destructor
  ~Fitter();
  
  vector<Helix> FitHelices(const vector<shared_ptr<Point>> &_points,
                           const Track &_track,
                           const Point &_eventVertex);
  
private:
  vector<shared_ptr<Point>> points;
  Track track;
  Point eventVertex;
  
  ///
  vector<Helix> GetSeeds(vector<vector<shared_ptr<Point>>> pointsByLayer);
  
  ///
  unique_ptr<Helix> FitSeed(const vector<shared_ptr<Point>> &points, int charge);
  
  ///
  void ExtendSeeds(vector<Helix> &helices,
                   const vector<vector<shared_ptr<Point>>> &pointsByLayer);
  
  ///
  bool MergeHelices(vector<Helix> &helices);
  
  ///
  void RefitHelix(Helix &helix);
  
  ///
  ROOT::Fit::Fitter* GetSeedFitter(const vector<shared_ptr<Point>> &points);
  
  /// Sets name, limits and starting value of a ROOT fitter's parameter
  void SetParameter(ROOT::Fit::Fitter *fitter, int i, string name, double start, double min, double max, bool fix=false);
  
  /// Sets and fixes a value for a ROOT fitter's parameter
  void FixParameter(ROOT::Fit::Fitter *fitter, int i, string name, double val);
};

#endif /* Fitter_hpp */
