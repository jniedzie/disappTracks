//
//  Fitter.hpp
//
//  Created by Jeremi Niedziela on 14/12/2018.
//

#ifndef Fitter_hpp
#define Fitter_hpp

#include "Helpers.hpp"
#include "PointsProcessor.hpp"
#include "Circle.hpp"
#include "HelixProcessor.hpp"
#include "FitterConfig.hpp"
#include "Track.hpp"

/// Provides a method to fit a helix to a collection of points.
class Fitter {
public:
  /// Default constructor
  Fitter(shared_ptr<FitterConfig> _config);
  
  /// Default destructor
  ~Fitter();
  
  /// Returns the best fitting helix for a set of points, given the origin along a track with known eta and phi
  /// \points Vector of points to which the helix will be fitted
  /// \trackTheta Theta angle of the track
  /// \trackPhi Phi angle of the track
  /// \drawCircles Optionally, graph with all candidate circles can be plotted
  unique_ptr<Helix> GetBestFittingHelix(const shared_ptr<vector<Point>> _points,
                                        const shared_ptr<Track> _track,
                                        bool drawCircles=false);
  
private:
  shared_ptr<FitterConfig> config;
  unique_ptr<PointsProcessor> pointsProcessor;
  unique_ptr<HelixProcessor> helixProcessor;
  
  shared_ptr<vector<Point>> points;
  shared_ptr<Track> track;
  
  
  /// Finds circles fitting points, matching all the other criteria specified in the config and for
  /// a specific signs of X and Y coordinates of the momentum and charge.
  /// Points, track theta and phi must be set before calling this method.
  /// \pxSign Momentum X coordinate sign to be tested
  /// \pySign Momentum Y coortinate sign to be tested
  /// \charge Charge of the pion
  vector<unique_ptr<Circle>> FitCirclesToPoints(int pxSign, int pySign);
  
  /// Finds all possible circle candidates for previously set collection of points and track parameters.
  vector<unique_ptr<Circle>> GetAllCirclesForPoints();
  
  unique_ptr<Helix> GetHelixFromCircle(const unique_ptr<Circle> &circle, double pz, int charge);
  
  /// Checks whether provided helix has better properties than the one specified as arguments. If so,
  /// it will not only return true, but also update maxNregularPoints and maxFractionRegularPoints
  bool IsHelixBetterThanBefore(const unique_ptr<Helix> &helix,
                               int &maxNregularPoints,
                               double &maxFractionRegularPoints);
  
  /// Sets name, limits and starting value of a ROOT fitter's parameter
  void SetParameter(ROOT::Fit::Fitter *fitter, int i, string name, double start, double min, double max, bool fix=false);
  
  /// Sets and fixes a value for a ROOT fitter's parameter
  void FixParameter(ROOT::Fit::Fitter *fitter, int i, string name, double val);
};

#endif /* Fitter_hpp */
