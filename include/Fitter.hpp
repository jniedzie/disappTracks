//
//  Fitter.hpp
//
//  Created by Jeremi Niedziela on 14/12/2018.
//

#ifndef Fitter_hpp
#define Fitter_hpp

#include "Helpers.hpp"
#include "PointsProcessor.hpp"
#include "CircleProcessor.hpp"
#include "ArcSetProcessor.hpp"
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
  
  /// Returns the best fitting helix for a set of points, given the origin along a track with known eta and phi
  /// \points Vector of points to which the helix will be fitted
  /// \trackTheta Theta angle of the track
  /// \trackPhi Phi angle of the track
  /// \drawCircles Optionally, graph with all candidate circles can be plotted
  unique_ptr<Helix> GetBestFittingHelix(const vector<shared_ptr<Point>> _points,
                                        const Track &_track,
                                        const Point &_vertex,
                                        bool drawCircles=false);
  
  vector<Helix> FitHelix(const vector<shared_ptr<Point>> &_points,
                         const Track &_track,
                         const Point &_eventVertex);
  
  HelixParams FitHelixParams(const vector<shared_ptr<Point>> &points, const Point &nextPoint,
                             const Point &origin, const Track &track,
                             const Point &eventVertex, EHelixParams iParam);
  
  vector<Helix> FitHelix2(const vector<shared_ptr<Point>> &_points,
                          const Track &_track,
                          const Point &_eventVertex);
  
  unique_ptr<Helix> FitSeed(const vector<shared_ptr<Point>> &points,
                            const Track &track,
                            const Point &eventVertex);
  
private:
  
  vector<shared_ptr<Point>> points;
  Track track;
  Point eventVertex;
  
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
  
  /// Returns fitter that fits circle in terms of px, py and L, based on track parameters and settings in config
  ROOT::Fit::Fitter* GetCirclesFitter(int pxSign, int pySign);
  
  /// Fits seeds to doublets of points and a track, setting first point of each triplet to the most
  /// probable value. Trashes solutions that are above given chi2 threshold.
  /// \param pointTriplets Input points - first one in each triplet doesn't matter, as it will be set by this method
  /// \param pxSign px component direction
  /// \param pySign py component direction
  /// \param chi2threshold Max allowed chi2
  vector<unique_ptr<Circle>> FitCirclesAndAdjustFirstPoint(TripletsVector &pointTriplets,
                                                           int pxSign, int pySign,
                                                           double chi2threshold);
  
  void PlotSeeds(       int iPad, const vector<unique_ptr<ArcSet2D>> &potentialPionTracks);
  void PlotTracks(      int iPad, const vector<unique_ptr<ArcSet2D>> &potentialPionTracks);
  void PlotGoodTracks(  int iPad, const vector<unique_ptr<ArcSet2D>> &potentialPionTracks);
  void PlotRadiiChi2(   int iPad, const vector<unique_ptr<ArcSet2D>> &potentialPionTracks);
  void PlotRadiiVsIter( int iPad, const vector<unique_ptr<ArcSet2D>> &potentialPionTracks);
  
  void PlotClusters(    int iPad, const vector<shared_ptr<Point>> &filteredPoints);
  
  void PlotBestTrack(   int iPad, const unique_ptr<ArcSet2D> &pionTrack);
  void PlotRadiiAngles( int iPad, const vector<double> &alphaVector);
  
  TGraph* GetDecayGraph();
  TCanvas *c1;
  
  ///
  ROOT::Fit::Fitter* GetHelixParamsFitter(range<double> rangeL);
  ROOT::Fit::Fitter* GetSeedFitter(range<double> rangeL);
};

#endif /* Fitter_hpp */
