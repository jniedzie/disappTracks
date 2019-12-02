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
  
  /// Returns a vector of helices fitting provided points
  /// \param _points Collection of input points
  /// \param _track Track along which the helix starts
  /// \param _eventVertex Primary vertex of the event
  Helices FitHelices(const Points &_points,
                     const Track &_track,
                     const Point &_eventVertex,
                     int _nTrackerLayers=-1);
  
private:
  Points points;                    ///< collection of all points in the event
  vector<Points> pointsByLayer;     ///< all points in event grouped by layer index
  vector<Points> pointsByDisk;      ///< all points in event grouped by disk index
  
  Track track;                      ///< track from which helix should start
  Point eventVertex;                ///< primary vertex of the event
  
  int nTrackLayers;                 ///< number of layers before helix should start
  int charge;                       ///< assumed charge of the helix
  
  int nDegreesOfFreedom = 6;        ///< number of free fit parameters
  Points fittingPoints; ///< collection of points to which we are fitting at the moment
  function<double(const double*)> chi2Function; ///< chi2 function
  
  // Initial parameters and their limits:
  double startL, minL, maxL;
  
  Helices PerformFittingCycle();
  
  /// Checks parameters of all combinations of points in layers close to the decay point
  /// and returns vector of 3-layer helices constructed from them
  Helices GetSeeds();
  
  /// Fits helix of given charge to the collection of points provided
  /// \return Resulting helix may be a nullptr if something went wrong (e.g. parameters
  /// of helix fitting provided points were outside of assumed limits.
  unique_ptr<Helix> FitSeed(const Points &middleHits,
                            const Points &lastHits);
  
  /// Attempts to extend provided seeds to following layers. If not possible, tries to turn
  /// back to the same layer. If that's also not possible, assigns missing hit if still
  /// allowed by the limits.
  void ExtendSeeds(Helices &helices);
  
  /// Merges similar helices as long as there's nothing left to merge
  void MergeHelices(Helices &helices);
  
  /// Performs single merging operation (accorging to merging_max_different_point
  /// and candidate_min_n_points parameters). All helices that don't meet merging
  /// conditions are left intact.
  bool LinkAndMergeHelices(Helices &helices);
  
  /// Refits helix params using points assigned to this helix
  void RefitHelix(Helix &helix);
  
  /// Removes helices that are below track_min_n_points threshold
  void RemoveShortHelices(Helices &helices);
  
  /// Creates a fitter for seeds, with the best guess of the initial parameters
  /// \return nullptr if parameters were outside of limits, but requested only good params
  ROOT::Fit::Fitter* GetSeedFitter();
  
  /// Sets name, limits and starting value of a ROOT fitter's parameter
  void SetParameter(ROOT::Fit::Fitter *fitter, int i, string name, double start, double min, double max, bool fix=false);
  
  /// Sets and fixes a value for a ROOT fitter's parameter
  void FixParameter(ROOT::Fit::Fitter *fitter, int i, string name, double val);
  
  /// Returns the most probable X0 and Y0 limits and starting values for given decay point
  /// and assumed charge
  void GetXYranges(const Point &trackPoint,
                   double &startX0, double &minX0, double &maxX0,
                   double &startY0, double &minY0, double &maxY0);
  
  /// Sets L limits and starting value based on the curent number of tracker layers
  void InitLparams();
  
  
  vector<set<int>> GetLayersAndDisks();
  
  pair<Points, Points> GetSeedPoints();
  
  unique_ptr<Helix> CreateHelixFromFitResult(const ROOT::Fit::FitResult &result);
  
  vector<set<int>> GetNextPointLayersAndDisks(const Helix &helix);
  
  Points GetPossibleNextPoints(const Helix &helix);
  Points GetPossibleTurningBackPoints(const Helix &helix);
  
  Points GetGoodNextPoints(const Helix &helix,
                           const Points &points,
                           const set<int> &nextPointLayers);
  
  Points GetGoodTurningBackPoints(const Helix &helix,
                                  const Points &points);
  
  unique_ptr<Helix> TryToExtendHelix(const Helix &helix, const Points &points, bool turnsBack);
};

#endif /* Fitter_hpp */
