//  MonitorsManager.hpp
//
//  Created by Jeremi Niedziela on 18/01/2019.

#ifndef MonitorsManager_hpp
#define MonitorsManager_hpp

#include "Helpers.hpp"
#include "HelixProcessor.hpp"
#include "Track.hpp"
#include "ConfigManager.hpp"

class MonitorsManager;
extern MonitorsManager monitorsManager;

/// Creates, populates, plots and saves hisograms monitoring fitting performance.
class MonitorsManager {
public:
  /// Default constructor
  /// It will automatically create all the monitors
  MonitorsManager();
  
  /// Default destructor
  ~MonitorsManager();
  
  /// Creates new canvas, draws all monitoring histograms and saves them to the output file
  void PlotAndSaveMonitors();
  
  /// Status of the fit
  enum EFitStatus {
    kFail,        ///< Didn't manage to fit any helix
    kSuccess,     ///< Fitted a helix, but their properites are not similar to the true one
    kFullSuccess  ///< Fitted a helix very similar to the true one
  };
  
  /// Fills defined monitors for a given set of true and fitted helix.
  void FillMonitors(const unique_ptr<Helix> &fittedHelix,
                    const Helix &trueHelix,
                    const shared_ptr<Track> &track);
  
  /// Checks and returns status of the fit given the true and the fitted helix
  EFitStatus GetFittingStatus(const Helix &fittedHelix, const Helix &trueHelix);
  
private:

  map<string, TH1D*> monitors1D;                     ///< 1D Monitors
  map<string, TH2D*> monitors2D;                     ///< 2D Monitors
  map<string, pair<TH1D*, TH1D*>> fractionMonitors;  ///< Monitors that contain numerator and denominator
};



#endif /* MonitorsManager_hpp */
