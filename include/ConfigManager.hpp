//
//  FitterConfig.hpp
//
//  Copyright © 2019 Jeremi Niedziela. All rights reserved.
//

#ifndef FitterConfig_hpp
#define FitterConfig_hpp

#include "Helpers.hpp"

struct ConfigManager;
extern ConfigManager config;

/**
 Wrapper on a config file that provides access to options from the code.
 Reads a config file in markdown format. For all options that are not specified in the config
 default values will be used.
 */
struct ConfigManager {
  /// Default constructor
  /// \param path Path to the config file
  ConfigManager(string path="");
  
  map<string, double> params; ///< All bool, int and double parameters from the config file
  
  string category;              ///< Name of the analysis category
  string secondaryCategory;     ///< Name of the secondary analysis category (e.g. Wmunu)
  string outputPath;            ///< Output file path
  
  vector<bool> runBackground; ///< Should run given backgorund sample (by EBackground enum)
  vector<bool> runSignal;     ///< Should run given signal sample (by ESignal enum)
  vector<bool> runData;       ///< Should run given data sample (by EData enum)
};

#endif /* FitterConfig_hpp */
