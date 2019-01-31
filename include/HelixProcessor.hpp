//
//  HelixProcessor.hpp
//
//  Created by Jeremi Niedziela on 21/01/2019.
//

#ifndef HelixProcessor_hpp
#define HelixProcessor_hpp

#include "Helpers.hpp"
#include "Helix.hpp"
#include "ConfigManager.hpp"
#include "PointsProcessor.hpp"
#include "Track.hpp"

class HelixProcessor {
public:
  /// Default contrustor
  HelixProcessor();
  
  /// Default destructor
  ~HelixProcessor();
  
  /// Checks if input and output helices are identical.
  /// \return Returns zero if identical, otherwise returns failure reason code
  vector<int> AreIdentical(const unique_ptr<Helix> &h1, const unique_ptr<Helix> &h2);
  
  /// Returns vector of points along helix trajectory that hit the tracker
  shared_ptr<vector<Point>> GetPointsHittingSilicon(const unique_ptr<Helix> &helix);
  
  /// Calculates number of regular points.
  /// Splits all points into lines along Z axis. For each line, checks all possible distances between points.
  /// For each possible distance, calculates number of regular points (for all points in the collection,
  /// within zRegularityTolarance) and finds a maximum number of such points.
  /// \helix Helix object that will be modified with an updated number of regular points
  /// \limit Stop calculation after reaching this number of regular points
  void CalculateNregularPoints(unique_ptr<Helix> &helix, int limit=inf);
  
  unique_ptr<Helix> GetRandomPionHelix(const shared_ptr<Track> &track);
  
  
  /// Link class variables to branches of a specified tree
  /// \param tree Tree from which helix parameters will be read
  void SetupBranches(TTree *tree);
  
  /// Returns a vector of helices with parameters read from tree previously set with SetupBranches(..)
  vector<shared_ptr<Helix>> GetHelicesFromTree();
  
private:
  static const int maxNhelices = 1000;   ///< Maximum supported number of helices per event
  int nHelices;                          ///< Number of helices in the current tree entry
  
  map<string, float[maxNhelices]> arrayValuesFloat; ///< Float per-helix variables in the current entry
  map<string, int[maxNhelices]>   arrayValuesInt;   ///< Int per-helix variables in the current entry
  
  unique_ptr<PointsProcessor> pointsProcessor;
};

#endif /* HelixProcessor_hpp */
