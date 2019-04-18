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

class HelixProcessor;
extern HelixProcessor helixProcessor;

class HelixProcessor {
public:
  /// Default contrustor
  HelixProcessor();
  
  /// Default destructor
  ~HelixProcessor();
  
  /// Checks if input and output helices are identical.
  /// \return Returns zero if identical, otherwise returns failure reason code
  vector<int> AreIdentical(const Helix &h1, const Helix &h2);
  
  /// Returns vector of points along helix trajectory that hit the tracker
  vector<shared_ptr<Point>> GetPointsHittingSilicon(const Helix &helix);
  
  /// Calculates number of regular points.
  /// Splits all points into lines along Z axis. For each line, checks all possible distances between points.
  /// For each possible distance, calculates number of regular points (for all points in the collection,
  /// within zRegularityTolarance) and finds a maximum number of such points.
  /// \param helix Helix object that will be modified with an updated number of regular points
  /// \param limit Stop calculation after reaching this number of regular points
  void CalculateNregularPoints(unique_ptr<Helix> &helix, int limit=inf);
  
  /// Generates random helix starting on the track (withing the decay region)
  /// \param track Track on which the helix will be generated
  /// \param pionHelix Returns generated helix
  void GetRandomPionHelix(const shared_ptr<Track> &track, Helix &pionHelix);
  
  
  /// Link class variables to branches of a specified tree
  /// \param tree Tree from which fitted helix parameters will be read
  void SetupBranchesForReading(TTree *tree);
  
  /// Returns a vector of fitted helices read from tree previously set with SetupBranchesForReading(..)
  vector<shared_ptr<Helix>> GetHelicesFromTree();
  
  /// Link class variables to branches of a specified tree
  /// \param tree Tree to which fitted helices parameters will be saved
  void SetupBranchesForWriting(TTree *tree);
  
  /// Writes all fitted helices in the vector to the tree previously set with SetupBranchesForWriting(...)
  void SaveHelicesToTree(vector<shared_ptr<Helix>> helices);
  
  double GetHelixToPointDistance(const unique_ptr<Helix> &helix, const shared_ptr<Point> &point);
  
  double GetChi2toPoints(const unique_ptr<Helix> &helix, const vector<shared_ptr<Point>> &points);
  
private:
  static const int maxNhelices = 1000;   ///< Maximum supported number of helices per event
  int nHelices;                          ///< Number of helices in the current tree entry
  
  map<string, float[maxNhelices]> arrayValuesFloat; ///< Float per-helix variables in the current entry
  map<string, int[maxNhelices]>   arrayValuesInt;   ///< Int per-helix variables in the current entry
  
  vector<string> arrayNamesFloat;     ///< Names or float per-track variables
  vector<string> arrayNamesInt;       ///< Names or int per-track variables
  
};

#endif /* HelixProcessor_hpp */
