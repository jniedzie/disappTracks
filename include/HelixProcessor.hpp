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
  
private:
  static const int maxNhelices = 1000;   ///< Maximum supported number of helices per event
  int nHelices;                          ///< Number of helices in the current tree entry
  
  map<string, float[maxNhelices]> arrayValuesFloat; ///< Float per-helix variables in the current entry
  map<string, int[maxNhelices]>   arrayValuesInt;   ///< Int per-helix variables in the current entry
  
  vector<string> arrayNamesFloat;     ///< Names or float per-track variables
  vector<string> arrayNamesInt;       ///< Names or int per-track variables
  
};

#endif /* HelixProcessor_hpp */
