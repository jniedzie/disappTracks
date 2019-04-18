//  LeptonProcessor.hpp
//
//  Created by Jeremi Niedziela on 31/01/2019.

#ifndef LeptonProcessor_hpp
#define LeptonProcessor_hpp

#include "Helpers.hpp"
#include "Lepton.hpp"

class LeptonProcessor;
extern LeptonProcessor leptonProcessor;

/// Class description
class LeptonProcessor {
public:
  /// Default constructor
  LeptonProcessor();
  
  /// Default destructor
  ~LeptonProcessor();
  
  /// Check if lepton passes selection criteria
  /// \param lepton Lepton to which cuts should be applied
  /// \param cut Lepton selection criteria to be checked
  bool IsPassingCut(const shared_ptr<Lepton> lepton, const LeptonCut &cut);
  
  /// Link class variables to branches of a specified tree
  /// \param tree Tree from which Lepton parameters will be read
  void SetupBranchesForReading(TTree *tree);
  
  /// Returns a vector of Lepton with parameters read from tree previously set with SetupBranches(..)
  vector<shared_ptr<Lepton>> GetLeptonsFromTree();
  
  /// Link class variables to branches of a specified tree
  /// \param tree Tree to which leptons parameters will be saved
  void SetupBranchesForWriting(TTree *tree);
  
  /// Writes all leptons in the vector to the tree previously set with SetupBranchesForWriting(...)
  void SaveLeptonsToTree(vector<shared_ptr<Lepton>> leptons);
  
private:
  static const int maxNleptons = 1000;   ///< Maximum supported number of Leptons per event
  int nLeptons;                          ///< Number of Leptons in the current tree entry
  
  map<string, float[maxNleptons]>  arrayValuesFloat;  ///< Float per-Lepton variables in the current entry
  map<string, int[maxNleptons]>    arrayValuesInt;    ///< Int per-Lepton variables in the current entry
  
  vector<string> arrayNamesFloat;     ///< Names or float per-track variables
  vector<string> arrayNamesInt;       ///< Names or int per-track variables
};


#endif /* LeptonProcessor_hpp */
