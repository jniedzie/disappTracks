//  JetProcessor.hpp
//
//  Created by Jeremi Niedziela on 31/01/2019.

#ifndef JetProcessor_hpp
#define JetProcessor_hpp

#include "Helpers.hpp"
#include "Jet.hpp"

class JetProcessor;
extern JetProcessor jetProcessor;

/// Class description
class JetProcessor {
public:
  /// Default constructor
  JetProcessor();
  
  /// Default destructor
  ~JetProcessor();
  
  /// Check if jet passes selection criteria
  /// \param jet Jet to which cuts should be applied
  /// \param cut Jet selection criteria to be checked
  bool IsPassingCut(const shared_ptr<Jet> jet, const JetCut &cut);
  
  /// Link class variables to branches of a specified tree
  /// \param tree Tree from which jet parameters will be read
  void SetupBranchesForReading(TTree *tree);
  
  /// Returns a vector of jet with parameters read from tree previously set with SetupBranches(..)
  vector<shared_ptr<Jet>> GetJetsFromTree();
  
  /// Link class variables to branches of a specified tree
  /// \param tree Tree to which jets parameters will be saved
  void SetupBranchesForWriting(TTree *tree);
  
  /// Writes all jets in the vector to the tree previously set with SetupBranchesForWriting(...)
  void SaveJetsToTree(vector<shared_ptr<Jet>> jets);
  
private:
  static const int maxNjets = 1000;   ///< Maximum supported number of jets per event
  int nJets;                          ///< Number of jets in the current tree entry
  int nJetsFwd;                       ///< Number of forward jets in the current tree entry
  
  map<string, float[maxNjets]>  arrayValuesFloat;  ///< Float per-jet variables in the current entry
  
  vector<string> arrayNamesFloat;     ///< Names or float per-track variables
};

#endif /* JetProcessor_hpp */
