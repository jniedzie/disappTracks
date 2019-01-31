//
//  JetProcessor.hpp
//
//  Created by Jeremi Niedziela on 31/01/2019.
//

#ifndef JetProcessor_hpp
#define JetProcessor_hpp

#include "Helpers.hpp"
#include "Jet.hpp"

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
  bool IsPassingCut(const shared_ptr<Jet> jet,
                    const unique_ptr<JetCut> &cut);
  
  /// Link class variables to branches of a specified tree
  /// \param tree Tree from which jet parameters will be read
  void SetupBranches(TTree *tree);
  
  /// Returns a vector of jet with parameters read from tree previously set with SetupBranches(..)
  vector<shared_ptr<Jet>> GetJetsFromTree();
  
private:
  static const int maxNjets = 1000;   ///< Maximum supported number of jets per event
  int nJets;                          ///< Number of jets in the current tree entry
  int nJetsFwd;                       ///< Number of forward jets in the current tree entry
  
  map<string, float[maxNjets]>  arrayValuesFloat;  ///< Float per-track variables in the current entry
};

#endif /* JetProcessor_hpp */
