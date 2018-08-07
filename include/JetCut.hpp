//
//  JetCut.hpp
//  xDisappTracks
//
//  Created by Jeremi Niedziela on 20/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#ifndef JetCut_hpp
#define JetCut_hpp

class JetCut {
public:
  enum ECut {
    kEmpty    = 1,
    kPt100GeV = 1 << 1, ///< pT ≥ 100 GeV
    kPt200GeV = 1 << 2, ///< pT ≥ 200 GeV
    kEta2p4   = 1 << 3, ///< |eta| < 2.4
    kChHEF0p1 = 1 << 4, ///< charged hadron energy fraction ≥ 0.1
    kNeHEF0p8 = 1 << 5, ///< neutral hadron energy fraction ≤ 0.8
  };
  
  JetCut(ECut cutType=kEmpty);
  ~JetCut();
  
  inline double GetMinPt(){return minPt;}
  inline double GetMaxPt(){return maxPt;}
  inline double GetMaxEta(){return maxEta;}
  inline double GetMinChargedHadronEnergyFraction(){return minChargedHadronEnergyFraction;}
  inline double GetMaxNeutralHadronEnergyFraction(){return maxNeutralHadronEnergyFraction;}
  
  inline void SetPtRange(double min, double max){minPt=min;maxPt=max;}
  inline void SetMaxEta(double max){maxEta = max;}
  inline void SetMinChargedHadronEnergyFraction(double min){minChargedHadronEnergyFraction = min;}
  inline void SetMaxNeutralHadronEnergyFraction(double max){maxNeutralHadronEnergyFraction = max;}
  
private:
  double minPt;   ///< min pT of the jet
  double maxPt;   ///< max pT of the jet
  double maxEta;  ///< max pseudorapidity
  double minChargedHadronEnergyFraction; ///< min charged hadron energy fraction
  double maxNeutralHadronEnergyFraction; ///< max neutral hadron energy fraction
  
};

#endif /* JetCut_hpp */
