//
//  JetCut.hpp
//
//  Created by Jeremi Niedziela on 20/07/2018.
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
    kFwdEta4p7= 1 << 6, ///< |eta| < 4.7 for forward jets
    kPt30GeV  = 1 << 7, ///< pT ≥ 30 GeV
  };
  
  JetCut(ECut cutType=kEmpty);
  ~JetCut();
  
  inline double GetMinPt(){return minPt;}
  inline double GetMaxPt(){return maxPt;}
  inline double GetMaxEta(){return maxEta;}
  inline double GetMaxEtaFwd(){return maxEtaFwd;}
  inline double GetMinChargedHadronEnergyFraction(){return minChargedHadronEnergyFraction;}
  inline double GetMaxNeutralHadronEnergyFraction(){return maxNeutralHadronEnergyFraction;}
  inline double GetMinTrackDeltaR(){return minTrackDeltaR;}
  
  inline void SetPtRange(double min, double max){minPt=min;maxPt=max;}
  inline void SetMaxEta(double max){maxEta = max;}
  inline void SetMaxEtaFwd(double max){maxEtaFwd = max;}
  inline void SetMinChargedHadronEnergyFraction(double min){minChargedHadronEnergyFraction = min;}
  inline void SetMaxNeutralHadronEnergyFraction(double max){maxNeutralHadronEnergyFraction = max;}
  inline void SetMinTrackDeltaR(double min){minTrackDeltaR = min;}
  
private:
  double minPt;   ///< min pT of the jet
  double maxPt;   ///< max pT of the jet
  double maxEta;  ///< max pseudorapidity
  double minChargedHadronEnergyFraction; ///< min charged hadron energy fraction
  double maxNeutralHadronEnergyFraction; ///< max neutral hadron energy fraction
  double maxEtaFwd;       ///< max pseudorapidity for forward jets
  double minTrackDeltaR;  ///< min allowed separation with any of the tracks
  
};

#endif /* JetCut_hpp */
