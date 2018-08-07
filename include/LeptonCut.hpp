//
//  LeptonCut.hpp
//  xDisappTracks
//
//  Created by Jeremi Niedziela on 07/08/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#ifndef LeptonCut_hpp
#define LeptonCut_hpp

class LeptonCut {
public:
  enum ECut {
    kEmpty    = 1,
    kTightID  = 1 << 1, ///< require tight identification
    kPt20GeV  = 1 << 2, ///< pT ≥ 20 GeV
    kIsolated = 1 << 3, ///< isolation ≤ 0.15
    
  };
  
  LeptonCut(ECut cutType=kEmpty);
  ~LeptonCut();
  
  // setters
  inline void SetPtRange(int min, int max){minPt=min;maxPt=max;}
  inline void SetMaxIsolation(double max){maxIsolation = max;}
  inline void SetRequireTightID(bool id){requireTightID = id;}
  
  // getters
  inline int GetMinPt(){return minPt;}
  inline int GetMaxPt(){return maxPt;}
  inline bool RequiresTightID(){return requireTightID;}
  inline double GetMaxIsolation(){return maxIsolation;}
  
private:
  double minPt;         ///< min pT of the jet
  double maxPt;         ///< max pT of the jet
  bool requireTightID;  ///< should tight ID be required
  double maxIsolation;  ///< max isolation value
  
};

#endif /* LeptonCut_hpp */
