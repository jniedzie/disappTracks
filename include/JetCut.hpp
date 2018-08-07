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
  };
  
  JetCut(ECut cutType=kEmpty);
  ~JetCut();
  
  inline int GetMinPt(){return minPt;}
  inline int GetMaxPt(){return maxPt;}
  
  inline void SetPtRange(int min, int max){minPt=min;maxPt=max;}
  
private:
  double minPt;   ///< min pT of the jet
  double maxPt;   ///< max pT of the jet
  
};

#endif /* JetCut_hpp */
