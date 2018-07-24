//
//  EventCut.hpp
//
//  Created by Jeremi Niedziela on 23/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#ifndef EventCut_hpp
#define EventCut_hpp

class EventCut {
public:
  enum ECut {
    kEmpty,
    kMet100GeV
  };
  
  EventCut(ECut cutType=kEmpty);
  ~EventCut();
  
  inline double GetMinMetPt(){return minMetPt;}
  
  inline void SetMinMetPt(double min){minMetPt=min;}
  
private:
  double minMetPt;      ///< min number of dedx clusters along the track
  
};

#endif /* EventCut_hpp */
