//
//  Jet.hpp
//  disappTracksTarget
//
//  Created by Jeremi Niedziela on 17/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#ifndef Jet_hpp
#define Jet_hpp

class Jet{
public:
  Jet();
  ~Jet();

  void Print();
  
  bool IsPassingCut(/*JetCut *cut*/);
private:
  
};

#endif /* Jet_hpp */
