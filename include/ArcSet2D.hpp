//
//  ArcSet2D.hpp
//  xDisappTracks
//
//  Created by Jeremi Niedziela on 18/03/2019.
//  Copyright © 2019 Jeremi Niedziela. All rights reserved.
//

#ifndef ArcSet2D_hpp
#define ArcSet2D_hpp

#include "Helpers.hpp"

#include "Circle.hpp"

class ArcSet2D {
public:
  ArcSet2D();
  ~ArcSet2D();
  
  void AddCircle(const unique_ptr<Circle> &circle, range<double> range);
  
  vector<TArc*> GetArcs();
  
private:
  vector<unique_ptr<Circle>> circles;
  vector<range<double>> circlesRanges;
  
};

#endif /* ArcSet2D_hpp */
