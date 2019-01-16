//
//  Fitter.hpp
//
//  Created by Jeremi Niedziela on 14/12/2018.
//

#ifndef Fitter_hpp
#define Fitter_hpp

#include "Helpers.hpp"
#include "Point.hpp"
#include "PointsProcessor.hpp"
#include "Circle.hpp"
#include "Helix.hpp"
#include "FitterConfig.hpp"

class Fitter {
public:
  Fitter();
  ~Fitter();
  
  
  unique_ptr<Helix> GetBestFittingHelix(vector<Point> allSimplePoints,
                                               shared_ptr<FitterConfig> config,
                                               double trackTheta, double trackPhi,
                                               bool drawCircles=false);
  
private:
  vector<unique_ptr<Circle>> FitCirclesToPoints(vector<Point> allSimplePoints,
                                                int pxSign, int pySign, int charge,
                                                shared_ptr<FitterConfig> config,
                                                double trackTheta, double trackPhi);
  
  unique_ptr<PointsProcessor> pointsProcessor;
  
  void SetParameter(ROOT::Fit::Fitter *fitter, int i, string name, double start, double min, double max, bool fix=false);
  void FixParameter(ROOT::Fit::Fitter *fitter, int i, string name, double val);
};

#endif /* Fitter_hpp */
