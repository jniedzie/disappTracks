//
//  Fitter.hpp
//
//  Created by Jeremi Niedziela on 14/12/2018.
//

#ifndef Fitter_hpp
#define Fitter_hpp

#include "Helpers.hpp"
#include "Point.hpp"
#include "Circle.hpp"
#include "Helix.hpp"
#include "FitterConfig.hpp"

class Fitter {
public:
  Fitter(int _nPar);
  ~Fitter();
  
  template<typename Func>
  void SetFitFunction(const Func &function){
    fitFunction = ROOT::Math::Functor(function, nPar);
    double pStart[nPar];
    fitter->SetFCN(fitFunction, pStart);
  }
  
  void SetParameter(int i, string name, double start, double min, double max, bool fix=false);
  
  void FixParameter(int i, string name, double val);
  
  bool RunFitting();
  
  const ROOT::Fit::FitResult& GetResult();
  
  static vector<unique_ptr<Circle>> FitCirclesToPoints(vector<Point> allSimplePoints,
                                                       int pxSign, int pySign, int charge,
                                                       shared_ptr<FitterConfig> config,
                                                       double trackTheta, double trackPhi);
  
  static unique_ptr<Helix> GetBestFittingHelix(vector<Point> allSimplePoints,
                                               shared_ptr<FitterConfig> config,
                                               double trackTheta, double trackPhi);
  
private:
  int nPar;
  
  ROOT::Math::Functor fitFunction;
  ROOT::Fit::Fitter *fitter;
};

#endif /* Fitter_hpp */
