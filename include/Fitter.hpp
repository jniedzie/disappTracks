//
//  Fitter.hpp
//
//  Created by Jeremi Niedziela on 14/12/2018.
//

#ifndef Fitter_hpp
#define Fitter_hpp

#include "Helpers.hpp"

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
  
private:
  int nPar;
  
  ROOT::Math::Functor fitFunction;
  ROOT::Fit::Fitter *fitter;
};

#endif /* Fitter_hpp */
