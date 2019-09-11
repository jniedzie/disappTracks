//
//  Logger.hpp
//  xDisappTracks
//
//  Created by Jeremi Niedziela on 11/09/2019.
//  Copyright © 2019 Jeremi Niedziela. All rights reserved.
//

#ifndef Logger_hpp
#define Logger_hpp

#include "Helpers.hpp"
#include "ConfigManager.hpp"

class Log
{
public:
  Log(int _level) : level(_level) {}
  
  template <class T>
  Log &operator<<(const T &v){
    if(config.params["verbosity_level"] >= level) cout << v;
    return *this;
  }
  
  ~Log() {}
  
private:
  int level;
};

#endif /* Logger_hpp */
