//
//  Display.hpp
//
//  Created by Jeremi Niedziela on 14/12/2018.
//

#ifndef Display_hpp
#define Display_hpp

#include "Helpers.hpp"
#include "Event.hpp"
#include "Jet.hpp"
#include "Track.hpp"

#include "Helix.hpp"
#include "Circle.hpp"
#include "Point.hpp"

class Display {
public:
  Display();
  ~Display();
  
  
  void DrawSimplePoints(vector<Point> points, const map<string,any> options);
  
  /// Draws a helix in given t parameter range
  /// \param helix Object of type Helix
  void DrawHelix(const unique_ptr<Helix> &helix, const map<string,any> options);
  
  void DrawEvent(shared_ptr<Event> event, const map<string,any> options);
  
private:
  
  double scale = 0.1;
  
  bool showUnderflowBins = false;
  bool showOverflowBins = true;
  
  bool showGeometry = false;
  
  double jetConeRadius = scale*0.4;
  
  double metRadius = scale * 2000;
  double metBoxSize = scale * 30;
  double metBoxAngularSize = 0.1;
  
  int geomTransparency = 90; // 30 - 100
  
  //     Underflow                                       Overflow
  const vector<int> dedxBinColors = { kGray,     kBlue, kCyan, kGreen, kYellow, kRed, kMagenta };

  TEvePointSetArray* PreparePointsEventDisplay(map<string,any> options);
  
  void DrawMET(double metPhi, double metTheta);
  
};

#endif /* Display_hpp */
