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

/// Display class creates an event display window and allows to visualize vector of points, a helix of a complete event
class Display {
public:
  /// Default constructor
  Display();
  
  /// Default destructor
  ~Display();
  
  /// Visualizes vector of points
  /// \param points Points to be visualized
  /// \param options Map containing plotting options
  void DrawSimplePoints(const vector<shared_ptr<Point>> points, map<string,any> options);
  
  /// Visualizes a helix
  /// \param helix Helix to be visualized
  /// \param options Map containing plotting options
  void DrawHelix(const Helix &helix, const map<string,any> options);
  
  void DrawShrinkingHelix(const Helix &helix, const map<string,any> options);
  
  /// Visualizes entire event
  /// \param event event to be visualized
  /// \param options Map containing plotting options
  void DrawEvent(const shared_ptr<Event> &event, const map<string,any> options);
  
private:
  double scale             = 0.1;           ///< Scale factor that will be applied to all visualized objects
  double jetConeRadius     = scale * 0.4;   ///< Determines size of jets' cones
  double metRadius         = scale * 2000;  ///< Distance at which the MET box will be drawn
  double metBoxSize        = scale * 30;    ///< Size of the MET box in radial coordinate
  double metBoxAngularSize = 0.1;           ///< Size of the MET box in the angular coordinate
  
  bool showUnderflowBins   = false;         ///< Should underflow dE/dx hits be shown by defdault
  bool showOverflowBins    = true;          ///< Should overflow dE/dx hits be shown by defdault
  
  int geomTransparency     = 90;            ///< Transparency of geometry (something between 30 - 100)
  
  //                                  Underflow                                   Overflow
  const vector<int> dedxBinColors = { kGray, kBlue, kCyan, kGreen, kYellow, kRed, kMagenta };

  /// Returns a point set that can be filled with XYZ points and their values and then visualized
  TEvePointSetArray* PreparePointsEventDisplay(map<string,any> options);
  
  /// Draws MET box at given phi-theta coordinates
  void DrawMET(double metPhi, double metTheta);
  
  void AddStripCluster(TEveElementList *stripClusters, const shared_ptr<Point> &point, map<string,any> options);
};

#endif /* Display_hpp */
