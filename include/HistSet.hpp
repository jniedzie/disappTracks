//
//  HistSet.hpp
//
//  Created by Jeremi Niedziela on 16/07/2018.
//

#ifndef HistSet_hpp
#define HistSet_hpp

#include "EventSet.hpp"

#include <TH1D.h>
#include <TCanvas.h>

#include <vector>

/// HistSet contains three vectors of histograms, for signal, background and data. It provides methods
/// that allow to automatically fill histograms from EventSet, according to settings in the Helpers.hpp file.
class HistSet
{
public:
  /// Default constructor
  /// \param _var Specifies which variable (e.g. Jet pT) the histograms will store
  HistSet(EVar _var);
  
  /// Default destructor
  ~HistSet();
  
  /// Fills all histograms from provided EventSet
  /// \param events EventSet from which histograms will be filled
  void FillFromEvents(shared_ptr<EventSet> events);
  
  /// Draws histograms for signal, background and data
  /// \param canvas Pointer to canvas in which histograms should be plotted
  /// \param pad Index of the pad in the canvas in which to draw the histograms
  void Draw(TCanvas *canvas, int pad);
  
  /// Draws per-layer plots. This will automatically create a new canvas for each variable.
  void DrawPerLayer();
  
private:
  vector<TH1D*> signal;     ///< Vector of histograms for signal events
  vector<TH1D*> background; ///< Vector of histograms for background events
  vector<TH1D*> data;       ///< Vector of histograms for data events
  
  vector<vector<TH1D*>> signalPerLayer;     ///< Array of per-layer histograms for signal [signalType][layer]
  vector<vector<TH1D*>> backgroundPerLayer; ///< Array of per-layer histograms for background [signalType][layer]
  vector<vector<TH1D*>> dataPerLayer;       ///< Array of per-layer histograms for data [signalType][layer]
  
  EVar var;                 ///< For which variable the histograms should be filled
  const char* title;        ///< Title of the variable
  int nBins;                ///< Number of bins for this variable
  double min;               ///< X axis minimum for this variable
  double max;               ///< X axis maximum for this variable
  
  /// Fills per-layer plots with data from privided EventSet
  /// \param events EventSet from which histograms will be filled
  void FillFromEventsPerLayer(shared_ptr<EventSet> events);
  
  /// Returns TLegend object with default position, size and header
  /// \return Returns already prepared TLegend object
  TLegend* GetLegend();
  
  /// Tells whether given histogram should be normalized or not
  bool ShouldNormalize();
  
  /// Fills given histogram type with data from the EventSet.
  /// \param hist Histogram which will be filled
  /// \param events EventSet from which histogram will be filled
  /// \param dataType Type of the data for which the histogram will be filled (signal/background/data)
  /// \param setIter Data set for which to fill the histogram
  /// \param iDetId Detector ID for per-layer histograms
  void Fill(TH1D* hist, shared_ptr<EventSet> events,
            EventSet::EDataType dataType, int setIter,
            int iDetId=-1);
};

#endif /* HistSet_hpp */
