//
//  Helpers.hpp
//
//  Created by Jeremi Niedziela on 13/06/2018.
//

#ifndef Helpers_h
#define Helpers_h

#pragma clang diagnostic push // save the current state
#pragma clang diagnostic ignored "-Wdocumentation" // turn ROOT's warnings
#pragma clang diagnostic ignored "-Wconversion"

// include ROOT headers
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TLorentzVector.h>
#include <THStack.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TFitter.h>
#include <TMath.h>
#include <TGraph2D.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <Fit/Fitter.h>
#include <Math/Functor.h>
#include <TEllipse.h>
#include <TArc.h>
#include <Fit/BinData.h>
#include <TSystem.h>
#include <TEveManager.h>
#include <TEveScene.h>
#include <TEvePointSet.h>
#include <TEveJetCone.h>
#include <TEveBox.h>
#include <TApplication.h>
#include <TGeoShape.h>
#include <TGeoTube.h>
#include <TEveGeoShape.h>
#include <TF3.h>
#include <TH3F.h>
#include <TLine.h>
#include <TEnv.h>
#include <TTreeReaderArray.h>
#include <TLorentzVector.h>

#pragma clang diagnostic pop // restores the saved state for diagnostics

#include <vector>
#include <iostream>
#include <map>
#include <numeric>
#include <algorithm>
#include <memory>
#include <utility>
#include <unordered_set>
#include <locale>
#include <variant>
#include <any>
#include <sstream>
#include <string>

using namespace std;

#define inf 99999999

//template<typename T, typename... Args>
//std::unique_ptr<T> make_unique(Args&&... args) {
//  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
//}

namespace xtracks {

const bool lxplus = false;

enum EHelixParams
{
  kMinS = 0,
  kMaxS,
  kMinR,
  kMaxR,
  kNhelixParams
};


enum EDataType{
  kBackground,
  kSignal,
  kData
};

// Plotting style
const double fillOpacity = 1.0;
const int fillStyleBack = 1000;
const int fillStyleSignal = 1000;
const int fillStyleData = 1000;

const vector<int> signalMarkers = {
  20, // wino m=300 cτ=3
  20, // wino m=300 cτ=10
  20, // wino m=300 cτ=30
  21, // wino m=500 cτ=10
  21, // wino m=500 cτ=20
  22, // wino m=650 cτ=10
  22, // wino m=650 cτ=20
  23, // wino m=800 cτ=10
  23, // wino m=800 cτ=20
  29, // wino m=1000 cτ=10
  29, // wino m=1000 cτ=20
};

const vector<vector<int>> backColors = {
  {230, 25 , 75 },  // QCD
  {60 , 180, 75 },  // Z->μμ + jets
  {0  , 130, 200},  // tops
  {245, 130, 48 },  // VV
  {145, 30 , 180},  // W->μν + jets
  {70 , 240, 240},  // Z->νν + jets
};

const vector<vector<int>> signalColors = {
  {128, 128, 0  },  // wino m=300 cτ=3
  {170, 110, 40 },  // wino m=300 cτ=10
  {128, 0  , 0  },  // wino m=300 cτ=30
  {170, 110, 40 },  // wino m=500 cτ=10
  {128, 0  , 0  },  // wino m=500 cτ=20
  {170, 110, 40 },  // wino m=650 cτ=10
  {128, 0  , 0  },  // wino m=650 cτ=20
  {170, 110, 40 },  // wino m=800 cτ=10
  {128, 0  , 0  },  // wino m=800 cτ=20
  {170, 110, 40 },  // wino m=1000 cτ=10
  {128, 0  , 0  },  // wino m=1000 cτ=20
  {170, 110, 40 },  //
  {128, 0  , 0  },  //
  {170, 110, 40 },  //
  {128, 0  , 0  },  //
};

// {140, 50 , 230} {170, 200, 195} {240, 50 , 100} {255, 225, 25 } {0  , 128, 128} {230, 190, 255} {170, 100, 195} {240, 50 , 230}

const vector<vector<int>> dataColors = {
  {200 , 10, 10},  // single electron (2017B)
  {200 , 10, 10},  // single electron (2018A)
  {200 , 10, 10},  // single muon CR (2018A)
};

enum EBackground{
  kQCD = 0,
  kZmumuJets,
  kTT,
  kVV,
  kWmunuJets,
  kZnunuJets
};

constexpr initializer_list<EBackground> backgrounds = {kQCD, kZmumuJets, kTT, kVV, kWmunuJets, kZnunuJets};

enum ESignal{
  kWino_M_300_cTau_3,
  kWino_M_300_cTau_10,
  kWino_M_300_cTau_30,
  kWino_M_500_cTau_10,
  kWino_M_500_cTau_20,
  kWino_M_650_cTau_10,
  kWino_M_650_cTau_20,
  kWino_M_800_cTau_10,
  kWino_M_800_cTau_20,
  kWino_M_1000_cTau_10,
  kWino_M_1000_cTau_20,
  kTaggerSignalNoPU,
  kTaggerBackgroundNoPU,
  kTaggerSignalWithPU,
  kTaggerBackgroundWithPU,
  kChargino300_1,
  kChargino300_10,
  kChargino400_1,
  kChargino500_1,
  kChargino500_10,
  kChargino600_10,
  kChargino700_10,
  kChargino700_30,
  kChargino800_10,
  kChargino800_30,
  kChargino900_30,
};

constexpr initializer_list<ESignal> signals = {
 kWino_M_300_cTau_3,
  kWino_M_300_cTau_10,
  kWino_M_300_cTau_30,
  kWino_M_500_cTau_10,
  kWino_M_500_cTau_20,
  kWino_M_650_cTau_10,
  kWino_M_650_cTau_20,
  kWino_M_800_cTau_10,
  kWino_M_800_cTau_20,
  kWino_M_1000_cTau_10,
  kWino_M_1000_cTau_20,
  kTaggerSignalNoPU,
  kTaggerBackgroundNoPU,
  kTaggerSignalWithPU,
  kTaggerBackgroundWithPU,
  kChargino300_1,
  kChargino300_10,
  kChargino400_1,
  kChargino500_1,
  kChargino500_10,
  kChargino600_10,
  kChargino700_10,
  kChargino700_30,
  kChargino800_10,
  kChargino800_30,
  kChargino900_30,
};

// Names of background, signal and data samples
const map<EBackground, string> backgroundTitle = {
  {kQCD       , "QCD"                       },
  {kZmumuJets , "Z#rightarrow#mu#mu + jets" },
  {kTT        , "tt"                        },
  {kVV        , "VV"                        },
  {kWmunuJets , "W#rightarrow#mu#nu + jets" },
  {kZnunuJets , "Z#rightarrow#nu#nu + jets" },
};

const map<EBackground, string> backgroundName = {
  {kQCD       , "QCD"       },
  {kZmumuJets , "Z_mumu"    },
  {kTT        , "tt"        },
  {kVV        , "VV"        },
  {kWmunuJets , "W_munu"    },
  {kZnunuJets , "Z_nunu"    },
};

const map<ESignal, string> signalTitle = {
  {kWino_M_300_cTau_3, "Wino m=300 c#tau=3"},
  {kWino_M_300_cTau_10, "Wino m=300 c#tau=10"},
  {kWino_M_300_cTau_30, "Wino m=300 c#tau=30"},
  {kWino_M_500_cTau_10, "Wino m=500 c#tau=10"},
  {kWino_M_500_cTau_20, "Wino m=500 c#tau=20"},
  {kWino_M_650_cTau_10, "Wino m=650 c#tau=10"},
  {kWino_M_650_cTau_20, "Wino m=650 c#tau=20"},
  {kWino_M_800_cTau_10, "Wino m=800 c#tau=10"},
  {kWino_M_800_cTau_20, "Wino m=800 c#tau=20"},
  {kWino_M_1000_cTau_10, "Wino m=1000 c#tau=10"},
  {kWino_M_1000_cTau_20, "Wino m=1000 c#tau=20"},
  {kTaggerSignalNoPU      , "Tagger signal no PU"       },
  {kTaggerBackgroundNoPU  , "Tagger background no PU"   },
  {kTaggerSignalWithPU    , "Tagger signal with PU"     },
  {kTaggerBackgroundWithPU, "Tagger background with PU" },
  {kChargino300_1,  "Chargino, 300 GeV, 1 cm"   },
  {kChargino300_10, "Chargino, 300 GeV, 10 cm"  },
  {kChargino400_1,  "Chargino, 400 GeV, 1 cm"   },
  {kChargino500_1,  "Chargino, 500 GeV, 1 cm"   },
  {kChargino500_10, "Chargino, 500 GeV, 10 cm"  },
  {kChargino600_10, "Chargino, 600 GeV, 10 cm"  },
  {kChargino700_10, "Chargino, 700 GeV, 10 cm"  },
  {kChargino700_30, "Chargino, 700 GeV, 30 cm"  },
  {kChargino800_10, "Chargino, 800 GeV, 10 cm"  },
  {kChargino800_30, "Chargino, 800 GeV, 30 cm"  },
  {kChargino900_30, "Chargino, 900 GeV, 30 cm"  },
  
};

const map<ESignal, string> signalName = {
  {kWino_M_300_cTau_3, "Wino_m300_ct3"},
  {kWino_M_300_cTau_10, "Wino_m300_ct10"},
  {kWino_M_300_cTau_30, "Wino_m300_ct30"},
  {kWino_M_500_cTau_10, "Wino_m500_ct10"},
  {kWino_M_500_cTau_20, "Wino_m500_ct20"},
  {kWino_M_650_cTau_10, "Wino_m650_ct10"},
  {kWino_M_650_cTau_20, "Wino_m650_ct20"},
  {kWino_M_800_cTau_10, "Wino_m800_ct10"},
  {kWino_M_800_cTau_20, "Wino_m800_ct20"},
  {kWino_M_1000_cTau_10, "Wino_m1000_ct10"},
  {kWino_M_1000_cTau_20, "Wino_m1000_ct20"},
  {kTaggerSignalNoPU      , "Tagger_signal_no_PU"       },
  {kTaggerBackgroundNoPU  , "Tagger_background_no_PU"   },
  {kTaggerSignalWithPU    , "Tagger_signal_with_PU"     },
  {kTaggerBackgroundWithPU, "Tagger_background_with_PU" },
  {kChargino300_1,  "Chargino_300_1"   },
  {kChargino300_10, "Chargino_300_10"  },
  {kChargino400_1,  "Chargino_400_1"   },
  {kChargino500_1,  "Chargino_500_1"   },
  {kChargino500_10, "Chargino_500_10"  },
  {kChargino600_10, "Chargino_600_10"  },
  {kChargino700_10, "Chargino_700_10"  },
  {kChargino700_30, "Chargino_700_30"  },
  {kChargino800_10, "Chargino_800_10"  },
  {kChargino800_30, "Chargino_800_30"  },
  {kChargino900_30, "Chargino_900_30"  },
};

const map<ESignal, string> signalShortName = {
  {kWino_M_300_cTau_3, "300_3"},
  {kWino_M_300_cTau_10, "300_10"},
  {kWino_M_300_cTau_30, "300_30"},
  {kWino_M_500_cTau_10, "500_10"},
  {kWino_M_500_cTau_20, "500_20"},
  {kWino_M_650_cTau_10, "650_10"},
  {kWino_M_650_cTau_20, "650_20"},
  {kWino_M_800_cTau_10, "800_10"},
  {kWino_M_800_cTau_20, "800_20"},
  {kWino_M_1000_cTau_10, "1000_10"},
  {kWino_M_1000_cTau_20, "1000_20"},
  {kTaggerSignalNoPU      , "Tagger_signal_no_PU"       },
  {kTaggerBackgroundNoPU  , "Tagger_background_no_PU"   },
  {kTaggerSignalWithPU    , "Tagger_signal_with_PU"     },
  {kTaggerBackgroundWithPU, "Tagger_background_with_PU" },
  {kChargino300_1,  "Chargino_300_1"   },
  {kChargino300_10, "Chargino_300_10"  },
  {kChargino400_1,  "Chargino_400_1"   },
  {kChargino500_1,  "Chargino_500_1"   },
  {kChargino500_10, "Chargino_500_10"  },
  {kChargino600_10, "Chargino_600_10"  },
  {kChargino700_10, "Chargino_700_10"  },
  {kChargino700_30, "Chargino_700_30"  },
  {kChargino800_10, "Chargino_800_10"  },
  {kChargino800_30, "Chargino_800_30"  },
  {kChargino900_30, "Chargino_900_30"  },
};

enum EData{
  kSignalRegion,
  kControlRegion
};

constexpr initializer_list<EData> datas = { kSignalRegion, kControlRegion };

const map<EData, string> dataTitle = {
  {kControlRegion, "CR"},
  {kSignalRegion, "SR"},
};

const vector<int> years = { 2017, 2018 };

const string commonDataSuffix = lxplus ? "treeProducerXtracks/" : "";

/*
// On lxplus for tagger studies
const map<int, string> baseDataPath = {
  {2017, "/afs/cern.ch/work/j/jniedzie/private/disapp_tracks/data/"},
  {2018, "/afs/cern.ch/work/j/jniedzie/private/disapp_tracks/data/"},
};
*/

// On lxplus for main analysis
const map<int, string> baseDataPath = {
  {2017, lxplus ? "/eos/cms/store/group/phys_exotica/xtracks/"  : "../data/"},
  {2018, lxplus ? "/eos/cms/store/group/phys_exotica/xtracks/"  : "../data/"},
};

// Path to trees with background, signal and data samples (also determines which samples will be merged)
const map<EBackground, map<int, pair<string, vector<string>>>> inFileNameBackground = {
  { kQCD, {
      { 2017, {lxplus ? "6Mar2019-Hadded/QCD_HT" : "MC-SR/QCD_HT",
        {"100to200/", "200to300/", "300to500/", "500to700/", "700to1000/", "1000to1500/", "1500to2000/", "2000toInf/"}}},
      { 2018, {"7Sep2019/Calibrated-MC-SR-2018-Hadded-2/QCD_HT",
        {"100to200/", "200to300/", "300to500/", "500to700/", "700to1000/", "1000to1500/", "1500to2000/", "2000toInf/"}}},
  }},
  { kZmumuJets, {
    { 2017, {lxplus ? "6Mar2019-Hadded/DYJetsM50_HT" : "MC-SR/DYJetsM50_HT",
      {"100to200/", "100to200e/", "200to400/", "200to400e/", "400to600/", "400to600e/", "600to800/", "800to1200/", "1200to2500/", "2500toInf/"}}},
    { 2018, {"7Sep2019/Calibrated-MC-SR-2018-Hadded-2/DYJetsToLL_M50_HT", // on lxplus
      {"100to200/", "200to400/", "400to600/", "400to600_ext2/", "600to800/", "800to1200/", "1200to2500/", "2500toInf/"}}},
  }},
  { kTT, {
    { 2017, {lxplus ? "6Mar2019-Hadded/" : "MC-SR/",
      {"TTHad/", "TTLep/", "TTSemi/", "T_tch/", "T_tWch/", "TBar_tch/", "TBar_tWch/"}}},
    { 2018, {"7Sep2019/Calibrated-MC-SR-2018-Hadded-2/",
      {"TTHad_pow/", "TTLep_pow/", "TTSemi_pow/", "T_tch/", "T_tWch_noFullyHad/", "TBar_tch/", "TBar_tWch_noFullyHad/", "T_sch_lep/" }}},
  }},
  { kVV, {
    { 2017, {lxplus ? "6Mar2019-Hadded/" : "MC-SR/",
      {"WW/", "WZ/", "ZZ/"}}},
    { 2018, {"7Sep2019/Calibrated-MC-SR-2018-Hadded-2/",
      {"WW/", "WZ/", "ZZ/"}}},
  }},
  { kWmunuJets, {
    { 2017, {lxplus ? "6Mar2019-Hadded/WJets_HT" : "MC-SR/WJets_HT",
      {"100to200/", "200to400/", "400to600/", "600to800/", "800to1200/", "1200to2500/", "2500toInf/"}}},
    { 2018, {"7Sep2019/Calibrated-MC-SR-2018-Hadded-2/WJetsToLNu_HT",
      {"100to200/", "200to400/", "400to600/", "600to800/", "800to1200/", "1200to2500/", "2500toInf/"}}},
  }},
  { kZnunuJets, {
    { 2017, {lxplus ? "6Mar2019-Hadded/ZvvJets_HT" : "MC-SR/ZvvJets_HT",
      {"100to200/", "200to400/", "400to600/", "600to800/", "800to1200/", "1200to2500/", "2500toInf/"}}},
    { 2018, {"7Sep2019/Calibrated-MC-SR-2018-Hadded-2/ZvvJets_HT",
      {"100to200/", "200to400/", "400to600/", "600to800/", "800to1200/", "1200to2500/", "2500toInf/"}}},
  }},
};

const map<ESignal, map<int, pair<string, vector<string>>>> inFileNameSignal = {
  { kWino_M_300_cTau_3, {
    { 2017, {lxplus ? "6Mar2019-Hadded/Wino_M_" : "SIG-SR/Wino_M_", {"300_cTau_3/"}}},
    { 2018, {lxplus ? "7Sep2019/Calibrated-SIG-SR-2018-Hadded/Wino_M_" : "SIG-SR/Wino_M_", {"300_cTau_3/"}}},
  }},
  { kWino_M_300_cTau_10, {
    { 2017, {lxplus ? "SIG-SR-old-Hadded/" : "SIG-SR/Wino_M_", {"300GeV10cm/"}}},
    { 2018, {lxplus ? "7Sep2019/Calibrated-SIG-SR-2018-Hadded/Wino_M_" : "SIG-SR/Wino_M_", {"300_cTau_10/"}}},
  }},
  { kWino_M_300_cTau_30, {
    { 2017, {lxplus ? "6Mar2019-Hadded/Wino_M_" : "SIG-SR/Wino_M_", {"300_cTau_30/"}}},
    { 2018, {lxplus ? "7Sep2019/Calibrated-SIG-SR-2018-Hadded/Wino_M_" : "SIG-SR/Wino_M_", {"300_cTau_30/"}}},
  }},
  { kWino_M_500_cTau_10, {
    { 2017, {lxplus ? "6Mar2019-Hadded/Wino_M_" : "SIG-SR/Wino_M_", {"500_cTau_10/"}}},
    { 2018, {lxplus ? "7Sep2019/Calibrated-SIG-SR-2018-Hadded/Wino_M_" : "SIG-SR/Wino_M_", {"500_cTau_10/"}}},
  }},
  { kWino_M_500_cTau_20, {
    { 2017, {lxplus ? "6Mar2019-Hadded/Wino_M_" : "SIG-SR/Wino_M_", {"500_cTau_20/"}}},
    { 2018, {lxplus ? "7Sep2019/Calibrated-SIG-SR-2018-Hadded/Wino_M_" : "SIG-SR/Wino_M_", {"500_cTau_20/"}}},
  }},
  { kWino_M_650_cTau_10, {
    { 2017, {lxplus ? "6Mar2019-Hadded/Wino_M_" : "SIG-SR/Wino_M_", {"650_cTau_10/"}}},
    { 2018, {lxplus ? "7Sep2019/Calibrated-SIG-SR-2018-Hadded/Wino_M_" : "SIG-SR/Wino_M_", {"650_cTau_10/"}}},
  }},
  { kWino_M_650_cTau_20, {
    { 2017, {lxplus ? "6Mar2019-Hadded/Wino_M_" : "SIG-SR/Wino_M_", {"650_cTau_20/"}}},
    { 2018, {lxplus ? "7Sep2019/Calibrated-SIG-SR-2018-Hadded/Wino_M_" : "SIG-SR/Wino_M_", {"650_cTau_20/"}}},
  }},
  { kWino_M_800_cTau_10, {
    { 2017, {lxplus ? "6Mar2019-Hadded/Wino_M_" : "SIG-SR/Wino_M_", {"800_cTau_10/"}}},
    { 2018, {lxplus ? "7Sep2019/Calibrated-SIG-SR-2018-Hadded/Wino_M_" : "SIG-SR/Wino_M_", {"800_cTau_10/"}}},
  }},
  { kWino_M_800_cTau_20, {
    { 2017, {lxplus ? "6Mar2019-Hadded/Wino_M_" : "SIG-SR/Wino_M_", {"800_cTau_20/"}}},
    { 2018, {lxplus ? "7Sep2019/Calibrated-SIG-SR-2018-Hadded/Wino_M_" : "SIG-SR/Wino_M_", {"800_cTau_20/"}}},
  }},
  { kWino_M_1000_cTau_10, {
    { 2017, {lxplus ? "6Mar2019-Hadded/Wino_M_" : "SIG-SR/Wino_M_", {"1000_cTau_10/"}}},
    { 2018, {lxplus ? "7Sep2019/Calibrated-SIG-SR-2018-Hadded/Wino_M_" : "SIG-SR/Wino_M_", {"1000_cTau_10/"}}},
  }},
  { kWino_M_1000_cTau_20, {
    { 2017, {lxplus ? "6Mar2019-Hadded/Wino_M_" : "SIG-SR/Wino_M_", {"1000_cTau_20/"}}},
    { 2018, {lxplus ? "7Sep2019/Calibrated-SIG-SR-2018-Hadded/Wino_M_" : "SIG-SR/Wino_M_", {"1000_cTau_20/"}}},
  }},
  { kTaggerSignalNoPU, {
    { 2017, {"taggerStudy/", {"signal/noPU/"}}},
    { 2018, {"taggerStudy/", {"signal/noPU/"}}},
  }},
  { kTaggerBackgroundNoPU, {
    { 2017, {"taggerStudy/", {"background/noPU/"}}},
    { 2018, {"taggerStudy/", {"background/noPU/"}}},
  }},
  { kTaggerSignalWithPU, {
    { 2017, {"taggerStudy/", {"signal/withPU/"}}},
    { 2018, {"taggerStudy/", {"signal/withPU/"}}},
  }},
  { kTaggerBackgroundWithPU, {
    { 2017, {"taggerStudy/", {"background/withPU/"}}},
    { 2018, {"taggerStudy/", {"background/withPU/"}}},
  }},
  { kChargino300_1, {
    { 2017, {"7Sep2019/Calibrated-SIG-SR-new-2018-Hadded/", {"Wino_300GeV1cm/"}}},
    { 2018, {"7Sep2019/Calibrated-SIG-SR-new-2018-Hadded/", {"Wino_300GeV1cm/"}}},
  }},
  { kChargino300_10, {
    { 2017, {"7Sep2019/Calibrated-SIG-SR-new-2018-Hadded/", {"Wino_300GeV10cm/"}}},
    { 2018, {"7Sep2019/Calibrated-SIG-SR-new-2018-Hadded/", {"Wino_300GeV10cm/"}}},
  }},
  { kChargino400_1, {
    { 2017, {"7Sep2019/Calibrated-SIG-SR-new-2018-Hadded/", {"Wino_400GeV1cm/"}}},
    { 2018, {"7Sep2019/Calibrated-SIG-SR-new-2018-Hadded/", {"Wino_400GeV1cm/"}}},
  }},
  { kChargino500_10, {
    { 2017, {lxplus ? "7Sep2019/Calibrated-SIG-SR-new-2018-Hadded/" : "SIG-SR/", {"Wino_500GeV10cm/"}}},
    { 2018, {"7Sep2019/Calibrated-SIG-SR-new-2018-Hadded/", {"Wino_500GeV10cm/"}}},
  }},
  { kChargino500_1, {
    { 2017, {"7Sep2019/Calibrated-SIG-SR-new-2018-Hadded/", {"Wino_500GeV1cm/"}}},
    { 2018, {"7Sep2019/Calibrated-SIG-SR-new-2018-Hadded/", {"Wino_500GeV1cm/"}}},
  }},
  { kChargino600_10, {
    { 2017, {"7Sep2019/Calibrated-SIG-SR-new-2018-Hadded/", {"Wino_600GeV10cm/"}}},
    { 2018, {"7Sep2019/Calibrated-SIG-SR-new-2018-Hadded/", {"Wino_600GeV10cm/"}}},
  }},
  { kChargino700_10, {
    { 2017, {"7Sep2019/Calibrated-SIG-SR-new-2018-Hadded/", {"Wino_700GeV10cm/"}}},
    { 2018, {"7Sep2019/Calibrated-SIG-SR-new-2018-Hadded/", {"Wino_700GeV10cm/"}}},
  }},
  { kChargino700_30, {
    { 2017, {"7Sep2019/Calibrated-SIG-SR-new-2018-Hadded/", {"Wino_700GeV30cm/"}}},
    { 2018, {"7Sep2019/Calibrated-SIG-SR-new-2018-Hadded/", {"Wino_700GeV30cm/"}}},
  }},
  { kChargino800_10, {
    { 2017, {"7Sep2019/Calibrated-SIG-SR-new-2018-Hadded/", {"Wino_800GeV10cm/"}}},
    { 2018, {"7Sep2019/Calibrated-SIG-SR-new-2018-Hadded/", {"Wino_800GeV10cm/"}}},
  }},
  { kChargino800_30, {
    { 2017, {"7Sep2019/Calibrated-SIG-SR-new-2018-Hadded/", {"Wino_800GeV30cm/"}}},
    { 2018, {"7Sep2019/Calibrated-SIG-SR-new-2018-Hadded/", {"Wino_800GeV30cm/"}}},
  }},
  { kChargino900_30, {
    { 2017, {"7Sep2019/Calibrated-SIG-SR-new-2018-Hadded/", {"Wino_900GeV30cm/"}}},
    { 2018, {"7Sep2019/Calibrated-SIG-SR-new-2018-Hadded/", {"Wino_900GeV30cm/"}}},
  }},
};

const map<EData, map<int, pair<string, vector<string>>>> inFileNameData = {
  { kSignalRegion, {
    { 2017, {lxplus ? "6Mar2019-Hadded/DATA-SR-2017/MET_Run2017" : "Data-SR/tree_MET_Run2017",
      {"B_31Mar2018/", "C_31Mar2018/", "D_31Mar2018/", "E_31Mar2018/", "F_31Mar2018/"}}},
    { 2018, {lxplus ? "7Sep2019/Calibrated-DATA-SR-2018-Hadded/MET_Run2018" : "Data-SR/MET_Run2018",
      {"A_17Sep2018/", "B_17Sep2018/", "C_17Sep2018/",
        "D_PromptReco_v2_1/", "D_PromptReco_v2_2/", "D_PromptReco_v2_3/", "D_PromptReco_v2_4/"}}},
  }},
  { kControlRegion, {
    { 2017, {"Data-CR/MET_Run2018",
      {"A/"}}},
    { 2018, {"Data-CR/MET_Run2018",
      {"A/"}}},
  }},
};

const map<ESignal, double> signalCrossSectionTwoTracks = { // (fb)
  { kWino_M_300_cTau_3  , 190   },
  { kWino_M_300_cTau_10 , 190   },
  { kWino_M_300_cTau_30 , 190   },
  { kWino_M_500_cTau_10 , 22    },
  { kWino_M_500_cTau_20 , 22    },
  { kWino_M_650_cTau_10 , 6.4   },
  { kWino_M_650_cTau_20 , 6.4   },
  { kWino_M_800_cTau_10 , 2.2   },
  { kWino_M_800_cTau_20 , 2.2   },
  { kWino_M_1000_cTau_10, 0.62  },
  { kWino_M_1000_cTau_20, 0.62  },
  { kTaggerSignalNoPU      , 190},
  { kTaggerBackgroundNoPU  , 190},
  { kTaggerSignalWithPU    , 190},
  { kTaggerBackgroundWithPU, 190},
  { kChargino300_1 , 190   },
  { kChargino300_10, 190   },
  { kChargino400_1 , 59    },
  { kChargino500_1 , 22    },
  { kChargino500_10, 22    },
  { kChargino600_10, 9.5   },
  { kChargino700_10, 4.4   },
  { kChargino700_30, 4.4   },
  { kChargino800_10, 2.2   },
  { kChargino800_30, 1.15  },
  { kChargino900_30, 1.15  },
};

const map<ESignal, double> signalCrossSectionOneTrack = { // (fb)
  { kWino_M_300_cTau_3  , 380   },
  { kWino_M_300_cTau_10 , 380   },
  { kWino_M_300_cTau_30 , 380   },
  { kWino_M_500_cTau_10 , 45    },
  { kWino_M_500_cTau_20 , 45    },
  { kWino_M_650_cTau_10 , 13    },
  { kWino_M_650_cTau_20 , 13    },
  { kWino_M_800_cTau_10 , 4.6   },
  { kWino_M_800_cTau_20 , 4.6   },
  { kWino_M_1000_cTau_10, 1.3   },
  { kWino_M_1000_cTau_20, 1.3   },
  { kTaggerSignalNoPU      , 380},
  { kTaggerBackgroundNoPU  , 380},
  { kTaggerSignalWithPU    , 380},
  { kTaggerBackgroundWithPU, 380},
  { kChargino300_1 , 387  },
  { kChargino300_10, 387  },
  { kChargino400_1 , 121  },
  { kChargino500_1 , 46   },
  { kChargino500_10, 46   },
  { kChargino600_10, 20   },
  { kChargino700_10, 9.5  },
  { kChargino700_30, 9.5  },
  { kChargino800_10, 4.8  },
  { kChargino800_30, 4.8  },
  { kChargino900_30, 2.5  },
};

template <class T>
class range
{
public:
  range(T _min=-inf, T _max=inf) : min(_min), max(_max){
    if(min > max){
      throw logic_error("You try to set min grater than max!");
    }
  }
  
  inline T GetMin() const {return min;}
  inline T GetMax() const {return max;}
  
  inline bool IsInside(T val) const {
    if(val >= min && val <= max) return true;
    else return false;
  }
  
  inline bool IsOutside(T val) const {
    if(val < min || val > max) return true;
    else return false;
  }
  
  inline void Print(){
    cout<<"( "<<min<<" -- "<<max<<" )";
  }
  
private:
  T min;
  T max;
};

// Constants for tracker layers
const int nLayers = 14;
const double pixelBarrelZsize = 265; // mm
const double trackerZsize = 2700; // mm
const double layerR[nLayers] = { 29, 68, 109, 160, 250, 340, 430, 520, 610, 696, 782, 868, 965, 1080 };
const double stripModuleZlength = 300;

const vector<range<double>> layerRanges = { // mm
  // pixel barrel
  range<double>(20, 50), // 0 (29)
  range<double>(50, 90), // 1 (68)
  range<double>(90, 145), // 2 (109)
  range<double>(145, 200), // 3 (160)
  
  // strips (TIB)
  range<double>(200, 300),  // 4 (250)
  range<double>(300, 380),  // 5 (340)
  range<double>(380, 460),  // 6 (430)
  range<double>(460, 555),  // 7 (520)
  
  // strips (TOB)
  range<double>(555, 650),  // 8 (610)
  range<double>(650, 740),  // 9 (696)
  range<double>(740, 825),  // 10 (782)
  range<double>(825, 920),  // 11 (868)
  range<double>(920, 1025), // 12 (965)
  range<double>(1025, 1110),// 13 (1080)
};

const vector<vector<range<double>>> diskRanges = { // mm
  // PXEC
  {range<double>(280, 350)}, // 0
  {range<double>(350, 440)}, // 1
  {range<double>(440, 600)}, // 2
  
  // TID
  {range<double>(600, 760),  // 3
    range<double>(760, 790),  //
    range<double>(790, 820),  //
    range<double>(820, 850)},  //
  
  {range<double>(850, 880),  // 4?
    range<double>(880, 905),  //
    range<double>(905, 920),  //
    range<double>(920, 960)},  //
  
  {range<double>(960, 980),   // 5?
    range<double>(980, 1020),  //
    range<double>(1020, 1050), //
    range<double>(1050, 1080), //
    range<double>(1080, 1200)}, //
  
  // TEC
  {range<double>(1200, 1390)},// 6
  {range<double>(1390, 1530)},// 7
  {range<double>(1530, 1670)},// 8
  {range<double>(1670, 1810)},// 9
  {range<double>(1810, 1950)},// 10
  {range<double>(1950, 2130)},// 11
  {range<double>(2130, 2250),// 12
    range<double>(2250, 2310)},//
  {range<double>(2310, 2450),// 13
    range<double>(2450, 2540)},//
  {range<double>(2540, 2670),// 14
    range<double>(2670, 2730)},//
};

inline size_t GetDisksArrayIndex(int index, int signZ)
{
  return diskRanges.size()+signZ*(index+1);
}

const double solenoidField = 3.8; // T


enum EVar{
  kCustom,
  
  // per event variables
  kNvertices,
  kNisoTracks,
  kNjets,
  kNjets30,
  kNjets30a,
  kMetSumEt,
  kMetPt,
  kMetMass,
  kMetEta,
  kMetPhi,
  kMetJetDphi,
  kNhelices,
  
  // per track variables
  kTrackNclusters,
  kTrackTotalDedx,
  kTrackDedxPerCluster,
  kTrackPt,
  kTrackEta,
  kTrackPhi,
  kTrackCaloEm,
  kTrackCaloHad,
  kTrackDxy,
  kTrackDz,
  kTrackCharge,
  kTrackMass,
  kTrackPid,
  kTrackMissingOuterTrackerHits,
  kTrackPixelHits,
  kTrackTrackerHits,
  kTrackRelativeIsolation,
  kTrackAbsoluteIsolation,
  kTrackMetDphi,
  kTrackDedxPerHit,
  kTrackMinDedx,
  kTrackTrackerLayers,
  
  // per jet variables
  kJetPt,
  kJetEta,
  kJetPhi,
  kJetTrackDr,
  kJetCHF,
  kJetNHF,
  
  // per helix variables
  kHelixX,
  kHelixY,
  kHelixZ,
  kHelixPx,
  kHelixPy,
  kHelixPz,
  kHelixCharge,
  
  
  // per track per layer variables
  kDedx,  ///< dE/dx per layer
  kSizeX, ///< X cluster size in each layer
  kSizeY  ///< Y cluster size in each layer
};

inline int BackColor(EBackground bck){
  return TColor::GetColor(backColors[bck][0],backColors[bck][1],backColors[bck][2]);
}

inline int SignalColor(ESignal sig){
  return TColor::GetColor(signalColors[sig][0],signalColors[sig][1],signalColors[sig][2]);
}

inline int DataColor(EData data){
  return TColor::GetColor(dataColors[data][0],dataColors[data][1],dataColors[data][2]);
}

//,,,{255, 215, 180},{0, 0, 128},{128, 128, 128},{255, 255, 255},{0, 0, 0}

//{2,63,165},{125,135,185},{190,193,212},{214,188,192},{187,119,132},{142,6,59},{74,111,227},{133,149,225},{181,187,227},{230,175,185},{224,123,145},{211,63,106},{17,198,56},{141,213,147},{198,222,199},{234,211,198},{240,185,141},{239,151,8},{15,207,192},{156,222,214},{213,234,231},{243,225,235},{246,196,225},{247,156,212}

//,{157,204,0},{194,0,136},{0,51,128},{255,164,5},{255,168,187},{66,102,0},{255,0,16},{94,241,242},{0,153,143},{224,255,102},{116,10,255},{153,0,0},{255,255,128},{255,255,0},{255,80,5}, {240,163,255}, {0,117,220}, {153,63,0}, {76,0,92}, {25,25,25},{0,92,49}, {43,206,72}, {255,204,153}, {128,128,128}, {148,255,181}, {143,124,0}


// title, nBins, min, max, logY
const map<EVar, tuple<string, int, double, double, bool>> settings =
{
  {kNvertices,    make_tuple("N good vertices",           20,   0,  80,   false)},
  {kNjets,        make_tuple("N jets",                    10,   0,  10,   false)},
  {kNisoTracks,   make_tuple("N iso tracks",              5,    0,  5,    false)},
  {kMetPt,        make_tuple("MET p_{T}",                 25,  200, 1000, true)},
  {kMetJetDphi,   make_tuple("#Delta#phi (jet, MET)",    25, -3.5, 3.5,  false)},
  {kTrackMetDphi, make_tuple("#Delta#phi (track, MET)",  25, -3.5, 3.5,  false)},

  {kMetSumEt,   make_tuple("MET sum Et",100,-20,5000, false)},
  {kMetMass,    make_tuple("MET mass",100,-10e-6,10e6, false)},
  {kMetEta,     make_tuple("MET eta",100,-3.5,3.5, false)},
  {kMetPhi,     make_tuple("MET phi",100,-3.5,3.5, false)},
  {kNhelices,   make_tuple("N fitted helices",10,0,10, false)},
  {kNjets30,    make_tuple("N jets with pt > 30, |eta|<2.4",15,0.0,15, false)},
  {kNjets30a,   make_tuple("N jets with pt > 30, |eta|<4.7",15,0.0,15, false)},

  {kTrackPt,                      make_tuple("Track p_{T} (GeV)",   25,  0,   1000, true)},
  {kTrackTrackerLayers,           make_tuple("N tracker layers",    20,  0,   20,   false)},
  {kTrackRelativeIsolation,       make_tuple("Relative isolation",  50,  0,   0.5,  true)},
  {kTrackDedxPerHit,              make_tuple("dE/dx per hit",       50,  0,   10,   false)},
  {kTrackMinDedx,                 make_tuple("min dE/dx"    ,       50,  0,   10,   false)},
  {kTrackCaloEm,                  make_tuple("EM calo energy",      50,  0,   10,   true)},
  {kTrackCaloHad,                 make_tuple("Hadron calo energy",  50,  0,   10,   true)},

  {kTrackEta,                     make_tuple("Track #eta",                    50, -3.0, 3.0,  false)},
  {kTrackPhi,                     make_tuple("Track #phi",50,-3.5,3.5, false)},
  {kTrackDxy,                     make_tuple("Displacement in XY",100,-0.02,0.02, false)},
  {kTrackDz,                      make_tuple("Displacement in Z",100,-0.02,0.02, false)},
  {kTrackCharge,                  make_tuple("Charge dist",100,-10,10, false)},
  {kTrackMass,                    make_tuple("Mass dist",500,0.0,0.25, false)},
  {kTrackPid,                     make_tuple("PDG PID",441,-220,220, false)},
  {kTrackMissingOuterTrackerHits, make_tuple("Missing outer tracker hits",20,0,20, false)},
  {kTrackPixelHits,               make_tuple("N pixel hits",10,0,10, false)},
  {kTrackTrackerHits,             make_tuple("N tracker hits",40,0,40, false)},
  {kTrackAbsoluteIsolation,       make_tuple("Absolute isolation in dR=0.3",50,0,1, false)},
  {kTrackNclusters,               make_tuple("N detIDs per track",  20,0,22, false)},
  {kTrackTotalDedx,               make_tuple("total dedx per track",50,0,140, false)},
  {kTrackDedxPerCluster,          make_tuple("total dedx per track / n clusters",50,0,14, false)},

  {kHelixX , make_tuple("Helix origin X",600,-300,300, false)},
  {kHelixY , make_tuple("Helix origin Y",600,-300,300, false)},
  {kHelixZ , make_tuple("Helix origin Z",600,-300,300, false)},

  {kHelixPx , make_tuple("Helix momentum X",2000,-1000,1000, false)},
  {kHelixPy , make_tuple("Helix momentum Y",2000,-1000,1000, false)},
  {kHelixPz , make_tuple("Helix momentum Z",2000,-1000,1000, false)},

  {kHelixCharge , make_tuple("Helix charge",3,-1,1, false)},

  {kJetPt       , make_tuple("Jet p_{T} (GeV)",       50,  0.0, 1000, true)},
  {kJetEta      , make_tuple("Jet #eta",              50, -3.0, 3.0,  false)},
  {kJetPhi      , make_tuple("Jet #phi",              50, -3.5, 3.5,  false)},
  {kJetTrackDr  , make_tuple("#DeltaR (jet, track)",  100, 0,   8,    true)},
  {kJetCHF      , make_tuple("Jet f_{CH}",            100, 0,   1.0,  true)},
  {kJetNHF      , make_tuple("Jet f_{NH}",            100, 0,   1.0,  true)},


  {kDedx , make_tuple("dedx",50,0,13, false)},
  {kSizeX , make_tuple("sizeX",10,0,13, false)},
  {kSizeY , make_tuple("sizeY",10,0,13, false)}
};

inline bool IsPerEventVariable(EVar var)
{
  if(var == kNvertices || var == kNisoTracks || var == kNjets || var == kNjets30 ||
     var == kNjets30a || var == kMetSumEt || var == kMetPt || var == kMetMass ||
     var == kMetEta || var == kMetPhi || var == kNhelices)
    return true;
  return false;
}

inline bool IsPerTrackVariable(EVar var)
{
  if(var == kTrackNclusters || var == kTrackTotalDedx || var == kTrackDedxPerCluster || var == kTrackPt ||
     var == kTrackEta || var == kTrackPhi || var == kTrackCaloEm || var == kTrackCaloHad ||
     var == kTrackDxy || var == kTrackDz || var == kTrackCharge || var == kTrackMass || var == kTrackPid ||
     var == kTrackMissingOuterTrackerHits || var == kTrackPixelHits || var == kTrackTrackerHits ||
     var == kTrackTrackerLayers || var == kTrackRelativeIsolation || var == kTrackAbsoluteIsolation ||
     var == kTrackMetDphi || var == kTrackDedxPerHit || var == kDedx || var == kTrackMinDedx)
    return true;
  return false;
}

inline bool IsPerJetVariable(EVar var)
{
  if(var == kJetPt || var == kJetEta || var ==  kJetPhi || var ==  kJetTrackDr || var ==  kJetCHF ||
     var ==  kJetNHF  || var == kMetJetDphi)
    return true;
  return false;
}

inline bool IsPerLayerVariable(EVar var)
{
  if(var == kDedx || var == kSizeX || var == kSizeY)
    return true;
  return false;
};

inline bool IsPerHelixVariable(EVar var)
{
  if(var == kHelixX || var == kHelixY || var == kHelixZ ||
     var == kHelixPx || var == kHelixPy || var == kHelixPz ||
     var == kHelixCharge)
    return true;
  return false;
};

template<typename ContainerType>
inline void EraseFast(ContainerType &container, size_t index)
{
  // ensure that we're not attempting to access out of the bounds of the container.
  assert(index < container.size());
  
  //Swap the element with the back element, except in the case when we're the last element.
  if (index + 1 != container.size())
    std::swap(container[index], container.back());
  
  //Pop the back of the container, deleting our old element.
  container.pop_back();
}

inline double GetRadiusInMagField(double px, double py, double B)
{
  return sqrt(pow(px,2)+pow(py,2))/(B*3)*10;
}

inline double GetPtInMagField(double radius, double B)
{
  return B*3*radius/10;
}

inline int RandInt(int min, int max)
{
  return min + static_cast<int>(rand()) /( static_cast<int>(RAND_MAX/(max-min)));
}

inline double RandDouble(double min, double max)
{
  return min + static_cast<double>(rand()) /( static_cast<double>(RAND_MAX/(max-min)));
}

inline int RandSign()
{
  if(static_cast<double>(rand())/static_cast<double>(RAND_MAX) < 0.5) return -1;
  return 1;
}

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

/// Calculate duration between two events
/// \param t0 Start time
/// \param t1 End time
/// \return Difference between events t0 and t1 in seconds
template<class T>
double duration(T t0,T t1)
{
  auto elapsed_secs = t1-t0;
  typedef std::chrono::duration<float> float_seconds;
  auto secs = std::chrono::duration_cast<float_seconds>(elapsed_secs);
  return secs.count();
}

/// Returns current time
inline std::chrono::time_point<std::chrono::steady_clock> now()
{
  return std::chrono::steady_clock::now();
}

inline double GetRofT(double R0, double a, double tMin, double t, int charge){
  return R0 + charge*sgn(tMin)*a*(t-tMin);
}

inline double GetSofT(double s0, double b, double tMin, double t, int charge){
  return s0 + charge*sgn(tMin)*b*(t-tMin);
}

template <typename T>
string to_string_with_precision(const T a_value, const int n = 6)
{
  ostringstream out;
  out.precision(n);
  out << fixed << a_value;
  return out.str();
}

inline string exec(const char* cmd)
{
  array<char, 128> buffer;
  string result;
  unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
  if (!pipe) throw runtime_error("popen() failed!");
  
  while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr){ result += buffer.data();}
  return result;
}

}

using namespace xtracks;

#endif /* Helpers_h */
