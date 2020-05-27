//  FitterConfig.cpp
//
//  Created by Jeremi Niedziela on 09/01/2019.

#include "ConfigManager.hpp"

ConfigManager config("init");

ConfigManager::ConfigManager(string path)
{
  for(EBackground iBck : backgrounds) {  runBackground.push_back(false);}
  for(ESignal iSig : signals)         {  runSignal.push_back(false);    }
  for(EData iData : datas)            {  runData.push_back(false);      }
 
  if(path=="init") return;
  
  ifstream infile(path);
  
  string line;
  while (getline(infile, line)){
    if(line.find("#") == 0) continue;
    if(line.find("*") == 0) continue;
    
    if(line.find_first_not_of(' ') == string::npos) continue;
    
    stringstream lineStream(line);
    string key, value;

    getline(lineStream, key, ':');
    getline(lineStream, value);
    
    key.erase(std::remove_if(key.begin(), key.end(), ::isspace), key.end());
    value.erase(std::remove_if(value.begin(), value.end(), ::isspace), value.end());
    
    if(key == "analysis_category")                category          = value;
    else if(key == "secondary_category")          secondaryCategory = value;
    else if(key == "output_path")   	            outputPath        = value;
    
    else if(key == "do_QCD")                      runBackground[kQCD]                 = stoi(value);
    else if(key == "do_Zmm")                      runBackground[kZmumuJets]           = stoi(value);
    else if(key == "do_tops")                     runBackground[kTT]                  = stoi(value);
    else if(key == "do_dibosons")                 runBackground[kVV]                  = stoi(value);
    else if(key == "do_Wmv")                      runBackground[kWmunuJets]           = stoi(value);
    else if(key == "do_Zvv")                      runBackground[kZnunuJets]           = stoi(value);
    
    else if(key == "do_300_3")                    runSignal[kWino_M_300_cTau_3]       = stoi(value);
    else if(key == "do_300_10")                   runSignal[kWino_M_300_cTau_10]      = stoi(value);
    else if(key == "do_300_30")                   runSignal[kWino_M_300_cTau_30]      = stoi(value);
    else if(key == "do_500_10")                   runSignal[kWino_M_500_cTau_10]      = stoi(value);
    else if(key == "do_500_20")                   runSignal[kWino_M_500_cTau_20]      = stoi(value);
    else if(key == "do_650_10")                   runSignal[kWino_M_650_cTau_10]      = stoi(value);
    else if(key == "do_650_20")                   runSignal[kWino_M_650_cTau_20]      = stoi(value);
    else if(key == "do_800_10")                   runSignal[kWino_M_800_cTau_10]      = stoi(value);
    else if(key == "do_800_20")                   runSignal[kWino_M_800_cTau_20]      = stoi(value);
    else if(key == "do_1000_10")                  runSignal[kWino_M_1000_cTau_10]     = stoi(value);
    else if(key == "do_1000_20")                  runSignal[kWino_M_1000_cTau_20]     = stoi(value);
    else if(key == "do_tagger_signal_noPU")       runSignal[kTaggerSignalNoPU]        = stoi(value);
    else if(key == "do_tagger_signal_noPU_pion_removed") runSignal[kTaggerSignalNoPUpionRemoved] = stoi(value);
    else if(key == "do_tagger_background_noPU")   runSignal[kTaggerBackgroundNoPU]    = stoi(value);
    else if(key == "do_tagger_signal_withPU")     runSignal[kTaggerSignalWithPU]      = stoi(value);
    else if(key == "do_tagger_background_withPU") runSignal[kTaggerBackgroundWithPU]  = stoi(value);
    else if(key == "do_chargino_300_1")           runSignal[kChargino300_1]           = stoi(value);
    else if(key == "do_chargino_300_10")          runSignal[kChargino300_10]          = stoi(value);
    else if(key == "do_chargino_300_30")          runSignal[kChargino300_30]          = stoi(value);
    else if(key == "do_chargino_400_1")           runSignal[kChargino400_1]           = stoi(value);
    else if(key == "do_chargino_500_1")           runSignal[kChargino500_1]           = stoi(value);
    else if(key == "do_chargino_500_10")          runSignal[kChargino500_10]          = stoi(value);
    else if(key == "do_chargino_500_30")          runSignal[kChargino500_30]          = stoi(value);
    else if(key == "do_chargino_500_10_noMETfilter") runSignal[kChargino500_10_noMETfilter]          = stoi(value);
    else if(key == "do_chargino_500_10_newGT") runSignal[kChargino500_10_newGT]          = stoi(value);
    else if(key == "do_chargino_600_10")          runSignal[kChargino600_10]          = stoi(value);
    else if(key == "do_chargino_700_10")          runSignal[kChargino700_10]          = stoi(value);
    else if(key == "do_chargino_700_30")          runSignal[kChargino700_30]          = stoi(value);
    else if(key == "do_chargino_800_10")          runSignal[kChargino800_10]          = stoi(value);
    else if(key == "do_chargino_800_30")          runSignal[kChargino800_30]          = stoi(value);
    else if(key == "do_chargino_900_1")           runSignal[kChargino900_1]           = stoi(value);
    else if(key == "do_chargino_900_10")          runSignal[kChargino900_10]          = stoi(value);
    else if(key == "do_chargino_900_30")          runSignal[kChargino900_30]          = stoi(value);
    
    else if(key == "do_SR")                       runData[kSignalRegion]              = stoi(value);
    else if(key == "do_CR")                       runData[kControlRegion]             = stoi(value);
    
    else params[key] =  stod(value);
  }
}
