const vector<string> backgroundNames = { "QCD", "Z_mumu", "tt", "VV", "W_munu", "Z_nunu"};
const vector<string> signalNames = { "Wino_m300_ct3", "Wino_m300_ct10", "Wino_m300_ct30", "Wino_m500_ct10", "Wino_m500_ct20", "Wino_m650_ct10", "Wino_m650_ct20", "Wino_m800_ct10", "Wino_m800_ct20", "Wino_m1000_ct10", "Wino_m1000_ct20" };

map<string, vector<int>> backColors = {
  {"QCD",     {230, 25 , 75 }},
  {"Z_mumu",  {60 , 180, 75 }},
  {"tt",      {0  , 130, 200}},
  {"VV",      {245, 130, 48 }},
  {"W_munu",  {145, 30 , 180}},
  {"Z_nunu",  {70 , 240, 240}},
};

map<string, vector<int>> signalColors = {
  {"Wino_m300_ct3"    , {128, 128, 0  }},
  {"Wino_m300_ct10"   , {170, 110, 40 }},
  {"Wino_m300_ct30"   , {128, 0  , 0  }},
  {"Wino_m500_ct10"   , {170, 110, 40 }},
  {"Wino_m500_ct20"   , {128, 0  , 0  }},
  {"Wino_m650_ct10"   , {170, 110, 40 }},
  {"Wino_m650_ct20"   , {128, 0  , 0  }},
  {"Wino_m800_ct10"   , {170, 110, 40 }},
  {"Wino_m800_ct20"   , {128, 0  , 0  }},
  {"Wino_m1000_ct10"  , {170, 110, 40 }},
  {"Wino_m1000_ct20"  , {128, 0  , 0  }},
};

inline int BackColor(string bck){
  return TColor::GetColor(backColors[bck][0],backColors[bck][1],backColors[bck][2]);
}

inline int SignalColor(string sig){
  return TColor::GetColor(signalColors[sig][0],signalColors[sig][1],signalColors[sig][2]);
}

void drawTriggerTurnOn()
{
  gStyle->SetOptStat(0);
  
  TFile *inFile = TFile::Open("../results/triggerTurnOn.root");
  
  map<string, TH1D*> histsBackground;
  map<string, TH1D*> histsSignal;
  
  bool first = true;
  
  TCanvas *c1 = new TCanvas("c1", "c1", 1500, 1000);
  c1->Divide(2,2);
  
  c1->cd(1);
  for(string name : backgroundNames){
    histsBackground[name] = (TH1D*)inFile->Get(("trigger_turnon_"+name).c_str());
    
    if(!histsBackground[name]) continue;
    
    histsBackground[name]->SetLineColor(BackColor(name));
    histsBackground[name]->Draw(first ? "" : "same");
    
    first = false;
  }
  
  c1->cd(3);
  first = true;
  for(string name : signalNames){
    histsSignal[name] = (TH1D*)inFile->Get(("trigger_turnon_"+name).c_str());
    
    if(!histsSignal[name]) continue;
    
    histsSignal[name]->SetLineColor(SignalColor(name));
    histsSignal[name]->Draw(first ? "" : "same");
    
    first = false;
  }
  
  c1->cd(2);
  
  TLegend *legend = new TLegend(0.5, 0.1, 0.9, 0.3);
  
  TH1D *histAllBackgrounds = (TH1D*)inFile->Get("trigger_turnon_backgrounds");
  histAllBackgrounds->SetLineColor(kBlue);
  histAllBackgrounds->GetXaxis()->SetTitle("MET p_{t} (GeV)");
  histAllBackgrounds->GetYaxis()->SetTitle("Trigger efficiency");
  histAllBackgrounds->Draw();
  legend->AddEntry(histAllBackgrounds, "Backgrounds", "le");
  
  TH1D *histAllSignals = (TH1D*)inFile->Get("trigger_turnon_signals");
  histAllSignals->SetLineColor(kRed);
  histAllSignals->Draw("same");
  legend->AddEntry(histAllSignals, "Signals", "le");
  
  legend->Draw("same");
  
  TLine *line = new TLine(200, 0, 200, 1);
  line->Draw("same");
  
  
}
