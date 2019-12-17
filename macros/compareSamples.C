
vector<tuple<string, int, int, int>> samples = {
//  {"300, 1"   , 300 , 1   , kRed      },
//  {"400, 1"   , 400 , 1   , kGreen    },
//  {"500, 1"   , 500 , 1   , kBlue     },
  {"500 GeV, 10 cm"  , 500 , 10  , kBlack     },
  {"700 GeV, 10 cm"  , 700 , 10  , kMagenta+2 },
  {"800 GeV, 10 cm"  , 800 , 10  , kGreen+2   },
  {"700 GeV, 30 cm"  , 700 , 30  , kOrange+2  }
};


map<string, tuple<int, string, int, double, int>> histSettings = {
//  hist name                       pad x axis title                                  canvas maxY rebin
  { "delta_phi_pion_chargino"     , {1, "#Delta#varphi^{gen} (#chi^{#pm}, #pi^{#pm})" , 0, 0.15, 1 } },
  { "delta_theta_pion_chargino"   , {2, "#Delta#theta^{gen} (#chi^{#pm}, #pi^{#pm})"  , 0, 0.15, 1 } },
  { "pion_pt"                     , {3, "Pion p_{t}^{gen} (MeV)"                      , 0, 0.14, 2 } },
  { "pion_pz"                     , {4, "Pion |p_{z}^{gen}| (MeV)"                    , 0, 0.14, 2 } },
  { "chargino_abs_eta"            , {5, "#chi^{#pm} |#eta^{gen}|"                     , 0, 0.15, 1 } },
  { "chargino_n_layers_gen"       , {6, "#chi^{#pm} N_{layers}^{gen}"                 , 0, 0.40, 1 } },

  { "pion_simhits_z"              , {1, "Pion sim hit z-coordinate dist. (mm)"        , 1, 0.06, 2 } },
  { "pion_simhits_range_z"        , {2, "Last pion sim hit z position (mm)"           , 1, 0.07, 2 } },
  { "pion_initial_radius"         , {3, "Initial pion helix radius from sim hits (mm)", 1, 0.16, 1 } },
  { "pion_final_radius"           , {4, "Final pion helix radius from sim hits (mm)"  , 1, 0.16, 1 } },

  { "noise_n_clusters"            , {1, "Number of noise hits"                        , 2, 0.10, 1 } },
  { "tracker_cluster_charge_TIB"  , {3, "Charge deposit in TIB (ADC counts)"          , 2, 0.06, 1 } },
  { "pion_cluster_charge_TIB"     , {3, "Charge deposit in TIB (ADC counts)"          , 2, 0.06, 2 } },
  { "tracker_cluster_charge_TOB"  , {4, "Charge deposit in TOB (ADC counts)"          , 2, 0.04, 1 } },
  { "pion_cluster_charge_TOB"     , {4, "Charge deposit in TOB (ADC counts)"          , 2, 0.04, 2 } },
  { "tracker_cluster_charge_TID"  , {5, "Charge deposit in TID (ADC counts)"          , 2, 0.09, 1 } },
  { "pion_cluster_charge_TID"     , {5, "Charge deposit in TID (ADC counts)"          , 2, 0.09, 2 } },
  { "tracker_cluster_charge_TEC"  , {6, "Charge deposit in TEC (ADC counts)"          , 2, 0.04, 1 } },
  { "pion_cluster_charge_TEC"     , {6, "Charge deposit in TEC (ADC counts)"          , 2, 0.04, 2 } },

  { "middle_seed_hit_delta_phi"   , {1, "Middle seed hit #Delta#varphi"               , 3, 0.19, 2 } },
  { "middle_seed_hit_delta_z"     , {2, "Middle seed hit #Delta Z (mm)"               , 3, 0.16, 3 } },
  { "last_seed_hit_delta_phi"     , {3, "Last seed hit #Delta#varphi"                 , 3, 0.09, 2 } },
  { "last_seed_hit_delta_z"       , {4, "Last seed hit #Delta Z (mm)"                 , 3, 0.40, 2 } },
  { "next_point_delta_phi"        , {5, "Next hit #Delta#varphi"                      , 3, 0.08, 2 } },
  { "next_point_delta_z"          , {6, "Next hit #Delta Z (mm)"                      , 3, 0.18, 3 } },
};

void prepareHist(TH1D *hist, string titleX, int color, double maxY, int rebin){
  hist->SetLineColor(color);
//  hist->SetFillColorAlpha(kViolet+2, 0.3);
//  hist->SetFillStyle(1000);
  gStyle->SetOptStat(0);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.12);
  hist->GetXaxis()->SetTitle(titleX.c_str());
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitle("entries");
  hist->GetYaxis()->SetTitleSize(0.05);
  
  hist->SetTitle("");
  
  hist->Rebin(rebin);
  hist->Scale(1.0/(rebin*hist->GetEntries()));
  
  hist->GetYaxis()->SetRangeUser(0, maxY);
  hist->Sumw2(false);
}

void compareSamples()
{
  map<string, TH1D*> hists;
  
  vector<TCanvas*> canvas = {
    new TCanvas("canvas0", "canvas0", 1000, 1200),
    new TCanvas("canvas1", "canvas1", 1000, 800),
    new TCanvas("canvas2", "canvas2", 1000, 1200),
    new TCanvas("canvas3", "canvas3", 1000, 1200),
  };
  canvas[0]->Divide(2,3);
  canvas[1]->Divide(2,2);
  canvas[2]->Divide(2,3);
  canvas[3]->Divide(2,3);
  
  vector<vector<int>> nDrawn;
  TLegend* leg = new TLegend(0.6, 0.5, 0.9, 0.9);
  
  for(int i=0; i<canvas.size(); i++){
    vector<int> vec;
    for(int j=0; j<7; j++){
      vec.push_back(0);
    }
    nDrawn.push_back(vec);
  }
  
  for(auto &[sampleName, mass, ct, color] : samples){
    TFile *inFile = TFile::Open(("../results/clusters_"+to_string(mass)+"_"+to_string(ct)+".root").c_str());
  
    bool first = true;
    for(auto &[name, settings] : histSettings){
      auto &[iPad, titleX, canvasIndex, maxY, rebin] = settings;
      
      TH1D *hist = (TH1D*)inFile->Get(name.c_str());
      if(!hist){
        cout<<"Hist: "<<name<<"\t not found!!"<<endl;
        continue;
      }
      prepareHist(hist, titleX, color, maxY, rebin);
      
      
      if(name.find("tracker_cluster_charge") != string::npos ||
         name.find("noise_n_clusters") != string::npos){
        
        hist->SetLineStyle(2);
      }
      else if(first){
        leg->AddEntry(hist, sampleName.c_str(), "l");
        first = false;
      }
      
      canvas[canvasIndex]->cd(iPad);
      
      if(nDrawn[canvasIndex][iPad-1] == 0){
        nDrawn[canvasIndex][iPad-1]++;
        hist->Draw();
        
      }
      else{
        nDrawn[canvasIndex][iPad-1]++;
        
        //        hist->SetFillColorAlpha(kGreen, 0.3);
        //        hist->SetLineColor(kGreen);
        hist->Draw("same");
      }
      
//      string legTitle = "";
//      if(name.find("pion") != string::npos) legTitle = "Pion";
//      else if(canvasIndex == 3) legTitle = "Noise (PU)";
//      else legTitle = "Noise";
//
//      if(nDrawn[canvasIndex][iPad-1] == 0){
//        nDrawn[canvasIndex][iPad-1]++;
//
//        hist->Draw();
//        leg[canvasIndex][iPad-1]->AddEntry(hist, legTitle.c_str(), "f");
//      }
//      else if(nDrawn[canvasIndex][iPad-1]==1){
//        nDrawn[canvasIndex][iPad-1]++;
//
//        hist->SetFillColorAlpha(kGreen, 0.3);
//        hist->SetLineColor(kGreen);
//        hist->Draw("same");
//        leg[canvasIndex][iPad-1]->AddEntry(hist, legTitle.c_str(), "f");
//        leg[canvasIndex][iPad-1]->Draw("same");
//      }
//
      
    }
  }
    
  for(int canvasIter=0; canvasIter<canvas.size(); canvasIter++){
    canvas[canvasIter]->cd(1); leg->Draw("same");
    canvas[canvasIter]->SaveAs(("../plots/plotsBySample_"+to_string(canvasIter)+".pdf").c_str());
  }
}
