vector<TFile*> files = {
  TFile::Open("../results/clustersSignalNoPU.root"),
  TFile::Open("../results/nCharginoLayers.root"),
  TFile::Open("../results/clustersSignalWithPU.root"),
};


map<string, tuple<int, string, int, int>> histSettings = {
//  hist name                       pad x axis title                            in file canvas
  { "delta_phi_pion_chargino"     , {1, "#Delta#varphi (#chi^{#pm}, #pi^{#pm})" , 0,     0} },
  { "nTrackerLayers"              , {2, "#chi N tracker layers"                 , 1,     0} },
  { "pion_pt"                     , {3, "Pion p_{t} (MeV)"                      , 0,     0} },

  { "pion_simhits_z"              , {1, "Hit z-coordinate dist. (mm)"           , 0,     1} },
  { "pion_simhits_range_z"        , {2, "Last hit z position (mm)"              , 0,     1} },
  { "pion_initial_radius"         , {3, "Initial pion helix radius (mm)"        , 0,     1} },
  { "pion_final_radius"           , {4, "Final pion helix radius (mm)"          , 0,     1} },
  
  { "noise_n_clusters"            , {1, "Number of noise hits"                  , 0,     2} },
  { "tracker_cluster_charge_TIB"  , {3, "Charge deposit in TIB (a.u.)"          , 0,     2} },
  { "pion_cluster_charge_TIB"     , {3, "Charge deposit in TIB (a.u.)"          , 0,     2} },
  { "tracker_cluster_charge_TOB"  , {4, "Charge deposit in TOB (a.u.)"          , 0,     2} },
  { "pion_cluster_charge_TOB"     , {4, "Charge deposit in TOB (a.u.)"          , 0,     2} },
  { "tracker_cluster_charge_TID"  , {5, "Charge deposit in TID (a.u.)"          , 0,     2} },
  { "pion_cluster_charge_TID"     , {5, "Charge deposit in TID (a.u.)"          , 0,     2} },
  { "tracker_cluster_charge_TEC"  , {6, "Charge deposit in TEC (a.u.)"          , 0,     2} },
  { "pion_cluster_charge_TEC"     , {6, "Charge deposit in TEC (a.u.)"          , 0,     2} },
  
  { "middle_seed_hit_delta_phi"   , {1, "Middle seed hit #Delta#varphi"         , 0,     3} },
  { "middle_seed_hit_delta_z"     , {2, "Middle seed hit #Delta Z (mm)"         , 0,     3} },
  { "last_seed_hit_delta_phi"     , {3, "Last seed hit #Delta#varphi"           , 0,     3} },
  { "last_seed_hit_delta_z"       , {4, "Last seed hit #Delta Z (mm)"           , 0,     3} },
  { "next_point_delta_phi"        , {5, "Next hit #Delta#varphi"                , 0,     3} },
  { "next_point_delta_z"          , {6, "Next hit #Delta Z (mm)"                , 0,     3} },
};

void prepareHist(TH1D *hist, string titleX){
  hist->SetLineColor(kViolet+2);
  hist->SetFillColorAlpha(kViolet+2, 0.3);
  hist->SetFillStyle(1000);
  gStyle->SetOptStat(0);
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.12);
  hist->GetXaxis()->SetTitle(titleX.c_str());
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitle("entries");
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->SetTitle("");
}

void drawPlots()
{
  map<string, TH1D*> hists;
  
  vector<TCanvas*> canvas = {
    new TCanvas("canvas0", "canvas0", 1000, 800),
    new TCanvas("canvas1", "canvas1", 1000, 800),
    new TCanvas("canvas2", "canvas2", 1000, 1200),
    new TCanvas("canvas3", "canvas3", 1000, 1200),
  };
  canvas[0]->Divide(2,2);
  canvas[1]->Divide(2,2);
  canvas[2]->Divide(2,3);
  canvas[3]->Divide(2,3);
  
  vector<vector<int>> nDrawn;
  vector<vector<TLegend*>> leg;
  
  for(int i=0; i<canvas.size(); i++){
    vector<int> vec;
    vector<TLegend*> vecLeg;
    for(int j=0; j<7; j++){
      vec.push_back(0);
      vecLeg.push_back(new TLegend(0.5, 0.7, 0.9, 0.9));
    }
    nDrawn.push_back(vec);
    leg.push_back(vecLeg);
  }
  
  for(auto &[name, settings] : histSettings){
    auto &[iPad, titleX, fileIndex, canvasIndex] = settings;
    
    TH1D *hist = (TH1D*)files[fileIndex]->Get(name.c_str());
    prepareHist(hist, titleX);
    
    canvas[canvasIndex]->cd(iPad++);
    if(name=="nTrackerLayers") hist->GetXaxis()->SetRangeUser(0, 20);
//    if(name=="pion_initial_radius") hist->Smooth();

    string legTitle = "";
    if(name.find("pion") != string::npos) legTitle = "Pion";
    else if(canvasIndex == 3) legTitle = "Noise (PU)";
    else legTitle = "Noise";
    
    if(nDrawn[canvasIndex][iPad-1] == 0){
      nDrawn[canvasIndex][iPad-1]++;
      
      hist->DrawNormalized();
      leg[canvasIndex][iPad-1]->AddEntry(hist, legTitle.c_str(), "f");
    }
    else if(nDrawn[canvasIndex][iPad-1]==1){
      nDrawn[canvasIndex][iPad-1]++;
      
      hist->SetFillColorAlpha(kGreen, 0.3);
      hist->SetLineColor(kGreen);
      hist->DrawNormalized("same");
      leg[canvasIndex][iPad-1]->AddEntry(hist, legTitle.c_str(), "f");
      leg[canvasIndex][iPad-1]->Draw("same");
    }
    
    if(name.find("tracker_cluster_charge") != string::npos ||
       name.find("noise_n_clusters") != string::npos){
      TH1D *hist2 = (TH1D*)files[2]->Get(name.c_str());
      prepareHist(hist2, titleX);
      
      hist2->SetFillColorAlpha(kRed, 0.3);
      hist2->SetLineColor(kRed);
      hist2->DrawNormalized("same");
      leg[canvasIndex][iPad-1]->AddEntry(hist2, (legTitle+" (PU)").c_str(), "f");
      leg[canvasIndex][iPad-1]->Draw("same");
    }
  }

  for(int canvasIter=0; canvasIter<canvas.size(); canvasIter++)
    canvas[canvasIter]->SaveAs(("../plots/genLevelPlots_"+to_string(canvasIter)+".pdf").c_str());
  
//  TH1D *pionRadius = new TH1D("pionRadius","pionRadius",100,0,1000);
//  TH1D *pionPt = (TH1D*)files[0]->Get("pion_pt");
//
//  for(int i=1; i<pionPt->GetNbinsX()+1; i++){
//    double pt = pionPt->GetXaxis()->GetBinCenter(i);
//    double val = pionPt->GetBinContent(i);
//    double radius = 10./(3*4) * pt;
//    pionRadius->Fill(radius, val);
//  }
//
//
//  canvas->cd(5);
//
//  pionRadius->Rebin(2);
//  pionRadius->Smooth();
//  pionRadius->Sumw2(false);
//  pionRadius->DrawNormalized("same");
}
