
map<string, tuple<int, string, int, int>> histSettings = {
//  hist name                   pad x axis title                            in file canvas
  { "delta_phi_pion_chargino" , {1, "#Delta#varphi (#chi^{#pm}, #pi^{#pm})" , 0,     0} },
  { "pion_pt"                 , {2, "Pion p_{t} (MeV)"                      , 0,     0} },
  { "nTrackerLayers"          , {3, "#chi N tracker layers"                 , 1,     0} },
  { "pion_range_z"            , {4, "Last hit z position (mm)"              , 0,     0} },
  { "pion_initial_radius"     , {5, "Initial pion helix radius (mm)"        , 0,     0} },
  { "pion_final_radius"       , {6, "Final pion helix radius (mm)"          , 0,     0} },
  { "noise_n_clusters"        , {1, "Number of noise hits"                  , 0,     1} },
  
  { "tracker_cluster_charge_TIB"  , {3, "Charge deposit in TIB (a.u.)"      , 0,     1} },
  { "pion_cluster_charge_TIB"     , {3, "Charge deposit in TIB (a.u.)"      , 0,     1} },
  
  { "tracker_cluster_charge_TOB"  , {4, "Charge deposit in TOB (a.u.)"      , 0,     1} },
  { "pion_cluster_charge_TOB"     , {4, "Charge deposit in TOB (a.u.)"      , 0,     1} },
  
  { "tracker_cluster_charge_TID"  , {5, "Charge deposit in TID (a.u.)"      , 0,     1} },
  { "pion_cluster_charge_TID"     , {5, "Charge deposit in TID (a.u.)"      , 0,     1} },
  
  { "tracker_cluster_charge_TEC"  , {6, "Charge deposit in TEC (a.u.)"      , 0,     1} },
  { "pion_cluster_charge_TEC"     , {6, "Charge deposit in TEC (a.u.)"      , 0,     1} },
  
  { "middle_seed_hit_delta_phi"   , {1, "Middle seed hit #Delta#varphi"     , 0,     2} },
  { "middle_seed_hit_delta_z"     , {2, "Middle seed hit #Delta Z (mm)"     , 0,     2} },
  { "last_seed_hit_delta_phi"     , {3, "Last seed hit #Delta#varphi"       , 0,     2} },
  { "last_seed_hit_delta_z"       , {4, "Last seed hit #Delta Z (mm)"       , 0,     2} },
  { "next_point_delta_phi"     , {5, "Next hit #Delta#varphi"       , 0,     2} },
  { "next_point_delta_z"       , {6, "Next hit #Delta Z (mm)"       , 0,     2} },
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
  vector<TFile*> files = {
    TFile::Open("../results/tmp.root"),
    TFile::Open("nCharginoLayers.root"),
  };
  
  map<string, TH1D*> hists;
  
  vector<TCanvas*> canvas = {
    new TCanvas("canvas0", "canvas0", 1000, 1500),
    new TCanvas("canvas1", "canvas1", 1000, 1500),
    new TCanvas("canvas2", "canvas2", 1000, 1500),
  };
  for(auto c : canvas) c->Divide(2,3);
  
  vector<vector<bool>> first;
  
  vector<vector<TLegend*>> leg;
  
  for(int i=0; i<canvas.size(); i++){
    vector<bool> vec;
    vector<TLegend*> vecLeg;
    for(int j=0; j<7; j++){
      vec.push_back(true);
      vecLeg.push_back(new TLegend(0.5, 0.7, 0.9, 0.9));
    }
    first.push_back(vec);
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
    else legTitle = "Noise";
    
    if(first[canvasIndex][iPad-1]){
      hist->DrawNormalized();
      first[canvasIndex][iPad-1] = false;
      leg[canvasIndex][iPad-1]->AddEntry(hist, legTitle.c_str(), "f");
    }
    else{
      hist->SetFillColorAlpha(kGreen, 0.3);
      hist->DrawNormalized("same");
      leg[canvasIndex][iPad-1]->AddEntry(hist, legTitle.c_str(), "f");
      leg[canvasIndex][iPad-1]->Draw("same");
    }
    
  }

  for(int canvasIter=0; canvasIter<canvas.size(); canvasIter++)
    canvas[canvasIter]->SaveAs(("genLevelPlots_"+to_string(canvasIter)+".pdf").c_str());
  
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
