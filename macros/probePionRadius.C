const double trackLengthStep = 1;
const int nTrackSteps = 1000;

const int nPhiSteps = 10000;
const double zMax = 500;

const double magField = 4;
const double trackerRadius = 1080;
const double trackLengthMax = 1080;


const bool logZ = false;

double getPtForRadius(double R){
  return (R*3*magField)/10;
}

double getUniformRand(double min, double max){
  return min + static_cast<double>(rand()) /( static_cast<double>(RAND_MAX/(max-min)));
}

TH1D *pionCharginoDeltaPhi;

double getRandFromDist(){
  return pionCharginoDeltaPhi->GetRandom();
}

double getX(double l, double R, double a, int sign=1){
  return (-pow(l,3) + l*R*R + 2*a*a*l*R*R - sign*sqrt(pow(l,4)*R*R + a*a*pow(l,4)*R*R - 2*l*l*pow(R,4) - 2*a*a*l*l*pow(R,4) + pow(R,6) + a*a*pow(R,6)))/(2*(-l*l + R*R + a*a*R*R));
}

inline int RandSign(){
  if(static_cast<double>(rand())/static_cast<double>(RAND_MAX) < 0.5) return -1;
  return 1;
}

pair<double, double> getMaxPionRadii(double trackLength, double phi){
  double lineCoefficient = tan(phi+TMath::Pi()/2.);
  
  double x1 = getX(trackLength, trackerRadius, lineCoefficient, +1);
  double x2 = getX(trackLength, trackerRadius, lineCoefficient, -1);
  
  double y1 = lineCoefficient*(x1 - trackLength);
  double y2 = lineCoefficient*(x2 - trackLength);
  
  double pionRadius1 = trackerRadius - sqrt(x1*x1 + y1*y1);
  double pionRadius2 = trackerRadius - sqrt(x2*x2 + y2*y2);
  
  return make_pair(pionRadius1, pionRadius2);
}

TH1D* getCharginoNlayersHist(){
  TFile *inFile = TFile::Open("nCharginoLayers.root");
  TH1D *hist = (TH1D*)inFile->Get("nTrackerLayers");
  return hist;
  /*
  TFile *inFile = TFile::Open("/Users/Jeremi/Documents/Fellow/disappearingTracks/pionSignal/tree_friend.root");
  
  TTree *tree = (TTree*)inFile->Get("CharginoAnalyzer/tree");
  if(!tree){
    cout<<"no tree"<<endl;
    exit(0);
  }
  vector<int> *nLayers = nullptr;
  tree->SetBranchAddress("chargino_nTrackerLayers", &nLayers);
  
  TH1D *hist = new TH1D("nTrackerLayers", "nTrackerLayers",100, 0, 100);
  
  for(int i=0; i<10000; i++){
    tree->GetEntry(i);
    for(int n : *nLayers){
      hist->Fill(n);
    }
  }
  
  hist->SaveAs("nCharginoLayers.root");
  
  return hist;*/
}

const int nLayers = 14;
const double layerR[nLayers] = { 29, 68, 109, 160, 250, 340, 430, 520, 610, 696, 782, 868, 965, 1080 };

void probePionRadius()
{
  gStyle->SetOptStat(0);
  
  TFile *inFile = TFile::Open("/Users/Jeremi/Documents/Fellow/disappearingTracks/disappTracks/results/tmp.root");
  pionCharginoDeltaPhi = (TH1D*)inFile->Get("delta_phi_pion_chargino");
  pionCharginoDeltaPhi->Smooth(2);
  
  TH1D *pionPt =  (TH1D*)inFile->Get("pion_pt");
  
  TH2D *hist = new TH2D("maxPt_vs_pionLength", "Uniform #Delta#varphi distribution",
                        600, 0, 1200, // χ length (mm)
                        500, 0, 1500);// pt (MeV)
  
  TH2D *histRealistic = new TH2D("maxPt_vs_pionLength_realistic", "Realistic #Delta#varphi distribution",
                                 600, 0, 1200, // χ length (mm)
                                 500, 0, 1500);// pt (MeV)
  
  const int nLines = 4;
  TGraph *phiLines[nLines];
  for(int i=0; i<nLines; i++){ phiLines[i] = new TGraph(); }
  double linesValues[nLines] = { 0, TMath::Pi()/5., TMath::Pi()/3., TMath::Pi()/2.};
  double linesColors[nLines] = { kBlack, kRed, kMagenta, kCyan };
  int iPoint[nLines] = {1};
  
  TCanvas *canvas = new TCanvas("canvas", "canvas", 1000, 1500);
  canvas->Divide(1,2);
  TCanvas *canvasPassing = new TCanvas("canvasPassing", "canvasPassing", 800, 600);
  
  TH1D *passingHist = new TH1D("passingHist", "passingHist", trackLengthMax/trackLengthStep, 0, trackLengthMax);

  TH1D *nLayersDist = getCharginoNlayersHist();
  
//  for(double trackLength=0; trackLength<trackLengthMax; trackLength+=trackLengthStep){
  
  double totalPassing = 0;
  double totalCount = 0;
  
  for(int trackStep=0; trackStep<nTrackSteps; trackStep++){
    int iLayer = nLayersDist->GetRandom();
    
    double trackLength;
    if(iLayer < nLayers-1){
      trackLength = getUniformRand(layerR[iLayer],layerR[iLayer+1]);
    }
    else{
      trackLength = layerR[nLayers-1];
    }
    
    double passing = 0;
    double counts = 0;
    
    for(int phiStep=0; phiStep<nPhiSteps; phiStep++){
      double phi = getUniformRand(0, TMath::Pi());
      auto pionRadii = getMaxPionRadii(trackLength, phi);
      
      hist->Fill(trackLength, getPtForRadius(pionRadii.first));
      hist->Fill(trackLength, getPtForRadius(pionRadii.second));
      
      for(int iLine=0; iLine<nLines; iLine++){
        if(fabs(phi-linesValues[iLine])<0.001){
          phiLines[iLine]->SetPoint(iPoint[iLine]++, trackLength, getPtForRadius(pionRadii.first));
          phiLines[iLine]->SetPoint(iPoint[iLine]++, trackLength, getPtForRadius(pionRadii.second));
        }
      }
      
      phi = getRandFromDist();
      pionRadii = getMaxPionRadii(trackLength, phi);
      
      double maxPt1 = getPtForRadius(pionRadii.first);
      double maxPt2 = getPtForRadius(pionRadii.second);
      
      histRealistic->Fill(trackLength, maxPt1);
      histRealistic->Fill(trackLength, maxPt2);
      
      double pt = pionPt->GetRandom();
      
      if(pt < maxPt1){
        passing += 1;
        totalPassing += 1;
      }
      if(pt < maxPt2){
        passing += 1;
        totalPassing += 1;
      }
      counts += 2;
      totalCount += 2;
    }
    passing /= counts;
    
    passingHist->Fill(trackLength, passing);
  }
  totalPassing /= totalCount;
  cout<<"Total looper probability: "<<totalPassing<<endl;
  
  canvasPassing->cd();
  passingHist->SetMarkerColor(kGreen+2);
  passingHist->SetMarkerSize(1.0);
  passingHist->SetMarkerStyle(20);
  passingHist->Rebin(20);
  passingHist->Scale(1/20.);
  passingHist->Sumw2(false);
  passingHist->GetXaxis()->SetTitle("#chi^{#pm} track length (mm)");
  passingHist->GetYaxis()->SetTitle("Looper probability (%)");
  passingHist->SetTitle("");
  passingHist->Draw("P");
  
  canvas->cd(1);
  hist->Draw("colz");
  hist->GetXaxis()->SetTitle("#chi^{#pm} track length (mm)");
  hist->GetYaxis()->SetTitle("Max pion p_{T} (MeV)");
  hist->GetZaxis()->SetRangeUser(0, zMax);
  gPad->SetLogz(logZ);
  
  canvas->cd(2);
  histRealistic->Draw("colz");
  histRealistic->GetXaxis()->SetTitle("#chi^{#pm} track length (mm)");
  histRealistic->GetYaxis()->SetTitle("Max pion p_{T} (MeV)");
  histRealistic->GetZaxis()->SetRangeUser(0, zMax);
  gPad->SetLogz(logZ);

  TLegend *leg = new TLegend(0.1, 0.7, 0.25, 0.9);

  
  
  for(int iLine=0; iLine<nLines; iLine++){
    phiLines[iLine]->SetMarkerColor(linesColors[iLine]);
    phiLines[iLine]->SetLineColor(linesColors[iLine]);
    phiLines[iLine]->SetMarkerStyle(20);
    phiLines[iLine]->SetMarkerSize(0.2);
    canvas->cd(1); phiLines[iLine]->Draw("sameP");
    canvas->cd(2); phiLines[iLine]->Draw("sameP");
  }
  
  leg->AddEntry(phiLines[0], "#Delta#varphi = 0", "l");
  leg->AddEntry(phiLines[1], "#Delta#varphi = #pi/4", "l");
  leg->AddEntry(phiLines[2], "#Delta#varphi = #pi/3", "l");
  leg->AddEntry(phiLines[3], "#Delta#varphi = #pi/2", "l");
  canvas->cd(1); leg->Draw("same");
  canvas->cd(2); leg->Draw("same");
  
  double minY = 0;
  
  double a = -getPtForRadius(trackerRadius/2)/(2*trackerRadius);
  double b = getPtForRadius(trackerRadius/2)/2;
  
  for(int iLayer=0; iLayer<nLayers-2; iLayer++){
    TLine *line = new TLine(layerR[iLayer], minY, layerR[iLayer], getPtForRadius((trackerRadius-layerR[iLayer])/2.));
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    
    TLatex *text = new TLatex(layerR[iLayer]-12, a*layerR[iLayer]+b, (""+to_string(iLayer)).c_str());
    text->SetTextSize(0.04);
    
    canvas->cd(1);
    line->Draw("same");
    text->Draw("same");
    canvas->cd(2);
    line->Draw("same");
    text->Draw("same");
    
    canvasPassing->cd();
    
    double y = passingHist->GetBinContent(passingHist->GetXaxis()->FindFixBin(layerR[iLayer]));
    
    line = new TLine(layerR[iLayer], 0, layerR[iLayer], y);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    
    text = new TLatex(layerR[iLayer]-12, y/2., to_string(iLayer).c_str());
    text->SetTextSize(0.04);
    line->Draw("same");
    text->Draw("same");
  }
  
//  canvas->SaveAs("max_pion_pt_vs_layer.pdf");
//  canvasPassing->SaveAs("looper_probability.pdf");
}
