const double Rpion = 200; // mmm
const double Rpixel = 160; // mm
const double trackerRadius = 1080;
const double magField = 3.8;

//double trackL = 300; // mm

double getPtForRadius(double R){
  return (R*3*magField)/10;
}

double getRadiusFortPt(double pt){
  return 10*pt/(3*magField);
}

double getPixelFraction(double trackL, double Rpion)
{
  if(trackL > Rpixel + Rpion) return 0;
  if(Rpion > Rpixel+trackL) return 0;
  double y = (Rpixel*Rpixel-Rpion*Rpion+trackL*trackL)/(2*trackL);
  double x = sqrt(Rpixel*Rpixel-y*y);
  double h = trackL-y;
  return atan2(x,h)/TMath::Pi();
}

void calcPixelFraction()
{
  TGraph *graph[4];
  for(int i=0; i<4; i++){
    graph[i] = new TGraph();
  }
  TGraph *graphFromDist = new TGraph();
  
  TFile *inFile = TFile::Open("/Users/Jeremi/Documents/Fellow/disappearingTracks/disappTracks/results/clustersSignalNoPU.root");
  TH1D *pionPtHist = (TH1D*)inFile->Get("pion_pt");
  pionPtHist->Smooth();
  
  double pionPt[4] = {200, 300, 400, 500};
  int colors[4] = { kRed, kOrange, kGreen, kBlue };
  const char* names[4] = {"200 MeV", "300 MeV", "400 MeV", "500 MeV"};
  
  int iPoint=0;
  for(double trackL=0; trackL<trackerRadius; trackL+=10){
    
    TH1D *hist = new TH1D(("hist"+to_string(trackL)).c_str(),
                          ("hist"+to_string(trackL)).c_str(),
                          1000, 0, 1);
    
    for(int i=0; i<100000; i++){
      hist->Fill(getPixelFraction(trackL, getRadiusFortPt(pionPtHist->GetRandom())));
    }
    graphFromDist->SetPoint(iPoint, trackL, hist->GetMean());
                 
    for(int i=0; i<4; i++){
      graph[i]->SetPoint(iPoint, trackL, getPixelFraction(trackL, getRadiusFortPt(pionPt[i])));
    }
    iPoint++;
  }
  
//  TLegend *legend = new TLegend(0.5, 0.7, 0.9, 0.9);
//
//  for(int i=0; i<4; i++){
//    graph[i]->SetMarkerStyle(20);
//    graph[i]->SetMarkerColor(colors[i]);
//    graph[i]->SetMarkerSize(0.7);
//
//
//    graph[i]->Draw(i==0 ? "AP" : "Psame");
//    graph[i]->GetXaxis()->SetTitle("Track length (mm)");
//    graph[i]->GetYaxis()->SetTitle("Fraction of track within pixel barrel");
//    graph[i]->GetYaxis()->SetRangeUser(0, 1);
//
//    legend->AddEntry(graph[i], names[i], "p");
//  }
  graphFromDist->SetMarkerSize(0.7);
  graphFromDist->SetMarkerColor(kRed);
  graphFromDist->SetMarkerStyle(20);
  graphFromDist->Draw("AP");
  graphFromDist->GetXaxis()->SetTitle("Track length (mm)");
  graphFromDist->GetYaxis()->SetTitle("Fraction of track within pixel barrel");
  graphFromDist->GetYaxis()->SetRangeUser(0, 0.3);
  
  
  const int nLayers = 14;
  const double layerR[nLayers] = { 29, 68, 109, 160, 250, 340, 430, 520, 610, 696, 782, 868, 965, 1080 };
  
  for(int iLayer=0; iLayer<nLayers; iLayer++){
    TLine *line = new TLine(layerR[iLayer], 0, layerR[iLayer], 0.3);
    line->SetLineColor(kBlack);
    line->SetLineStyle(2);
    
    TLatex *text = new TLatex(layerR[iLayer]-12, 0.2, (""+to_string(iLayer)).c_str());
    text->SetTextSize(0.04);
    
    line->Draw("same");
    text->Draw("same");
    
    
    //     canvasPassing->cd();
    //
    //     double y = passingHist->GetBinContent(passingHist->GetXaxis()->FindFixBin(layerR[iLayer]));
    //
    //     line = new TLine(layerR[iLayer], 0, layerR[iLayer], y);
    //     line->SetLineColor(kBlack);
    //     line->SetLineStyle(2);
    //
    //     text = new TLatex(layerR[iLayer]-12, y/2., to_string(iLayer).c_str());
    //     text->SetTextSize(0.04);
    //     line->Draw("same");
    //     text->Draw("same");
  }
//  legend->Draw();
  
}

