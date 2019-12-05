
void drawTrackerHeatMap()
{
  TFile *inFile = TFile::Open("../results/tracksHeatMap.root");
  TH2D *heatMapAll = (TH2D*)inFile->Get("all");
  
  TH2D *heatMap3layers = (TH2D*)inFile->Get("3layers");
  TH2D *heatMap4layers = (TH2D*)inFile->Get("4layers");
  TH2D *heatMap5layers = (TH2D*)inFile->Get("5layers");
  TH2D *heatMap6layers = (TH2D*)inFile->Get("6layers");
  
  heatMap3layers->Add(heatMap4layers);
  heatMap3layers->Add(heatMap5layers);
  heatMap3layers->Add(heatMap6layers);
  
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->cd();
  
  gStyle->SetOptStat(0);
  
  heatMapAll->SetTitle("Tracker heat map");
  heatMapAll->GetXaxis()->SetTitleSize(0.05);
  heatMapAll->GetXaxis()->SetTitleOffset(0.9);
  heatMapAll->GetYaxis()->SetTitleSize(0.05);
  heatMapAll->GetYaxis()->SetTitleOffset(0.9);
  heatMapAll->DrawNormalized("colz");
  
  c1->SaveAs("../plots/trackerHeatMap.pdf");
  
  
//  heatMap3layers->Rebin2D(8,8);
//  heatMap3layers->Draw("colz");
  
}
