
void drawCharginoRecoPlots()
{
  auto inFile = TFile::Open("../results/tmp.root");
  auto nLayers = (TH1D*)inFile->Get("layers_gen_vs_rec");
  
  TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
  canvas->cd();
  
  gStyle->SetOptStat(0);
  nLayers->GetXaxis()->SetTitle("N_{rec}^{layers}");
  nLayers->GetYaxis()->SetTitle("N_{gen}^{layers}");
  nLayers->SetTitle("");
  
  gPad->SetBottomMargin(0.14);
  gPad->SetLeftMargin(0.12);
  
  nLayers->GetXaxis()->SetTitleSize(0.05);
  nLayers->GetYaxis()->SetTitleSize(0.05);
  
  nLayers->Draw("colzText");
 
  auto line = new TLine(0, 0, 10, 10);
  line->SetLineColor(kGreen);
  line->SetLineWidth(3);
  line->SetLineStyle(2);
  
  line->Draw("same");
  
  canvas->SaveAs("charginoReco.pdf");
}
