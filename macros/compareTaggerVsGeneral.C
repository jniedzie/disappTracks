
void compareTaggerVsGeneral()
{
  TFile *taggerFile = TFile::Open("../results/tagger_hist.root");
  TFile *generalFile = TFile::Open("../results/general_tracks_hist.root");
  
  TH1D *taggerHist = (TH1D*)taggerFile->Get("fractionTrueHitsHist");
  TH1D *generalHist = (TH1D*)generalFile->Get("histFractionPionHitsOnLooper");
  
  gStyle->SetOptStat(0);
  
  TLegend *legend = new TLegend(0.6, 0.75, 0.9, 0.9);
  
  taggerHist->Rebin(4);
  taggerHist->SetTitle("");
  taggerHist->GetXaxis()->SetTitle("Fraction of track hits that are true pion hits");
//  taggerHist->Sumw2(false);
  taggerHist->SetLineColor(kGreen+2);
  taggerHist->SetFillColorAlpha(kGreen+2, 0.1);
  taggerHist->SetMarkerStyle(20);
  taggerHist->SetMarkerSize(1.0);
  taggerHist->SetMarkerColor(kGreen+2);
  taggerHist->DrawNormalized("");
  
  legend->AddEntry(taggerHist, "Looper tagger", "F");
  
  gPad->SetLogy(true);
  

  generalHist->Rebin(4);
//  generalHist->Sumw2(false);
  generalHist->SetLineColor(kBlue);
  generalHist->SetFillColorAlpha(kBlue, 0.1);
  generalHist->SetMarkerStyle(20);
  generalHist->SetMarkerSize(1.0);
  generalHist->SetMarkerColor(kBlue);
  generalHist->DrawNormalized("same");
  
  legend->AddEntry(generalHist, "Official loopers", "F");
  
  legend->Draw();
  
}
