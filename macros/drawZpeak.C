const double Zmass = 91.1876; // GeV
const double Zwidth = 2.4952; // GeV

void drawZline(int shift=0, int style=1)
{
  TLine *zPeakLine = new TLine(Zmass+shift*Zwidth, 0, Zmass+shift*Zwidth, 550);
  zPeakLine->SetLineColor(kBlue);
  zPeakLine->SetLineWidth(1);
  zPeakLine->SetLineStyle(style);
  zPeakLine->Draw();
}

void drawZpeak()
{
  gStyle->SetOptStat(0.0000);
  
  TFile *inFile = TFile::Open("../results/tracksHeatMap.root");
  TH1D *invMass = (TH1D*)inFile->Get("invMass #mu#mu");
  
  TCanvas *c1 = new TCanvas("c1","c1",800, 600);
  c1->cd();
  
  invMass->GetXaxis()->SetRangeUser(50, 130);
  invMass->GetYaxis()->SetRangeUser(0 , 550);
  invMass->SetLineColor(kOrange+2);
  invMass->SetFillColorAlpha(kOrange, 0.3);
  invMass->SetTitle("");
  invMass->Draw();
  
  drawZline(0);
  drawZline(-3, 2);
  drawZline(3, 2);
  
  c1->Update();
  c1->SaveAs("../plots/zPeak.pdf");
}

