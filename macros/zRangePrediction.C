
double lengthZ = 2500;    // mm
double zStep = 1;  // mm

int nHelices = 2500;

vector<double> layerYpos = { 30, 70, 110, 160, 250, 340, 430, 520, 610, 700, 780, 868, 965, 1080 };
double yTolerance = 1;

enum ECase {
  kStoppingShrinking,
  kStoppingNotShrinking,
  kNotStoppingShrinking,
  kNotStoppingNotShrinking,
  kNcases
};

map<ECase, vector<double>> params = {
  //                           y         R0        a           s0           b         stop
  {kStoppingShrinking,        {110, 420, 300, 500, -100 , -10 , -50, -500  , 0 , 20 , 0.0005} },
  {kStoppingNotShrinking,     {110, 420, 300, 500, 0    , 0   , -50,   0   , 0 , 0  , 0.0005} },
  {kNotStoppingShrinking,     {110, 420, 300, 500, -100 , -10 , -50, -500  , 0 , 20 , 0.   } },
  {kNotStoppingNotShrinking,  {110, 420, 300, 500, 0    , 0   , -50,   0   , 0 , 0  , 0.   } },
};

map<ECase, tuple<string, int, int>> histParams = {
  //                          title                           color     style
  {kStoppingShrinking,        {"stoppipng shrinking",         kRed      , 1} },
  {kStoppingNotShrinking,     {"stoppipng not shrinking",     kRed      , 2} },
  {kNotStoppingShrinking,     {"not stoppipng shrinking",     kGreen+2  , 1} },
  {kNotStoppingNotShrinking,  {"not stoppipng not shrinking", kGreen+2  , 2} },
};

double getYofZ(double z, double R0, double a, double s0, double b, int sign = -1)
{
  if(fabs(a)<0.0001 || fabs(b)<0.0001) return R0*cos(z/s0);
  
  double t = (-s0 + sign*sqrt(s0*s0+4*z*b))/(2*b);
  return (R0 + a*t)*cos(t);
}

double RandDouble(double min, double max)
{
  return min + static_cast<double>(rand()) /( static_cast<double>(RAND_MAX/(max-min)));
}

void zRangePrediction()
{
  TH1D *hitZpos[kNcases];
  TH1D *lastHitZpos[kNcases];
  
  cout<<"Filling histograms"<<endl;
  
  for(int iCase=0; iCase<kNcases; iCase++){
    string title = "HitZpos"+to_string(iCase);
    hitZpos[iCase]     = new TH1D(title.c_str(), title.c_str(), 50, 0, 3500);
    lastHitZpos[iCase] = new TH1D(("last"+title).c_str(), ("last"+title).c_str(), 50, 0, 3500);
    
    for(int iHelix=0; iHelix<nHelices; iHelix++){
      int iPoint=0;
      
      ECase iiCase = (ECase)iCase;
      
      double yy0  = RandDouble(params[iiCase][0], params[iiCase][1]);
      double R0   = RandDouble(params[iiCase][2], params[iiCase][3]);
      double a    = RandDouble(params[iiCase][4], params[iiCase][5]);
      double s0   = RandDouble(params[iiCase][6], params[iiCase][7]);
      double b    = RandDouble(params[iiCase][8], params[iiCase][9]);
      
      double lastZ = -1;
      
      for(double z=0; z<lengthZ; z+=zStep){
        if(RandDouble(0, 1) < params[iiCase][10]) break;
        double y = yy0 + getYofZ(z, R0, a, s0, b);
        
        
        for(double layerY : layerYpos){
          if(fabs(layerY-y) < yTolerance){
            hitZpos[iCase]->Fill(z);
            lastZ = z;
          }
        }
      }
      
      lastHitZpos[iCase]->Fill(lastZ);
    }
  }
  
  cout<<"Plotting"<<endl;
  
  gStyle->SetOptStat(0);
  
  TCanvas *c1 = new TCanvas("c1","c1", 1000, 1500);
  c1->Divide(1,2);
  
  TLegend *leg = new TLegend(0.5, 0.7, 0.9, 0.9);
  
  hitZpos[0]->GetXaxis()->SetTitle("Hit z (mm)");
  hitZpos[0]->GetYaxis()->SetTitle("entries");
  hitZpos[0]->SetTitle("Hit z-coordinate distribution");
  
  lastHitZpos[0]->GetXaxis()->SetTitle("Last hit z (mm)");
  lastHitZpos[0]->GetYaxis()->SetTitle("entries");
  lastHitZpos[0]->SetTitle("Last hit z-coordinate distribution");
  
  for(int iCase=0; iCase<kNcases; iCase++){
    auto &[title, color, style] = histParams[(ECase)iCase];
    
    c1->cd(1);
    hitZpos[iCase]->Print();
    hitZpos[iCase]->SetLineColor(color);
    hitZpos[iCase]->SetLineStyle(style);
    hitZpos[iCase]->Draw(iCase==0 ? "" : "same");
    hitZpos[iCase]->GetYaxis()->SetRangeUser(0, 2.2*nHelices);
    leg->AddEntry(hitZpos[iCase], title.c_str(), "l");
    
    c1->cd(2);
    lastHitZpos[iCase]->SetLineColor(color);
    lastHitZpos[iCase]->SetLineStyle(style);
    lastHitZpos[iCase]->Draw(iCase==0 ? "" : "same");
    lastHitZpos[iCase]->GetYaxis()->SetRangeUser(0, 0.8*nHelices);
  }
  
  c1->cd(1);
  leg->Draw("same");
  
  c1->SaveAs("../plots/z_range_prediction.pdf");
}
