// ct in cm (*0.03333 to convert to ns)
const double ctMin = 1;
const double ctMax = 30;
const double ctStep = 0.1;

// mass in GeV
const double massMin = 300;
const double massMax = 900;
const double massStep = 1;

bool doExtrapolation = false;

vector<tuple<int, int, double>> getRvalues(string fileName){
  vector<tuple<int, int, double>> rValues;
  int mass, ct;
  double r;
  
  ifstream inFile(fileName);
  while(inFile >> mass >> ct >> r){
    if(r==0) continue;
    rValues.push_back({mass, ct, (ct==1 ? 1000 : 1)*r});
  }
  inFile.close();
  
  return rValues;
}

TGraph2D* getRgraph(const vector<tuple<int, int, double>> &rValues)
{
  TGraph2D *graph = new TGraph2D();
  
  int iPoint=0;
  
  for(auto &[mass, ct, r] : rValues) graph->SetPoint(iPoint++, mass,  ct  , r);
  
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(1.0);
  graph->SetMarkerColor(kRed);
  
  return graph;
}

void getLimitsExtrapolate(string inputPath, string outPath)
{
  auto rValues = getRvalues(inputPath);
  ofstream outFile(outPath);
  
  
  map<int, TGraph*> graphs;
  map<int, TF1*> funs;
  map<int, int> iPoints;
  
  map<int, int> colors = {
    {1  , kRed    },
    {10 , kGreen  },
    {30 , kBlue   },
  };
  
  TLegend *legend = new TLegend(0.5, 0.7, 0.9, 0.9);
  
  for(auto &[mass, ct, r] : rValues){
    if(graphs.count(ct) == 0){
      graphs[ct] = new TGraph();
      graphs[ct]->SetMarkerSize(1.0);
      graphs[ct]->SetMarkerStyle(20);
      graphs[ct]->SetMarkerColor(colors[ct]);
      legend->AddEntry(graphs[ct], ("c#tau="+to_string(ct)+" cm").c_str(), "p");
    }
    if(!funs.count(ct)){
      funs[ct] = new TF1(("fun_"+to_string(ct)).c_str(), "[0]*exp([1]*x)", 0, 1000);
      funs[ct]->SetParameter(0, 0.03);
      funs[ct]->SetParLimits(0, 0, 1);
      funs[ct]->SetParameter(1, 0.01);
//      funs[ct]->SetParameter(2, 0);
      funs[ct]->SetLineColor(colors[ct]);
    }
    cout<<iPoints[ct]<<"\t"<<mass<<"\t"<<r<<endl;
    graphs[ct]->SetPoint(iPoints[ct]++, mass, r);
    
  }
  
  TCanvas *canvas = new TCanvas("canvas", "canvas", 1000, 1000);
  canvas->Divide(2, 2);
  
  map<int, TGraph*> limits;
  
  canvas->cd(1);
  
  for(int ct : {1, 10, 30}){
    if(ct==1){
      graphs[ct]->Draw("AP");
      gPad->SetLogy();
      graphs[ct]->GetYaxis()->SetRangeUser(1e-2, 1e7);
      graphs[ct]->GetYaxis()->SetTitle("r");
      graphs[ct]->GetXaxis()->SetTitle("mass (GeV)");
    }
    else{
      graphs[ct]->Draw("Psame");
    }
    
    graphs[ct]->Fit(funs[ct]);
    funs[ct]->Draw("same");
    
    double massLimit = funs[ct]->GetX(1.0);
    
    if(!limits.count(ct)){
      limits[ct] = new TGraph();
      limits[ct]->SetMarkerStyle(20);
      limits[ct]->SetMarkerSize(1.0);
      limits[ct]->SetMarkerColor(colors[ct]);
    }
    
    if(isnormal(massLimit)){
      limits[ct]->SetPoint(0, massLimit, ct);
      
      double tau = ct * 0.0333333;
      outFile<<massLimit<<"\t"<<tau<<endl;
    }
    
  }
  legend->Draw();
  
  outFile.close();
  
  canvas->cd(2);
  for(auto &[ct, graph] : limits) graph->Draw(ct==1 ? "AP" : "Psame");
  limits[1]->GetXaxis()->SetLimits(0, 1000);
  limits[1]->GetYaxis()->SetRangeUser(0, 35);
  limits[1]->GetXaxis()->SetTitle("mass (GeV)");
  limits[1]->GetYaxis()->SetTitle("ct (cm)");
  
  
  canvas->cd(3);
  
  TF2 *fun = new TF2("fun","[0]*exp([1]*x) + [2]*pow(y, 2) + [3]*y + [4]",0, 1000, 0, 30);
  fun->SetParameter(0, 34.9384);
  fun->SetParLimits(0, 0, 1);
  fun->SetParameter(1, 0.00859244);
  fun->SetParameter(2, 31.5623);
  //  fun->SetParLimits(2, 0, 0);
  fun->SetParameter(3, -982.408);
  fun->SetParameter(4, 1066.55);
  
  TGraph2D *graph = getRgraph(rValues);
  //  graph->Draw("AP");
  
  graph->Fit(fun);
  
  fun->Draw("surf1");
  graph->Draw("sameP");
  
}



