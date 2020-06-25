// ct in cm (*0.03333 to convert to ns)
const double ctMin = 1;
const double ctMax = 30;
const double ctStep = 0.1;

// mass in GeV
const double massMin = 300;
const double massMax = 900;
const double massStep = 1;

bool doExtrapolation = true;

map<int, map<int, double>> getTupleFromFile(string fileName){
  map<int, map<int, double>> rValues;
  int mass, ct;
  double r;
    
  ifstream inFile(fileName);
  while(inFile >> mass >> ct >> r){
    if(r==0) continue;
    rValues[mass][ct] = r;
    
    if(ct==1)                 rValues[mass][ct] *= 1000;
    if(mass == 900 && ct==3)  rValues[mass][ct] *= 1000;
    
    
  }
  inFile.close();
  
  return rValues;
}


map<int, double> crossSection = {
// mass x-sec (fb)
  {300  , 387  },
  {400  , 121  },
  {500  , 46   },
  {600  , 20   },
  {650  , 13   },
  {700  , 9.5  },
  {800  , 4.8  },
  {900  , 2.5  },
  {1000 , 1.3  },
};

map<int, int> colors = {
  {1  , kRed    },
  {3  , kBlack  },
  {10 , kGreen  },
  {30 , kBlue   },
};

void getLimitsFromR(string inputPath, string outPath)
{
  auto rValues3x3 = getTupleFromFile(inputPath);
  ofstream outFile(outPath);
  
  /*
  TCanvas *canvas = new TCanvas("canvas", "canvas", 2000, 1000);
  canvas->Divide(2, 2);
  
  TGraph *rateCt1   = new TGraph();
  TGraph *rateCt30  = new TGraph();

  int iPointCt1=0;
  int iPointCt30=0;
  
  cout<<"Known points: "<<endl;
  for(auto &[mass, ctAndR] : rValues3x3){
    for(auto &[ct, r] : ctAndR){
      cout<<"m: "<<mass<<"\tct: "<<ct<<"\tr: "<<r<<endl;
      if(ct == 1)   rateCt1->SetPoint(iPointCt1++, mass, r);
      if(ct == 30)  rateCt30->SetPoint(iPointCt30++, mass, r);
    }
  }
    
  canvas->cd(1);

  rateCt1->SetMarkerStyle(20);
  rateCt1->SetMarkerSize(1.0);
  rateCt1->SetMarkerColor(kRed);
  rateCt1->Draw("ap");
  TF1 *funCt1 = new TF1("funCt1","[0]*exp([1]*x)+[2]",0, 1000);
  funCt1->SetParameter(0, 1);
  funCt1->SetParameter(1, -0.01);
  funCt1->SetParameter(2, 0);
  funCt1->SetLineColor(kRed);
  rateCt1->Fit(funCt1);
  funCt1->Draw("same");

  rateCt30->SetMarkerStyle(20);
  rateCt30->SetMarkerSize(1.0);
  rateCt30->SetMarkerColor(kGreen);
  rateCt30->Draw("psame");
  TF1 *funCt30  = new TF1("funCt30","[0]*exp([1]*x)+[2]",0, 1000);
  funCt30->SetParameter(0, 1);
  funCt30->SetParameter(1, -0.01);
  funCt30->SetLineColor(kGreen);
  rateCt30->Fit(funCt30);
  funCt30->Draw("same");
  
  cout<<"\n\nRate for ct=1, m=300: "<<funCt1->Eval(300)<<endl;
  cout<<"Rate for ct=1, m=500: "<<funCt1->Eval(500)<<endl;
  cout<<"Rate for ct=1, m=900: "<<funCt1->Eval(900)<<endl;
  cout<<"Rate for ct=30, m=300: "<<funCt30->Eval(300)<<endl;
  cout<<"Rate for ct=30, m=900: "<<funCt30->Eval(900)<<endl;
  
  canvas->cd(2);
  
  TGraph *rateM500 = new TGraph();

  rateM500->SetPoint(0, 1,  funCt1->Eval(500));
  rateM500->SetPoint(1, 10, rValues3x3[500][10]);
  rateM500->SetPoint(2, 30, funCt30->Eval(500));
  
  rateM500->SetMarkerStyle(20);
  rateM500->SetMarkerSize(1.0);
  rateM500->SetMarkerColor(kGreen);
  rateM500->Draw("ap");
  
  TF1 *funMass500 = new TF1("funMass500", "[0]+[1]*x+[2]*pow(x, 2)", 0, 30);
  funMass500->SetParameter(0, 0);
  funMass500->SetParameter(1, 1);
  funMass500->SetParameter(2, 1);
  
  rateM500->Fit(funMass500);
  funMass500->Draw("same");

  cout<<"\n\nRate for ct=10, m=500: "<<funMass500->Eval(10)<<endl;
  
  canvas->cd(3);
  
  TGraph *rateM300 = new TGraph();

  rateM300->SetPoint(0, 1, rValues3x3[300][1]);
//  rateM300->SetPoint(1, 10, rValues3x3[500][10]);
  rateM300->SetPoint(1, 30, rValues3x3[300][30]);
  
  rateM300->SetMarkerStyle(20);
  rateM300->SetMarkerSize(1.0);
  rateM300->SetMarkerColor(kGreen);
  rateM300->Draw("ap");
  
  TF1 *funMass300 = new TF1("funMass300", "[0]+[1]*x+[2]*pow(x, 2)", 0, 30);
  funMass300->SetParameter(0, 0);
  funMass300->SetParameter(1, 1);
  funMass300->FixParameter(2, 10*funMass500->GetParameter(2));
  
  rateM300->Fit(funMass300);
  funMass300->Draw("same");
  
  cout<<"\n\nRate for ct=1, m=300: "<<funMass300->Eval(1)<<endl;
  cout<<"Rate for ct=30, m=300: "<<funMass300->Eval(30)<<endl;
  
  TGraph *rateM900 = new TGraph();
  
  rateM900->SetPoint(0, 1, rValues3x3[900][1]);
  //  rateM900->SetPoint(1, 10, rValues3x3[500][10]);
  rateM900->SetPoint(1, 30, rValues3x3[900][30]);
  
  rateM900->SetMarkerStyle(20);
  rateM900->SetMarkerSize(1.0);
  rateM900->SetMarkerColor(kGreen);
  rateM900->Draw("ap");
  
  TF1 *funMass900 = new TF1("funMass900", "[0]+[1]*x+[2]*pow(x, 2)", 0, 30);
  funMass900->SetParameter(0, 0);
  funMass900->SetParameter(1, 1);
  funMass900->FixParameter(2, 10*funMass500->GetParameter(2));
  
  rateM900->Fit(funMass900);
  funMass900->Draw("same");
  
  cout<<"\n\nRate for ct=1, m=900: "<<funMass900->Eval(1)<<endl;
  cout<<"Rate for ct=30, m=900: "<<funMass900->Eval(30)<<endl;
  
  TGraph2D *rGraph = new TGraph2D();
  
  int iPoint2D=0;
  
  TGraph *limits = new TGraph();
  limits->SetMarkerStyle(20);
  limits->SetMarkerSize(1.0);
  limits->SetMarkerColor(kGreen);
  int iPointLimits=0;
  return;
  for(double ct=ctMin; ct<=ctMax+ctStep; ct+=ctStep){
    TGraph *graphForCt = new TGraph();

    graphForCt->SetPoint(0, 300, funMass300->Eval(ct));
    graphForCt->SetPoint(1, 500, funMass500->Eval(ct));
    graphForCt->SetPoint(2, 900, funMass900->Eval(ct));

    TF1 *funCt = new TF1("funCt","[0]*exp([1]*x)+[2]",0, 1000);
    funCt->SetParameter(0, 1);
    funCt->SetParameter(1, -0.01);
    funCt->SetParameter(2, 0);
    graphForCt->Fit(funCt);

    
    
    int closestValue = 999999;
    int bestMass = -1;
    
    for(int mass=massMin; mass<=massMax; mass+=massStep){
      rGraph->SetPoint(iPoint2D++, mass, ct, graphForCt->Eval(mass));
      
      if(fabs(graphForCt->Eval(mass)-1.0) < fabs(closestValue-1)){
        closestValue = graphForCt->Eval(mass);
        bestMass=mass;
      }
    }
    
//    double mLimit = log(1./funCt->GetParameter(0))/funCt->GetParameter(1);
//    cout<<"\n\nmLimit: "<<mLimit<<"\n\n"<<endl;
//    limits->SetPoint(iPointLimits++, mLimit, ct);
    limits->SetPoint(iPointLimits++, bestMass, ct);
    
    if(fabs(ct-10)<ctStep/2.){
      
      cout<<"\n\nFun 300, ct="<<ct<<": "<<funMass300->Eval(ct)<<endl;
      
      canvas->cd(3);
      graphForCt->SetMarkerStyle(20);
      graphForCt->SetMarkerSize(1.0);
      graphForCt->Draw("ap");
      funCt->Draw("same");
    }

  }
  
  canvas->cd(4);
  rGraph->Draw("colz");
  limits->Draw("psame");
  
  return;
  */
  if(doExtrapolation){
  
    TCanvas *canvas = new TCanvas("canvas", "canvas", 1000, 500);
    canvas->Divide(2, 1);
    
    TLegend *legend = new TLegend(0.1, 0.7, 0.5, 0.9);
    
    map<int, TGraph*> graphs;
    map<int, TGraph*> rPoints;
    
    canvas->cd(1);
    
    TGraph *xSecVsMass = new TGraph();
    int iPoint=0;
    for(auto &[mass, xsec] : crossSection) xSecVsMass->SetPoint(iPoint++, mass, xsec);
    
    xSecVsMass->SetMarkerSize(1.0);
    xSecVsMass->SetMarkerStyle(20);
    xSecVsMass->SetMarkerColor(kRed);
    
    TF1 *fun = new TF1("fun", "[0]*exp([1]*x)", 0, 1000);
    fun->SetParameter(0, 14000);
    fun->SetParameter(1, -0.012);
    xSecVsMass->Fit(fun);
    
    xSecVsMass->Draw("AP");
    
    canvas->cd(2);
    for(int ct : {1, 3, 10, 30}){
      graphs[ct] = new TGraph();
      iPoint=0;
      for(auto &[mass, xsec] : crossSection){
        graphs[ct]->SetPoint(iPoint++, mass, rValues3x3[300][ct] * crossSection[300]/xsec);
      }
      
      double xSecLimit = rValues3x3[300][ct] * crossSection[300];
      double massLimit = fun->GetX(xSecLimit);
      
      if(!isnormal(massLimit)) massLimit = 0;
      
      double tau = ct * 0.0333333;
      
      cout<<"ct: "<<ct<<"\tr=1 for x-sec: "<<xSecLimit<<endl;
      cout<<"corresponding mass: "<<massLimit<<endl;
      
      outFile<<massLimit<<"\t"<<tau<<endl;
      
      rPoints[ct] = new TGraph();
      iPoint=0;
      
      for(auto &[mass, ctAndR] : rValues3x3){
        for(auto &[ctt, r] : ctAndR){
          if(ctt != ct) continue;
          rPoints[ct]->SetPoint(iPoint++, mass, r);
        }
      }
      
      graphs[ct]->SetLineColor(colors[ct]);
      graphs[ct]->SetLineWidth(2);
      graphs[ct]->Draw(ct==1 ? "AL" : "Lsame");
      legend->AddEntry(graphs[ct], ("extrapolation c#tau="+to_string(ct)+" cm").c_str(), "L");
      gPad->SetLogy();
      
      
      rPoints[ct]->SetMarkerStyle(20);
      rPoints[ct]->SetMarkerSize(0.5);
      rPoints[ct]->SetMarkerColor(colors[ct]);
      rPoints[ct]->Draw("Psame");
      legend->AddEntry(rPoints[ct], ("r value c#tau="+to_string(ct)+" cm").c_str(), "P");
    }
    graphs[1]->GetYaxis()->SetRangeUser(1e-2, 1e6);
    graphs[1]->GetXaxis()->SetTitle("m (GeV)");
    graphs[1]->GetYaxis()->SetTitle("r");
    legend->Draw();
  }
  else{
    
    TGraph2D *rGraph = new TGraph2D();
    
    int iPoint=0;
    for(auto &[mass, ctAndR] : rValues3x3){
      for(auto &[ct, r] : ctAndR){
        if(r==0) r = 1E-10;
        rGraph->SetPoint(iPoint++, mass, ct, r);
      }
    }
    
    TGraph *limitGraph    = new TGraph();
    TGraph *badLimitGraph = new TGraph();
    iPoint=0;
    int iPointBad=0;
    
    for(double ct=ctMin; ct<=ctMax; ct+=ctStep){

      double bestMassForCt = 0;
      double closestValue = 999999;
      double tau = ct * 0.0333333;
      
      for(double mass=massMin; mass<=massMax; mass+=massStep){

        double value = rGraph->Interpolate(mass, ct);

        if(fabs(value-1) < fabs(closestValue-1)){
          closestValue = value;
          bestMassForCt = mass;
        }
      }
      
      if(fabs(closestValue-1) > 0.1){
        cout<<"WARNING -- value of r closest to 1.0 for ct: "<<ct<<" (cm) is far from it: "<<closestValue<<endl;
        badLimitGraph->SetPoint(iPointBad++, bestMassForCt, ct);
      }
      else{
        limitGraph->SetPoint(iPoint++, bestMassForCt, ct);
      }
      outFile<<bestMassForCt<<"\t"<<tau<<endl;
    }

    rGraph->Draw("colz");
    rGraph->GetXaxis()->SetRangeUser(massMin, massMax);
    rGraph->GetYaxis()->SetRangeUser(ctMin, ctMax);

    limitGraph->SetMarkerColor(kGreen);
    limitGraph->SetMarkerStyle(20);
    limitGraph->SetMarkerSize(0.5);
    limitGraph->Draw("Psame");
    
    badLimitGraph->SetMarkerColor(kRed);
    badLimitGraph->SetMarkerStyle(20);
    badLimitGraph->SetMarkerSize(0.5);
    badLimitGraph->Draw("Psame");
  }
  
  outFile.close();
}
