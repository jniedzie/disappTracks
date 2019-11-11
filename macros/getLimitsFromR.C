// ct in cm (*0.03333 to convert to ns)
const double ctMin = 0;
const double ctMax = 35;
const double ctStep = 0.1;

// mass in GeV
const double massMin = 0;
const double massMax = 1100;
const double massStep = 1;

bool doExtrapolation = true;

map<int, map<int, double>> getTupleFromFile(string fileName){
  map<int, map<int, double>> rValues;
  int mass, ct;
  double r;
    
  ifstream inFile(fileName);
  while (inFile >> mass >> ct >> r){ rValues[mass][ct] = r; }
  inFile.close();
  
  return rValues;
}

// (fb), two track (for one track, multiply by 2
map<int, double> crossSection = {{300, 190}, {500, 22}, {650 ,6.4}, {800, 2.2}, {1000, 0.62}};

void getLimitsFromR(string inputPath, string outPath)
{
  auto rValues3x3 = getTupleFromFile(inputPath);
  ofstream outFile(outPath);
  
  if(doExtrapolation){
    for(int ct : {3, 30}){
      for(int mass : {500, 650, 800, 1000}){
        rValues3x3[mass][ct] = rValues3x3[300][ct] * crossSection[300]/crossSection[mass];
      }
    }
    for(int ct : {3, 10, 20, 30}){
      cout<<"ct: "<<ct<<endl;
      
      TGraph *rVsMassCt = new TGraph();
      
      int iPoint=0;
      for(int mass : {300, 500, 650, 800, 1000}){
        rVsMassCt->SetPoint(iPoint++, mass, rValues3x3[mass][ct]);
      }
      
      TF1 *fun3 = new TF1("fun3","[0]*exp([1]*x)",0, 1000);
      fun3->SetParameter(0, 1);
      fun3->SetParameter(1, -0.01);
      rVsMassCt->Fit(fun3);
      
      double minDiff = 9999;
      int bestMass = 0;
      for(int mass = 300; mass < 1000; mass++){
        double diff = fabs(1-fun3->Eval(mass));
        if(diff < minDiff){
          minDiff = diff;
          bestMass = mass;
        }
      }
      
      double tau = ct * 0.0333333;
      outFile<<bestMass<<"\t"<<tau<<endl;
      cout<<"Best mass: "<<bestMass<<endl;
    }
  }
  else{
    TH2D *rMap = new TH2D("rMap", "rMap",
                          (massMax-massMin)/massStep, massMin, massMax,
                          (ctMax-ctMin)/ctStep+1, ctMin, ctMax);
    
    TGraph *gridGraph = new TGraph();
    gridGraph->SetMarkerStyle(4);
    gridGraph->SetMarkerSize(1.0);
    gridGraph->SetMarkerColor(kRed);
    
    int iPoint=0;
    for(auto &[mass, ctAndR] : rValues3x3){
      for(auto &[ct, r] : ctAndR){
        if(r==0) r = 1E-10;
        rMap->Fill(mass, ct, r);
        gridGraph->SetPoint(iPoint++, mass, ct);
      }
    }
    
    TGraph2D *rGraph3x3 = new TGraph2D(rMap);
    TH2D *rMapInterpolated = new TH2D(*rMap);
    rMapInterpolated->Reset();
    
    for(double ct=ctMin; ct<ctMax; ct+=ctStep){
      
      double bestMassForCt = 0;
      double closestValue = 999999;
      
      for(double mass=massMin; mass<massMax; mass+=massStep){
        double value = rGraph3x3->Interpolate(mass, ct);
        rMapInterpolated->Fill(mass, ct, value);
        
        if(fabs(value-1) < fabs(closestValue-1)){
          closestValue = value;
          bestMassForCt = mass;
        }
      }
      double tau = ct * 0.0333333;
      outFile<<bestMassForCt<<"\t"<<tau<<endl;
    }
    rMapInterpolated->Draw("colz");
    gridGraph->Draw("Psame");
  }
  
  outFile.close();
}
