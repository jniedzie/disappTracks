// ct in cm (*0.03333 to convert to ns)
const double ctMin = 0;
const double ctMax = 35;
const double ctStep = 0.1;

// mass in GeV
const double massMin = 0;
const double massMax = 1100;
const double massStep = 1;

vector<tuple<double, double, double>> getTupleFromFile(string fileName){
  vector<tuple<double, double, double>> rValues;
  int mass, ct;
  double r;
    
  ifstream inFile(fileName);
  while (inFile >> mass >> ct >> r){ rValues.push_back(make_tuple(mass, ct, r)); }
  inFile.close();
  
  return rValues;
}

void getLimitsFromR(string inputPath, string outPath)
{
  auto rValues3x3 = getTupleFromFile(inputPath);
  
  TH2D *rMap = new TH2D("rMap", "rMap",
                        (massMax-massMin)/massStep, massMin, massMax,
                        (ctMax-ctMin)/ctStep+1, ctMin, ctMax);
  
  TGraph *gridGraph = new TGraph();
  gridGraph->SetMarkerStyle(4);
  gridGraph->SetMarkerSize(1.0);
  gridGraph->SetMarkerColor(kRed);
  
  int iPoint=0;
  for(auto &[mass, ct, r] : rValues3x3){
    if(r==0) r = 1E-10;
    rMap->Fill(mass, ct, r);
    gridGraph->SetPoint(iPoint++, mass, ct);
  }
  
  TGraph2D *rGraph3x3 = new TGraph2D(rMap);
  TH2D *rMapInterpolated = new TH2D(*rMap);
  rMapInterpolated->Reset();
  
  ofstream outFile(outPath);
  
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
  outFile.close();
  
  rMapInterpolated->Draw("colz");
  gridGraph->Draw("Psame");
}
