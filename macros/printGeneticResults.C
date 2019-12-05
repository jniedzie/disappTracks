#include <fstream>
#include <sstream>
#include <string>

vector<vector<string>> getIterations(string inPath)
{
  ifstream infile(inPath);

  
  string line;
  
  vector<vector<string>> allIterations;
  
  vector<string> thisIteration;
  
  while(getline(infile, line)){
    
    if(line.find("Starting iteration") != string::npos){
      allIterations.push_back(thisIteration);
      thisIteration.clear();
    }
    thisIteration.push_back(line);
  }
  return allIterations;
}

vector<string> cleanIteration(vector<string> input)
{
  vector<string> output;
  
  for(string line : input){
    if(line.find("======") != string::npos) continue;
    if(line.find("------") != string::npos) continue;
    if(line.find("Starting iteration") != string::npos) continue;
    
    output.push_back(line);
  }
  
  return output;
}

vector<double> getTimesInIteration(vector<string> iteration)
{
  vector<double> times;
  for(string line : iteration){
    if(line.find("Time per event") == string::npos) continue;
    istringstream iss(line);
    string trash;
    double time;
    if (!(iss >> trash >> trash >> trash >> time >> trash)) { cout<<"FAIL!"<<endl; break; }
    times.push_back(time);
  }
  return times;
}

vector<double> getDistancesInIteration(vector<string> iteration)
{
  vector<double> distances;
  for(string line : iteration){
    if(line.find("Distance to") == string::npos) continue;
    istringstream iss(line);
    string trash;
    double distance;
    if (!(iss >> trash >> trash >> trash >> distance)) { cout<<"FAIL!"<<endl; break; }
    distances.push_back(distance);
  }
  return distances;
}

double getBestDistanceInIteration(vector<string> iteration)
{
  vector<double> distances;
  for(string line : iteration){
    if(line.find("fitness") == string::npos) continue;
    istringstream iss(line);
    string trash;
    double fitness;
    if (!(iss >> trash >> trash >> fitness >> trash)) { cout<<"FAIL!"<<endl; break; }
    return 1-fitness;
  }
  return -2;
}

void printGeneticResults(string inPath)
{
  vector<vector<string>> textByIteration = getIterations(inPath);
  
  int nIterations = textByIteration.size()-1;
  
  TH2D *execTimeHist = new TH2D("execTimeHist", "execTimeHist", nIterations, 0, nIterations, 30, 0, 3.0);
  TH2D *distanceHist = new TH2D("distanceHist", "distanceHist", nIterations, 0, nIterations, 50, -0.2, 0.6);
  TGraph *distanceEvolution = new TGraph();
  distanceEvolution->SetMarkerStyle(20);
  distanceEvolution->SetMarkerSize(1.0);
  distanceEvolution->SetMarkerColor(kRed);
  
  int iGeneration=0;
  bool first = true;
  
  for(vector<string> thisIteration : textByIteration){
    if(first){ first = false; continue; }
//    cout<<"Iter "<<iGeneration<<endl;
    
    thisIteration = cleanIteration(thisIteration);
    
    vector<double> times = getTimesInIteration(thisIteration);
    vector<double> distances = getDistancesInIteration(thisIteration);
    double bestDistance = getBestDistanceInIteration(thisIteration);
    
    for(double time : times) execTimeHist->Fill(iGeneration, time);
    for(double distance : distances) distanceHist->Fill(iGeneration, distance);
    distanceEvolution->SetPoint(iGeneration, iGeneration, bestDistance);
    
//    for(string line : thisIteration) cout<<line<<endl;
    
    iGeneration++;
  }

  bool bestStarted = false;
  cout<<"Best parameters in the latest iteration:"<<endl;
  for(string line : textByIteration.back()){
    if(bestStarted) cout<<line<<endl;
    if(line.find("BEST PARAMS") != string::npos) bestStarted = true;
  }
  
  TCanvas *c1 = new TCanvas("c1","c1", 1000, 1000);
  c1->Divide(2,2);
  
  c1->cd(1); execTimeHist->Draw("colz");
  c1->cd(2); distanceHist->Draw("colz");
  c1->cd(3); distanceEvolution->Draw("AP");
}
