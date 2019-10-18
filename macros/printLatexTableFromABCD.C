
string histName = "abcd_plots_3x3_4-layers_tagSim_noPU_1000_10.root";

vector<tuple<int,int>> samples = {
  {300, 3}, {300, 10}, {300, 30},
  {500, 10}, {500, 20},
  {650, 10}, {650, 20},
  {800, 10}, {800, 20},
  {1000, 10}, {1000, 20},
};

template <typename T>
string to_string_with_precision(const T a_value, const int n = 6)
{
  ostringstream out;
  out.precision(n);
  out << fixed << a_value;
  return out.str();
}

string getTableBeginning(int nColumns, vector<tuple<double, double>> rangesX)
{
  string text = "\\begin{table}[t!]\n"
  "\\centering\n"
  "\\begin{tabular}{l";
  for(int i=0; i<nColumns; i++) text += "|l";
  text += "}\n"
  "\\hline\n"
  "dE/dx -- MET ";
  
  for(auto &[down, up] : rangesX) text += "& "+to_string_with_precision(down, 1) + " - " + to_string_with_precision(up, 1);
  text += "  \\\\\n";
  
  return text;
}

string getHeaderForSample(int mass, int ct, int nColumns){
  return "\\hline \\multicolumn{"+to_string(nColumns+1)+"}{c}{Mass: "+to_string(mass)+" GeV,~c$\\tau$: "+to_string(ct)+"~cm} \\\\ \\hline";
}

string getTableEnding()
{
  return "\\end{tabular}\n"
  "\\caption{Number of events for each signal sample for each value of dE/dx and MET p_{t}.}\n"
  "\\label{tab:n_events}\n"
  "\\end{table}\n";
}

void printLatexTableFromABCD()
{
  TFile *inFile = TFile::Open(("../results/"+histName).c_str());
  
    cout<<"\n\nCopy this to latex:\n\n"<<endl;
  
  TH2D *hist = (TH2D*)inFile->Get("background");
  int nColumns = hist->GetNbinsX();
  
  vector<tuple<double, double>> rangesX;
  for(int binX=1; binX<=nColumns; binX++){
    rangesX.push_back({hist->GetXaxis()->GetBinLowEdge(binX), hist->GetXaxis()->GetBinUpEdge(binX)});
  }
  
  cout<<getTableBeginning(nColumns, rangesX);
  
  for(auto &[mass, ct] : samples){
    hist = (TH2D*)inFile->Get(("Wino_m"+to_string(mass)+"_ct"+to_string(ct)).c_str());
    
    cout<<getHeaderForSample(mass, ct, nColumns)<<endl;
    
//    cout<<endl;
    
    for(int binY=1; binY<=hist->GetNbinsY(); binY++){
      for(int binX=1; binX<=hist->GetNbinsX(); binX++){
        if(binX==1) cout<<"\n"<<hist->GetYaxis()->GetBinLowEdge(binY)<<" - "<<hist->GetYaxis()->GetBinUpEdge(binY)<<" ";

        
        double nEvents = hist->GetBinContent(binX, binY);
        cout<<"\t&\t"<<to_string_with_precision(nEvents, 1);
      }
      cout<<"\\\\";
    }
    cout<<endl;
  }
  
  cout<<getTableEnding()<<endl;
  
}




