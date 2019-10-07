const double yAxisMin = 1E-2;
const double yAxisMax = 4E2;
const double xAxisMin = 50;
const double xAxisMax = 1000;

double logLabelSize = 0.10;
double labelSize = 0.045;
double titleSizeY = 0.12;
double titleSizeX = 0.050;

bool drawBySample = false;
const string binning = "3x3";
const string category = "5-6-layers";

bool drawByBinning = false;
bool drawByCategory = true;

bool sortbysec(const tuple<double, double>& a, const tuple<double, double>& b){ return (get<1>(a) < get<1>(b)); }

TGraph* GetGraphFromTxt(const char *fileName, double scale=1.0){
  
  cout<<"File:"<<fileName<<endl;
  ifstream inFile(fileName);
  TGraph *result = new TGraph();
  double x,y;
  int i=0;
  
  vector<tuple<double, double>> points;
  
  while (inFile >> x >> y){ points.push_back(make_tuple(x,scale*y)); }
  sort(points.begin(), points.end(), sortbysec);
  for(auto &[x, y] : points) result->SetPoint(i++,x,scale*y);
  
  inFile.close();
  return result;
}

vector<tuple<string, int, int, int, string>> otherGraphParams = {
// inFileName                                            color     width   style    title
  {"exo_16_044_observed"      , kBlack   , 2, 1, "CMS EXO-44-016 (observed)"            },
  {"exo_19_010_expected"      , kBlack   , 2, 2, "CMS EXO-19-010 (expected)"            },
  {"atlas_observed"           , kRed     , 2, 1, "ATLAS JHEP 06 (2018) 022 (observed)"  },
//  {"exo_16_044_expected"      , kBlack   , 2, 2, "EXO-44-016 (expected)"                },
//  {"exo_16_044_expected_68p"  , kGreen   , 1, 1, "EXO-44-016 (68% expected)"            },
};

vector<tuple<string, int, int, string>> graphParamsByBinning = { // best in each binning
// inFileName                                            color     width   style    title
  {"2x2_3-layers_notag_800_10"  , kCyan   , 1, "2x2, 3 layers, no tagger"},
  {"2x2_4-layers_notag_800_20"  , kCyan   , 2, "2x2, 4 layers, no tagger"},
  {"2x3_3-layers_notag_800_20"  , kMagenta, 1, "2x3, 3 layers, no tagger"},
  {"2x3_4-layers_notag_800_10"  , kMagenta, 2, "2x3, 4 layers, no tagger"},
  {"3x2_3-layers_notag_500_10"  , kBlue   , 1, "3x2, 3 layers, no tagger"},
  {"3x2_4-layers_notag_1000_10" , kBlue   , 2, "3x2, 4 layers, no tagger"},
  {"3x3_3-layers_notag_500_10"  , kOrange , 1, "3x3, 3 layers, no tagger"},
  {"3x3_4-layers_notag_1000_10" , kOrange , 2, "3x3, 4 layers, no tagger"},
};

vector<tuple<string, int, int, int, string>> graphParamsByCategory = { // best in each category
// inFileName                                            color     width   style    title
  {"2-tracks"                     , kBlue      , 2, 1, "2 tracks"  },
  {"3x3_3-layers_notag_500_10"    , kGreen     , 2, 1, "3x3, 3 layers, no tagger"   },
  {"3x3_4-layers_notag_1000_10"   , kMagenta+2 , 2, 1, "3x3, 4 layers, no tagger"   },
  {"3x3_5-6-layers_notag_650_10"  , kOrange+2  , 2, 1, "3x3, 5-6 layers, no tagger" },
  {"allcategories_notag"          , kCyan+1    , 2, 1, "all categories, no tagger"  },
//  {"allcategories_notag_run2" , kBlue     , 3, 2, "3x3, 3+4 layers, no tagger, Run 2" },
};

vector<tuple<string, int, int>> graphParamsBySample = {
  {"300_3"      , kMagenta  , 1 },
  {"300_10"     , kMagenta  , 2 },
  {"300_30"     , kMagenta  , 3 },
  {"500_10"     , kCyan     , 1 },
  {"500_20"     , kCyan     , 2 },
  {"650_10"     , kMagenta+2, 1 },
  {"650_20"     , kMagenta+2, 2 },
  {"800_10"     , kCyan+2   , 1 },
  {"800_20"     , kCyan+2   , 2 },
  {"1000_10"    , kOrange   , 1 },
  {"1000_20"    , kOrange   , 2 },
};

void setFirstGraphOptions(TGraph *graph)
{
  graph->GetXaxis()->SetTitle("m_{#chi} (GeV)");
  graph->GetXaxis()->SetTitleSize(titleSizeX);
  graph->GetXaxis()->SetTitleOffset(0.85);
  graph->GetXaxis()->SetLabelSize(labelSize);
  graph->GetXaxis()->SetLimits(xAxisMin, xAxisMax);
  
  graph->GetYaxis()->SetTitle("#tau (ns)");
  graph->GetYaxis()->SetTitleSize(titleSizeX);
  graph->GetYaxis()->SetTitleOffset(0.85);
  graph->GetYaxis()->SetLabelSize(labelSize);
  
  graph->GetYaxis()->SetRangeUser(yAxisMin,yAxisMax);
}

void drawLines()
{
  double min = yAxisMin/0.033333;
  double max = yAxisMax/0.033333;
  
  TGaxis *axis = new TGaxis(xAxisMax, yAxisMin,
                            xAxisMax, yAxisMax,
                            yAxisMin/0.033333, yAxisMax/0.033333, 510, "GL+");
  axis->SetLineColor(kBlack);
  axis->SetTitle("c#tau (cm)");
  axis->SetTitleOffset(1.1);
  axis->SetLabelColor(kBlack);
  axis->Draw();
  
  TLine *line1cm = new TLine(xAxisMin, 1*0.0333333, xAxisMax, 1*0.03333333);
  line1cm->SetLineColor(kBlack);
  line1cm->SetLineStyle(2);
  line1cm->SetLineWidth(1.0);
  line1cm->Draw("same");
  
  TLine *line3cm = new TLine(xAxisMin, 3*0.0333333, xAxisMax, 3*0.03333333);
  line3cm->SetLineColor(kBlack);
  line3cm->SetLineStyle(2);
  line3cm->SetLineWidth(1.0);
  line3cm->Draw("same");
  
  TLine *line30cm = new TLine(xAxisMin, 30*0.0333333, xAxisMax, 30*0.03333333);
  line30cm->SetLineColor(kBlack);
  line30cm->SetLineStyle(2);
  line30cm->SetLineWidth(1.0);
  line30cm->Draw("same");
}

void drawLimits()
{
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  
  map<string, TGraph*> graphs;
  string prefix = "limitsData";
  string thisAnalysisPrefix = "cms_short_disappearing_";
  
  TLegend *leg = new TLegend(0.10, 0.62, 0.50, 0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  
  gPad->SetLogy();
  
  bool first = true;
  for(auto &[path, color, width, style, title] : otherGraphParams){

    graphs[path] = GetGraphFromTxt(Form("%s/%s.txt",prefix.c_str(), path.c_str()));
    
    TGraph *graph = graphs[path];
    graph->SetLineColor(color);
    graph->SetLineWidth(width);
    graph->SetLineStyle(style);
    graph->Draw(first ? "AL" : "Lsame");
    
    if(first){
      setFirstGraphOptions(graph);
      drawLines();
      first = false;
    }

    leg->AddEntry(graph, title.c_str(), "L");
  }
  
  if(drawBySample){
    for(auto &[sampleName, color, style] : graphParamsBySample){
      string path = binning+"_"+category+"_notag_"+sampleName;
      graphs[path] = GetGraphFromTxt(Form("%s/%s%s.txt",prefix.c_str(), thisAnalysisPrefix.c_str(), path.c_str()));
      
      TGraph *graph = graphs[path];
      graph->SetLineColor(color);
      graph->SetLineWidth(2);
      graph->SetLineStyle(style);
      graph->Draw("Lsame");
      leg->AddEntry(graph, (binning+", "+category+", no tagger, "+sampleName).c_str(), "L");
    }
  }

  if(drawByBinning){
    for(auto &[path, color, style, title] : graphParamsByBinning){
      graphs[path] = GetGraphFromTxt(Form("%s/%s%s.txt",prefix.c_str(), thisAnalysisPrefix.c_str(), path.c_str()));
      
      TGraph *graph = graphs[path];
      graph->SetLineColor(color);
      graph->SetLineWidth(2);
      graph->SetLineStyle(style);
      graph->Draw("Lsame");
      leg->AddEntry(graph, title.c_str(), "L");
    }
  }
  
  if(drawByCategory){
    for(auto &[path, color, width, style, title] : graphParamsByCategory){
      graphs[path] = GetGraphFromTxt(Form("%s/%s%s.txt",prefix.c_str(), thisAnalysisPrefix.c_str(), path.c_str()));
      
      TGraph *graph = graphs[path];
      graph->SetLineColor(color);
      graph->SetLineWidth(width);
      graph->SetLineStyle(style);
      graph->Draw("Lsame");
      leg->AddEntry(graph, title.c_str(), "L");
    }
  }
  
  leg->Draw("same");
  
  TLegend *leg_cms_logo = new TLegend(0.10,0.71,0.20,1.15);
  leg_cms_logo->SetTextColor(kBlack);
  leg_cms_logo->SetFillStyle(0);
  leg_cms_logo->SetBorderSize(0);
  leg_cms_logo->SetTextFont(62);
  leg_cms_logo->AddEntry(new TGraph(),"CMS","");
  leg_cms_logo->Draw();
  
}
