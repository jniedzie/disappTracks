const double yAxisMin = 8E-2;
const double yAxisMax = 3.3;
const double xAxisMin = 300;
const double xAxisMax = 1000;

double logLabelSize = 0.10;
double labelSize = 0.045;
double titleSizeY = 0.12;
double titleSizeX = 0.050;

bool drawBySample      = false;
const string binning   = "3x3";
const string category  = "4-layers";
const string sampleTag = "_Chargino_";

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
//  {"atlas_observed"           , kRed     , 2, 1, "ATLAS JHEP 06 (2018) 022 (observed)"  },
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
//  {"2-tracks"                     , kBlue      , 2, 1, "2 tracks"  },
//  {"2-tracks_extended", kBlue      , 2, 1, "2 tracks"  },
//  {"2-tracks_with2016", kBlue      , 2, 2, "2 tracks, with 2016"  },
  
//  {"3x3_3-layers_notag_500_10"    , kGreen  , 2, 1, "3x3, 3 layers (old)"   },
//  {"3x3_4-layers_notag_1000_10"   , kRed    , 2, 1, "3x3, 4 layers (old)"   },
//  {"3x3_5-6-layers_notag_650_10"  , kBlue   , 2, 1, "3x3, 5-6 layers (old)" },
//  {"allcategories_notag"          , kCyan+1    , 2, 1, "all categories, no tagger"  },
//  {"allcategories_notag_run2" , kBlue     , 3, 2, "3x3, 3+4 layers, no tagger, Run 2" },
  
  
  {"new_interpolated"           , kGreen      , 2, 2, "CMG weights, exponential fit"      },
  {"new_extrapolated"           , kRed        , 2, 2, "CMG weights, x-sec extrapolation"  },
  
  {"new_interpolated_ohio"      , kGreen      , 2, 1, "EXO-19-010 weights, added ct=3 cm, exponential fit"      },
  {"new_interpolated_ohio"      , kGreen+2    , 2, 1, "EXO-19-010 weights, added ct=3 cm, added 600-10 and 800-30, exponential fit"      },
  {"new_extrapolated_ohio"      , kRed        , 2, 1, "EXO-19-010 weights, added ct=3 cm, x-sec extrapolation"  },
  {"new_root_interpolated_ohio" , kBlue       , 2, 1, "EXO-19-010 weights, added ct=3 cm, root interpolate"     },
  
  {"new_interpolate_ohio_fullRun2"  , kGreen  , 4, 3, "EXO-19-010 weights, added ct=3 cm, exponential fit, full Run 2" },
  {"extrapolated_ohio_fullRun2"     , kRed    , 4, 3, "EXO-19-010 weights, added ct=3 cm, x-sec extrapolation, full Run 2" },
  {"root_interpolate_ohio_fullRun2" , kBlue   , 4, 3, "EXO-19-010 weights, added ct=3 cm, root interpolate, full Run 2"     },
  
  
//  {"3x3_5-6-layers_800_10"      , kBlue    , 2, 1, "3x3, 5-6 layers, (new, int)"   },
//  {"3x3_4-layers_800_20"      , kRed     , 2, 1, "3x3, 4 layers, (new, int)"   },
//  {"3x3_3-layers_1000_20"     , kGreen   , 2, 1, "3x3, 3 layers, (new, int)"   },
//  {"3x3_4-layers_Chargino_500_1"    , kMagenta+2 , 2, 1, "3x3, 4 layers (new)"   },
//  {"3x3_4-layers_interpolated"      , kViolet+2  , 2, 1, "3x3, 4 layers (new, int)"   },
//  {"3x3_4-layers_extrapolated"      , kViolet+2  , 2, 2, "3x3, 4 layers (new, ext)"   },
//  {"3x3_5-6-layers_Chargino_900_10" , kOrange+2  , 2, 1, "3x3, 5-6 layers (new)" },
//  {"3x3_5-6-layers_interpolated"    , kCyan+2    , 2, 1, "3x3, 5-6 layers (new, int)" },
//  {"3x3_5-6-layers_extrapolated"    , kCyan+2    , 2, 2, "3x3, 5-6 layers (new, ext)" },
  {"allcategories_notag"          , kCyan+1    , 2, 1, "old"  },
//  {"allcategories_notag_run2" , kBlue     , 3, 2, "3x3, 3+4 layers, no tagger, Run 2" },
  
  
//  {"3x3_3-layers_noTag_1000_20_ext"   , kGreen     , 2, 1, "3x3, 3 layers, no tagger" },
//  {"3x3_4-layers_noTag_800_20_ext"    , kMagenta+2 , 2, 1, "3x3, 4 layers, no tagger"   },
//  {"3x3_5-6-layers_noTag_800_10_ext"  , kOrange+2  , 2, 1, "3x3, 5-6 layers, no tagger" },
//  {"3x3_noTag_ext"                    , kCyan+1    , 2, 1, "all categories, no tagger"  },
//
//  {"3x3_3-layers_noTag_1000_20"   , kGreen     , 2, 2, "3x3, 3 layers, no tagger" },
//  {"3x3_4-layers_noTag_800_20"    , kMagenta+2 , 2, 2, "3x3, 4 layers, no tagger"   },
//  {"3x3_5-6-layers_noTag_800_10"  , kOrange+2  , 2, 2, "3x3, 5-6 layers, no tagger" },
//  {"3x3_noTag"                    , kCyan+1    , 2, 2, "all categories, no tagger"  },
  
//  {"3x3_3-layers_tagSim_1000_20"   , kGreen     , 2, 2, "3x3, 3 layers, tagSim" },
//  {"3x3_4-layers_tagSim_800_20"    , kMagenta+2 , 2, 2, "3x3, 4 layers, tagSim"   },
//  {"3x3_5-6-layers_tagSim_800_10"  , kOrange+2  , 2, 2, "3x3, 5-6 layers, tagSim" },
//  {"3x3_tagSim"                    , kCyan+1    , 2, 2, "all categories, tagSim"  },
  
  
//  {"3x3_4-layers_LH_noTag_fix_800_20", kMagenta+2 , 2, 2, "3x3, 4 layers, no tagger, LH, fixed"},
  
  
//  {"3x3_3-layers_LH_notag_fix_1000_20_extended"   , kGreen  , 2, 1, "3x3, 3 layers (old)"   },
//  {"3x3_4-layers_LH_notag_fix_800_20_extended"    , kRed    , 2, 1, "3x3, 4 layers (old)"   },
//  {"3x3_5-6-layers_LH_notag_fix_800_10_extended"  , kBlue   , 2, 1, "3x3, 5-6 layers (old)" },
//  {"3x3_LH_notag_fix_extended"                    , kCyan+1    , 2, 2, "all categories, no tagger"  },
//
//  {"3x3_3-layers_LH_tagSimPU_fix_1000_20_extended"  , kGreen     , 2, 1, "3x3, 3 layers, tagSim PU"   },
//  {"3x3_4-layers_LH_tagSimPU_fix_800_20_extended"   , kMagenta+2 , 2, 1, "3x3, 4 layers, tagSim PU"   },
//  {"3x3_5-6-layers_LH_tagSimPU_fix_800_10_extended" , kOrange+2  , 2, 1, "3x3, 5-6 layers, tagSim PU" },
//  {"3x3_LH_tagSimPU_extended"                       , kCyan+1    , 2, 1, "all categories, tagSim PU"  },
  
//  {"3x3_3-layers_tagSim_noPU_500_10"    , kGreen     , 2, 2, "3 layers, tagger sim, no PU"    },
//  {"3x3_4-layers_tagSim_noPU_1000_10"   , kMagenta+2 , 2, 2, "4 layers, tagger sim, no PU"    },
//  {"3x3_5-6-layers_tagSim_noPU_650_10"  , kOrange+2  , 2, 2, "5-6 layers, tagger sim, no PU"  },
//  {"tagSim_noPU"                        , kCyan+1    , 2, 2, "all categories, tagger sim, no PU (no tagger for 4 layers"  },
};

vector<tuple<string, int, int>> graphParamsBySample = {
//  {"300_3"      , kMagenta  , 1 },
//  {"300_10"     , kMagenta  , 2 },
//  {"300_30"     , kMagenta  , 3 },
//  {"500_10"     , kCyan     , 1 },
//  {"500_20"     , kCyan     , 2 },
//  {"650_10"     , kMagenta+2, 1 },
//  {"650_20"     , kMagenta+2, 2 },
//  {"800_10"     , kCyan+2   , 1 },
//  {"800_20"     , kCyan+2   , 2 },
//  {"1000_10"    , kOrange   , 1 },
//  {"1000_20"    , kOrange   , 2 },
  
  {"300_1"      , kMagenta  , 1 },
  {"300_10"     , kMagenta+2, 1 },
  {"300_30"     , kMagenta+2, 2 },
  {"500_1"      , kCyan     , 1 },
  {"500_10"     , kCyan+2   , 1 },
  {"500_30"     , kCyan+2   , 2 },
  {"900_1"      , kGreen    , 1 },
  {"900_10"     , kGreen+2  , 1 },
  {"900_30"     , kGreen+2  , 2 },

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
  
  TLegend *leg = new TLegend(0.50, 0.15, 0.95, 0.5);
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
      string path = binning+"_"+category+sampleTag+sampleName;
      graphs[path] = GetGraphFromTxt(Form("%s/%s%s.txt",prefix.c_str(), thisAnalysisPrefix.c_str(), path.c_str()));
      
      TGraph *graph = graphs[path];
      graph->SetLineColor(color);
      graph->SetLineWidth(2);
      graph->SetLineStyle(style);
      graph->Draw(first ? "AL" : "Lsame");
      leg->AddEntry(graph, (binning+", "+category+", no tagger, "+sampleName).c_str(), "L");
      if(first){
        setFirstGraphOptions(graph);
        drawLines();
        first = false;
      }
    }
  }

  if(drawByBinning){
    for(auto &[path, color, style, title] : graphParamsByBinning){
      graphs[path] = GetGraphFromTxt(Form("%s/%s%s.txt",prefix.c_str(), thisAnalysisPrefix.c_str(), path.c_str()));
      
      TGraph *graph = graphs[path];
      graph->SetLineColor(color);
      graph->SetLineWidth(2);
      graph->SetLineStyle(style);
      graph->Draw(first ? "AL" : "Lsame");
      leg->AddEntry(graph, title.c_str(), "L");
      if(first){
        setFirstGraphOptions(graph);
        drawLines();
        first = false;
      }
    }
  }
  
  if(drawByCategory){
    for(auto &[path, color, width, style, title] : graphParamsByCategory){
      graphs[path] = GetGraphFromTxt(Form("%s/%s%s.txt",prefix.c_str(), thisAnalysisPrefix.c_str(), path.c_str()));
      
      TGraph *graph = graphs[path];
      graph->SetLineColor(color);
      graph->SetLineWidth(width);
      graph->SetLineStyle(style);
      graph->Draw(first ? "AL" : "Lsame");
      leg->AddEntry(graph, title.c_str(), "L");
      if(first){
        setFirstGraphOptions(graph);
        drawLines();
        first = false;
      }
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
