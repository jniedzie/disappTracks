
const double yAxisMin = 1E-2;
const double yAxisMax = 4E2;
const double xAxisMin = 50;
const double xAxisMax = 1000;

double logLabelSize = 0.10;
double labelSize = 0.045;
double titleSizeY = 0.12;
double titleSizeX = 0.050;

bool sortbysec(const tuple<double, double>& a,
               const tuple<double, double>& b)
{
    return (get<1>(a) < get<1>(b));
}

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

vector<tuple<string, int, int, int, double, string, string>> graphParams = {
// inFileName                                   color     width   style   opacity   options title
  {"exo_16_044_observed"                        , kBlack   , 2, 1, 0.4, "L", "CMS EXO-44-016 (observed)"                     },
  {"atlas_observed"                             , kRed     , 2, 1, 0.4, "L", "ATLAS JHEP 06 (2018) 022 (observed)"           },
  {"cms_short_disappearing_2x2_nocat_notag"     , kViolet+1, 3, 2, 0.4, "L", "CMS (2x2 histogram, no categories, no tagger)" },
  {"cms_short_disappearing_3x3_nocat_notag"     , kViolet+1, 3, 1, 0.4, "L", "CMS (3x3 histogram, no categories, no tagger)" },
  {"cms_short_disappearing_3x3_3layers_notag"   , kGreen+1, 3, 1, 0.4, "L", "CMS (3x3 histogram, 3 layers, no tagger)"       },
  {"cms_short_disappearing_3x3_4layers_notag"   , kGreen+1, 3, 2, 0.4, "L", "CMS (3x3 histogram, 4 layers, no tagger)"       },
  {"cms_short_disappearing_3x3_3+4layers_notag" , kGreen+1, 3, 3, 0.4, "L", "CMS (3x3 histogram, 3+4 layers, no tagger)"     },
//  {"exo_16_044_expected"    , kBlack    , 2     , 2  , 0.4   , "L", "EXO-44-016 (expected)" },
//  {"exo_16_044_expected_68p", kGreen   , 1     , 1  , 1.0   , "L*", "EXO-44-016 (68% expected)" },
};

void drawLimits()
{
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  
  map<string, TGraph*> graphs;
  string prefix = "limitsData";
  
  TLegend *leg = new TLegend(0.10, 0.62, 0.50, 0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  
  gPad->SetLogy();
  
  bool first = true;
  for(auto &[path, color, width, style, opacity, options, title] : graphParams){

    graphs[path] = GetGraphFromTxt(Form("%s/%s.txt",prefix.c_str(), path.c_str()));
    
    TGraph *graph = graphs[path];
    graph->SetLineColor(color);
    graph->SetLineWidth(width);
    graph->SetLineStyle(style);
    graph->SetFillColorAlpha(color,opacity);
  
    graph->Draw(first ? ("A"+options).c_str() : (options+"same").c_str());
      
    if(first){
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
      
      double min = yAxisMin/0.033333;
      double max = yAxisMax/0.033333;
      
      cout<<"Min: "<<min<<"\tmax: "<<max<<endl;
      
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
      
      
      first = false;
    }

    leg->AddEntry(graph, title.c_str(), options.c_str());
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
