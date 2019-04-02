#include "EventSet.hpp"
#include "Helpers.hpp"
#include "ConfigManager.hpp"

#include <TApplication.h>

string configPath = "configs/analysis.md";

void SetupDefaultCuts(unique_ptr<EventCut>  &eventCut,
                      unique_ptr<TrackCut>  &trackCut,
                      unique_ptr<JetCut>    &jetCut,
                      unique_ptr<LeptonCut> &leptonCut);

int main(int argc, char* argv[])
{
  TApplication *theApp = new TApplication("App", &argc, argv);
  
  // All events with initial cuts only
  config = make_unique<ConfigManager>(configPath);
  auto events = make_shared<EventSet>();
  events->LoadEventsFromFiles("");
  
  auto eventCut  = make_unique<EventCut>();
  auto trackCut  = make_unique<TrackCut>();
  auto jetCut    = make_unique<JetCut>();
  auto leptonCut = make_unique<LeptonCut>();
  
  SetupDefaultCuts(eventCut, trackCut, jetCut, leptonCut);
  
  // choose category:
  eventCut->SetNtracks(range<int>(2,2));
//  events->ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  
  cout<<"\n\nBefore cuts:"<<endl;
  events->PrintYields();

  // on top of default cuts, apply:
//  trackCut->SetNlayers(range<int>(0,10));
//  events->ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  
//  cout<<"\n\nAfter cuts:"<<endl;
//  events->PrintYields();
  
  if(config->drawStandardPlots) events->DrawStandardPlots();
  if(config->drawPerLayerPlots) events->DrawPerLayerPlots();
 
  theApp->Run();
  return 0;
}



void SetupDefaultCuts(unique_ptr<EventCut>  &eventCut,
                      unique_ptr<TrackCut>  &trackCut,
                      unique_ptr<JetCut>    &jetCut,
                      unique_ptr<LeptonCut> &leptonCut)
{
  // Remove bad jets
  jetCut->SetChargedHadronEnergyFraction(range<double>(0.01,0.99));
  jetCut->SetNeutralHadronEnergyFraction(range<double>(0.01,0.99));
  jetCut->SetPt(range<double>(30.0, inf));
  
  // Remove bad tracks
  trackCut->SetNmissingInnerPixel(range<int>(0, 0));
  trackCut->SetNmissingMiddleTracker(range<int>(0, 0));
  trackCut->SetRelativeIsolation(range<double>(0.0, 0.5));
  trackCut->SetNlayers(range<int>(2, inf));
  trackCut->SetEta(range<double>(-2.1, 2.1));
  
  // Check MET properties
  eventCut->SetMetNoMuPt(range<double>(200,inf));
  eventCut->SetRequireMetNoMuTrigger(true);
  eventCut->SetRequirePassingAllFilters(true);
  eventCut->SetJetMetDeltaPhi(range<double>(0.5,inf));
  
  // Check leading jet properties
  eventCut->SetLeadingJetPt(range<double>(100,inf));
  eventCut->SetLeadingJetEta(range<double>(-2.4,2.4));
  eventCut->SetLeadingJetNeHEF(range<double>(-inf,0.8));
  eventCut->SetLeadingJetChHEF(range<double>(0.1,inf));
  
  // Check number of objects after cuts
  eventCut->SetNtracks(range<int>(1,inf));
  eventCut->SetNjets(range<int>(1,inf));
  eventCut->SetNmuons(range<int>(0,0));
  eventCut->SetNtaus(range<int>(0,0));
  eventCut->SetNleptons(range<int>(0,0));
}
