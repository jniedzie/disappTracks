//  analyzeClusters.cpp
//
//  Created by Jeremi Niedziela on 10/05/2019.

#include "Helpers.hpp"
#include "EventSet.hpp"

string configPath = "configs/eventDisplay.md";
string cutLevel = "after_L1/3layers/";//after_L1/";

xtracks::EDataType dataType = xtracks::kSignal;
int setIter = kWino_M_300_cTau_10;
int iEvent = 3;

struct ComparePointByZ{
  bool operator() (const shared_ptr<Point> &p1, const shared_ptr<Point> &p2){
    return (p1->GetZ() < p2->GetZ());
  }
};

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
 
  config = ConfigManager(configPath);
  
  EventSet events;
  events.LoadEventsFromFiles(cutLevel);
  
  //  auto event = events.GetEvent(dataType, searchRun, searchLumi, searchEvent);
  
  for(int iEvent=0; iEvent<events.size(dataType, setIter); iEvent++){
    auto event = events.At(dataType, setIter, iEvent);
    
    if(!event){
      cout<<"eventDisplay -- event not found"<<endl;
      exit(0);
    }
    
    event->LoadAdditionalInfo();
    
    auto eventVertex = event->GetVertex();
    auto hits = event->GetPionSimHits();
    sort(hits.begin(), hits.end(), ComparePointByZ());
    
    shared_ptr<Point> farthestPoint;
    if(fabs(hits.front()->GetZ()) > fabs(hits.back()->GetZ()))  farthestPoint = hits.front();
    else                                                        farthestPoint = hits.back();
    
    double zRange = fabs(eventVertex->GetZ() - farthestPoint->GetZ());
    
    cout<<"Z range:"<<zRange<<endl;
    cout<<"Last point in "<<farthestPoint->GetSubDetName()<<endl;
  }
  
  theApp.Run();
}
