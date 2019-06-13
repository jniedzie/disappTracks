//
//  HelixProcessor.cpp
//
//  Created by Jeremi Niedziela on 21/01/2019.
//

#include "HelixProcessor.hpp"

HelixProcessor helixProcessor = HelixProcessor();

HelixProcessor::HelixProcessor()
{
  arrayNamesFloat = {
    "helix_x",
    "helix_y",
    "helix_z",
    "helix_px",
    "helix_py",
    "helix_pz"
  };
  
  arrayNamesInt = {
    "helix_charge"
  };
}

HelixProcessor::~HelixProcessor()
{
  
}

vector<int> HelixProcessor::AreIdentical(const Helix &h1, const Helix &h2)
{
  vector<int> reasons;
  
  if(fabs(h1.GetOrigin().GetX() - h2.GetOrigin().GetX()) > config.toleranceX) reasons.push_back(1);
  if(fabs(h1.GetOrigin().GetY() - h2.GetOrigin().GetY()) > config.toleranceY) reasons.push_back(2);
  if(fabs(h1.GetOrigin().GetZ() - h2.GetOrigin().GetZ()) > config.toleranceZ) reasons.push_back(3);
  if(fabs(h1.GetMomentum()->GetX() - h2.GetMomentum()->GetX()) > config.tolerancePx) reasons.push_back(4);
  if(fabs(h1.GetMomentum()->GetY() - h2.GetMomentum()->GetY()) > config.tolerancePy) reasons.push_back(5);
  if(fabs(h1.GetMomentum()->GetZ() - h2.GetMomentum()->GetZ()) > config.tolerancePz) reasons.push_back(6);
  if(h1.GetCharge() != h2.GetCharge()) reasons.push_back(7);
  
  return reasons;
}

vector<shared_ptr<Helix>> HelixProcessor::GetHelicesFromTree()
{
  vector<shared_ptr<Helix>> helices = vector<shared_ptr<Helix>>();
  
  for(int iHelix=0;iHelix<nHelices;iHelix++){
    auto origin   = Point(arrayValuesFloat["helix_x"][iHelix],
                          arrayValuesFloat["helix_y"][iHelix],
                          arrayValuesFloat["helix_z"][iHelix]);
    
    auto momentum = make_unique<Point>(arrayValuesFloat["helix_px"][iHelix],
                                       arrayValuesFloat["helix_py"][iHelix],
                                       arrayValuesFloat["helix_pz"][iHelix]);
    
    auto helix    = make_shared<Helix>(origin, momentum, arrayValuesInt["helix_charge"][iHelix]);
    
    helices.push_back(helix);
  }

  return helices;
}

void HelixProcessor::SaveHelicesToTree(vector<shared_ptr<Helix>> helices)
{
  nHelices = (int)helices.size();
  
  for(int iHelix=0;iHelix<nHelices;iHelix++){
    if(!helices[iHelix]) continue;
    
    arrayValuesFloat["helix_x"][iHelix]    = helices[iHelix]->GetOrigin().GetX();
    arrayValuesFloat["helix_y"][iHelix]    = helices[iHelix]->GetOrigin().GetY();
    arrayValuesFloat["helix_z"][iHelix]    = helices[iHelix]->GetOrigin().GetZ();
    arrayValuesFloat["helix_px"][iHelix]   = helices[iHelix]->GetMomentum()->GetX();
    arrayValuesFloat["helix_py"][iHelix]   = helices[iHelix]->GetMomentum()->GetY();
    arrayValuesFloat["helix_pz"][iHelix]   = helices[iHelix]->GetMomentum()->GetZ();
    arrayValuesInt["helix_charge"][iHelix] = helices[iHelix]->GetCharge();
  }
}

void HelixProcessor::SetupBranchesForReading(TTree *tree)
{
  // check if there is a branch with helices in the tree
  if(!tree->GetBranchStatus("nFittedHelices")){
    nHelices = 0;
    return;
  }
  else{
    tree->SetBranchAddress("nFittedHelices", &nHelices);
  }
  
  for(string name : arrayNamesFloat){
    if(!tree->GetBranchStatus(name.c_str())){
      cout<<"WARNING -- no branch named "<<name<<"!!"<<endl;
      continue;
    }
    tree->SetBranchAddress(name.c_str(), &arrayValuesFloat[name]);
  }
  
  for(string name : arrayNamesInt){
    if(!tree->GetBranchStatus(name.c_str())){
      cout<<"WARNING -- no branch named "<<name<<"!!"<<endl;
      continue;
    }
    tree->SetBranchAddress(name.c_str(), &arrayValuesInt[name]);
  }
}

void HelixProcessor::SetupBranchesForWriting(TTree *tree)
{
  tree->Branch("nFittedHelices", &nHelices, "nFittedHelices/I");
  
  for(string name : arrayNamesFloat){
    tree->Branch(name.c_str(), &arrayValuesFloat[name], Form("%s[nFittedHelices]/F", name.c_str()));
  }
  
  for(string name : arrayNamesInt){
    tree->Branch(name.c_str(), &arrayValuesInt[name], Form("%s[nFittedHelices]/I", name.c_str()));
  }
}

bool HelixProcessor::GetIntersectionWithLayer(const Helix &helix, int layerIndex, Point &pA, Point &pB)
{
  auto p1 = helix.GetLastPoints().front(); // here taking any element, probably doesn't matter which one
  
  double Rh = helix.GetRadius(p1->GetT());
  double Rl = layerR[layerIndex];
  
  double x0 = helix.GetOrigin().GetX();
  double y0 = helix.GetOrigin().GetY();
  
  double R = sqrt(x0*x0+y0*y0);
  
  double x_a = x0/2 + (Rl*Rl-Rh*Rh)/(2*R*R)*x0 + y0/2*sqrt(2*(Rl*Rl+Rh*Rh)/(R*R)-pow((Rl*Rl-Rh*Rh)/(R*R), 2)-1);
  double y_a = y0/2 + (Rl*Rl-Rh*Rh)/(2*R*R)*y0 - x0/2*sqrt(2*(Rl*Rl+Rh*Rh)/(R*R)-pow((Rl*Rl-Rh*Rh)/(R*R), 2)-1);
  
  pA = Point(x_a, y_a, 0);
  
  double x_b = x0/2 + (Rl*Rl-Rh*Rh)/(2*R*R)*x0 - y0/2*sqrt(2*(Rl*Rl+Rh*Rh)/(R*R)-pow((Rl*Rl-Rh*Rh)/(R*R), 2)-1);
  double y_b = y0/2 + (Rl*Rl-Rh*Rh)/(2*R*R)*y0 + x0/2*sqrt(2*(Rl*Rl+Rh*Rh)/(R*R)-pow((Rl*Rl-Rh*Rh)/(R*R), 2)-1);
  
  pB = Point(x_b, y_b, 0);
  
  if(!isnormal(x_a) || !isnormal(x_b) || !isnormal(y_a) || !isnormal(y_b)) return false;
  
  return true;
}

bool HelixProcessor::IsPointCloseToHelixInLayer(const Helix &helix, const Point &point, int layer)
{
  auto lastPoints         = helix.GetLastPoints();
  auto secondToLastPoints = helix.GetSecontToLastPoints();
  
  bool goodXY=false;
  Point pA, pB;
  GetIntersectionWithLayer(helix, layer, pA, pB);
  
  vector<Point> newPointsEstimated;
  
  for(auto &secondToLastPoint : secondToLastPoints){
    
    if(pointsProcessor.distanceXY(pA, *secondToLastPoint) <
       pointsProcessor.distanceXY(pB, *secondToLastPoint)) newPointsEstimated.push_back(pB);
    else                                                   newPointsEstimated.push_back(pA);
    
  }
  
  for(Point newPointEstimeted : newPointsEstimated){
    double distanceXYtoPoint = pointsProcessor.distanceXY(newPointEstimeted, point);
    if(distanceXYtoPoint <= 150){
      goodXY=true;
      break;
    }
  }
  return goodXY;
}
