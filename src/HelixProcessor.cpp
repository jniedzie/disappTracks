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
  
  if(fabs(h1.GetOrigin().GetX() - h2.GetOrigin().GetX()) > config.params["tolerance_x"]) reasons.push_back(1);
  if(fabs(h1.GetOrigin().GetY() - h2.GetOrigin().GetY()) > config.params["tolerance_y"]) reasons.push_back(2);
  if(fabs(h1.GetOrigin().GetZ() - h2.GetOrigin().GetZ()) > config.params["tolerance_z"]) reasons.push_back(3);
  if(fabs(h1.GetMomentum().GetX() - h2.GetMomentum().GetX()) > config.params["tolerance_px"]) reasons.push_back(4);
  if(fabs(h1.GetMomentum().GetY() - h2.GetMomentum().GetY()) > config.params["tolerance_py"]) reasons.push_back(5);
  if(fabs(h1.GetMomentum().GetZ() - h2.GetMomentum().GetZ()) > config.params["tolerance_pz"]) reasons.push_back(6);
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
    
    auto momentum = Point(arrayValuesFloat["helix_px"][iHelix],
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
    arrayValuesFloat["helix_px"][iHelix]   = helices[iHelix]->GetMomentum().GetX();
    arrayValuesFloat["helix_py"][iHelix]   = helices[iHelix]->GetMomentum().GetY();
    arrayValuesFloat["helix_pz"][iHelix]   = helices[iHelix]->GetMomentum().GetZ();
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
  // here taking any element, probably doesn't matter which one
  size_t lastPointIndex = helix.GetLastPointsIndices().front();
  
  double Rh = helix.GetRadius(helix.GetPointT(lastPointIndex));
  double Rl = layerRanges[layerIndex].GetMin();
  
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

bool HelixProcessor::IsPointCloseToHelixInLayer(const Helix &helix, const Point &point, int layer,
                                                bool closeToPoint)
{
  auto lastPoints = helix.GetLastPoints();
  
  bool goodXY=false;
  Point pA, pB;
  GetIntersectionWithLayer(helix, layer, pA, pB);
  
  vector<Point> newPointsEstimated;
  
  for(auto &lastPoint : lastPoints){
    
    bool pAcloser = pointsProcessor.distanceXY(pA, *lastPoint) < pointsProcessor.distanceXY(pB, *lastPoint);
    
    if(closeToPoint){
      if(pAcloser) newPointsEstimated.push_back(pA);
      else         newPointsEstimated.push_back(pB);
    }
    else{
      if(pAcloser) newPointsEstimated.push_back(pB);
      else         newPointsEstimated.push_back(pA);
    }
  }
  
  for(Point newPointEstimeted : newPointsEstimated){
    double distanceXYtoPoint = pointsProcessor.distanceXY(newPointEstimeted, point);
    
    if(distanceXYtoPoint <= config.params["next_point_max_delta_xy"]){
      goodXY=true;
      break;
    }
  }
  return goodXY;
}

shared_ptr<Point> HelixProcessor::GetPointCloseToHelixInLayer(const Helix &helix, int layer)
{
  Point pA, pB;
  GetIntersectionWithLayer(helix, layer, pA, pB);
  
  auto lastPoint = helix.GetLastPoints().front();
  shared_ptr<Point> newPoint;
  
  bool pAcloser = pointsProcessor.distanceXY(pA, *lastPoint) < pointsProcessor.distanceXY(pB, *lastPoint);
  
  if(pAcloser) newPoint = make_shared<Point>(pA);
  else         newPoint = make_shared<Point>(pB);
  
  newPoint->SetLayer(layer);
  
  return newPoint;
}

size_t HelixProcessor::GetNcommonPoints(const Helix &helix1, const Helix &helix2)
{
  Points commonPoints;
  auto points1 = helix1.GetPoints();
  auto points2 = helix2.GetPoints();
  
  for(auto &p1 : points1){
    for(auto &p2 : points2){
      if(*p1 == *p2) commonPoints.push_back(p1);
    }
  }
  
  return commonPoints.size();
}

double HelixProcessor::GetAvgNhits(Helices helices)
{
  if(helices.size()==0) return 0;
  double avgHits = 0;
  for(auto helix : helices) avgHits += helix.GetNpoints();
  avgHits /= helices.size();
  return avgHits;
}

int HelixProcessor::GetMaxNhits(Helices helices)
{
  int maxNhits = 0;
  for(auto helix : helices){
    if((int)helix.GetNpoints() > maxNhits) maxNhits = (int)helix.GetNpoints();
  }
  return maxNhits;
}

int HelixProcessor::GetAvgNlayers(Helices helices)
{
  double avgLayers = 0;
  for(auto helix : helices) avgLayers += helix.GetNlayers();
  avgLayers /= helices.size();
  return avgLayers;
}

int HelixProcessor::GetMaxNlayers(Helices helices)
{
  size_t maxNlayers = 0;
  
  for(auto helix : helices){
    size_t nLayers = helix.GetNlayers();
    if(nLayers > maxNlayers) maxNlayers = nLayers;
  }
  return (int)maxNlayers;
}

double HelixProcessor::GetAvgLength(Helices helices)
{
  if(helices.size()==0) return 0;
  double avgLength = 0;
  
  for(auto helix : helices){
    double length = fabs(helix.GetTmax() - helix.GetTmin());
    avgLength += length;
  }
  avgLength /= helices.size();
  return avgLength;
}

double HelixProcessor::GetMaxLength(Helices helices)
{
  if(helices.size()==0) return 0;
  double maxLength = -inf;
  
  for(auto helix : helices){
    double length = fabs(helix.GetTmax() - helix.GetTmin());
    if(length > maxLength) maxLength = length;
  }
  return maxLength;
}

double HelixProcessor::GetMinChi2(Helices helices)
{
  if(helices.size()==0) return 0;
  double minChi2 = inf;
  
  for(auto helix : helices){
    if(helix.GetChi2() < minChi2) minChi2 = helix.GetChi2();
  }
  return minChi2;
}

double HelixProcessor::GetMinChi2overNhits(Helices helices)
{
  if(helices.size()==0) return 0;
  double minChi2 = inf;
  
  for(auto helix : helices){
    if(helix.GetChi2()/helix.GetNpoints() < minChi2) minChi2 = helix.GetChi2()/helix.GetNpoints();
  }
  
  return minChi2;
}

bool HelixProcessor::DidTurnBack(Helices helices)
{
  if(helices.size()==0) return false;
  
  for(auto helix : helices){
    if(helix.GetFirstTurningPointIndex()>0) return true;
  }
  
  return false;
}
