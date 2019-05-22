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

vector<shared_ptr<Point>> HelixProcessor::GetPointsHittingSilicon(const Helix &helix)
{
  // this method may reaquire updating after sign of the pion charge changed...
  vector<shared_ptr<Point>> points;
  unique_ptr<Point> origin = make_unique<Point>(helix.origin);
  
  double dh = sqrt(pow(origin->GetX(),2)+pow(origin->GetY(),2));
  double Rl, C, delta;
  double x1,y1,x2,y2,z1,z2,t1,t2;
  
  for(int iLayer=0;iLayer<config.nTrackerLayers;iLayer++){
    Rl = layerR[iLayer];
    C = (Rl*Rl+dh*dh-helix.radius*helix.radius)/2.;
    
    delta = 4*origin->GetY()*origin->GetY()*C*C-4*dh*dh*(C*C-Rl*Rl*origin->GetX()*origin->GetX());
    if(delta < 0) continue;
    
    y1 = (2*origin->GetY()*C+sqrt(delta))/(2*dh*dh);
    y2 = (2*origin->GetY()*C-sqrt(delta))/(2*dh*dh);
    
    x1 = (C-y1*origin->GetY())/origin->GetX();
    x2 = (C-y2*origin->GetY())/origin->GetX();
    
    t1 = atan2((y1-origin->GetY()),(x1-origin->GetX()));
    t2 = atan2(y2-origin->GetY(),x2-origin->GetX());
    
    if(helix.charge < 0){
      t1 = TMath::Pi()/2. - t1;
      t2 = TMath::Pi()/2. - t2;
    }
    
    double nCycles = fabs(helix.GetNcycles());
    int signZ = sgn(helix.momentum->GetZ());
    double tShift = helix.tShift;
    double slopeAbs = helix.slopeAbs;
    
    for(double n=0;n<nCycles;n+=1){
      if(n>0
         || (signZ > 0 && t1 > tShift)
         || (signZ < 0 && t1 < tShift)
         ){
        z1 = origin->GetZ() + slopeAbs*(t1 + signZ*n*2*TMath::Pi());
        
        points.push_back(make_shared<Point>(x1, y1, z1));
      }
      if(n>0
         || (signZ > 0 && t2 > tShift)
         || (signZ < 0 && t2 < tShift)
         ){
        z2 = origin->GetZ() + slopeAbs*(t2 + signZ*n*2*TMath::Pi());
        points.push_back(make_shared<Point>(x2, y2, z2));
      }
    }
    
  }
  return points;
}

void HelixProcessor::GetRandomPionHelix(const shared_ptr<Track> &track, Helix &pionHelix)
{
  unique_ptr<Point> pionVector = make_unique<Point>(RandSign()*RandDouble(config.minPx, config.maxPx),
                                                    RandSign()*RandDouble(config.minPy, config.maxPy),
                                                    RandSign()*RandDouble(config.minPz, config.maxPz));
  int pionCharge = RandSign();
  
  pionHelix = Helix(track->GetDecayPoint(), pionVector, pionCharge);
  vector<shared_ptr<Point>> pionPoints = GetPointsHittingSilicon(pionHelix);
  for(auto &p : pionPoints){p->SetIsPionHit(true);}
  pionHelix.SetPoints(pionPoints);
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

double HelixProcessor::GetHelixToPointDistance(const unique_ptr<Helix> &helix, const shared_ptr<Point> &point)
{
  int zSign = sgn(helix->momentum->GetZ());
  
  double t = TMath::Pi() - atan2(-(point->GetX() - helix->origin.GetX()),
                                  (point->GetY() - helix->origin.GetY()) );
  
//  if(helix->charge < 0) t = TMath::Pi()/2. - t;
  
  double x = helix->origin.GetX();
  double y = helix->origin.GetY();
  double z = helix->origin.GetZ() + helix->slopeAbs*t;
  
  if(helix->charge * zSign < 0){
    x += helix->radius*cos(t);
    y += helix->radius*sin(t);
  }
  else{
    x += helix->radius*sin(t);
    y += helix->radius*cos(t);
  }
  
  int nCycles = round(fabs(point->GetZ() - z) / (helix->slopeAbs * 2 * TMath::Pi()));
  z += zSign * nCycles * helix->slopeAbs * 2 * TMath::Pi();
  
  double stripModuleLength = 200;
  
  return pointsProcessor.distanceWithUncertainZ(make_shared<Point>(x,y,z), point, stripModuleLength);
}

double HelixProcessor::GetChi2toPoints(const unique_ptr<Helix> &helix, const vector<shared_ptr<Point>> &points)
{
  double chi2 = 0;
  
  for(auto &point : points){
    chi2 += pow(GetHelixToPointDistance(helix, point), 2);
  }
  
  return inf;
}
