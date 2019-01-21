//
//  HelixProcessor.cpp
//
//  Created by Jeremi Niedziela on 21/01/2019.
//

#include "HelixProcessor.hpp"


HelixProcessor::HelixProcessor(const shared_ptr<FitterConfig> &_config) :
config(_config),
pointsProcessor(make_unique<PointsProcessor>())
{

}

HelixProcessor::~HelixProcessor()
{
  
}

vector<int> HelixProcessor::AreIdentical(const unique_ptr<Helix> &h1, const unique_ptr<Helix> &h2)
{
  vector<int> reasons;
  shared_ptr<FitterConfig> config = h1->config;
  
  if(fabs(h1->GetOrigin()->GetX() - h2->GetOrigin()->GetX()) > config->GetToleranceX()) reasons.push_back(1);
  if(fabs(h1->GetOrigin()->GetY() - h2->GetOrigin()->GetY()) > config->GetToleranceY()) reasons.push_back(2);
  if(fabs(h1->GetOrigin()->GetZ() - h2->GetOrigin()->GetZ()) > config->GetToleranceZ()) reasons.push_back(3);
  if(fabs(h1->GetMomentum()->GetX() - h2->GetMomentum()->GetX()) > config->GetTolerancePx()) reasons.push_back(4);
  if(fabs(h1->GetMomentum()->GetY() - h2->GetMomentum()->GetY()) > config->GetTolerancePy()) reasons.push_back(5);
  if(fabs(h1->GetMomentum()->GetZ() - h2->GetMomentum()->GetZ()) > config->GetTolerancePz()) reasons.push_back(6);
  if(h1->GetCharge() != h2->GetCharge()) reasons.push_back(7);
  
  return reasons;
}

vector<Point> HelixProcessor::GetPointsHittingSilicon(const unique_ptr<Helix> &helix)
{
  vector<Point> points;
  unique_ptr<Point> origin = make_unique<Point>(*(helix->origin));
  
  double dh = sqrt(pow(origin->GetX(),2)+pow(origin->GetY(),2));
  double Rl, C, delta;
  double x1,y1,x2,y2,z1,z2,t1,t2;
  
  for(int iLayer=0;iLayer<nPixelLayers;iLayer++){
    Rl = layerR[iLayer];
    C = (Rl*Rl+dh*dh-helix->radius*helix->radius)/2.;
    
    delta = 4*origin->GetY()*origin->GetY()*C*C-4*dh*dh*(C*C-Rl*Rl*origin->GetX()*origin->GetX());
    if(delta < 0) continue;
    
    y1 = (2*origin->GetY()*C+sqrt(delta))/(2*dh*dh);
    y2 = (2*origin->GetY()*C-sqrt(delta))/(2*dh*dh);
    
    x1 = (C-y1*origin->GetY())/origin->GetX();
    x2 = (C-y2*origin->GetY())/origin->GetX();
    
    t1 = atan2((y1-origin->GetY()),(x1-origin->GetX()));
    t2 = atan2(y2-origin->GetY(),x2-origin->GetX());
    
    if(helix->charge < 0){
      t1 = TMath::Pi()/2. - t1;
      t2 = TMath::Pi()/2. - t2;
    }
    
    double nCycles = fabs(helix->GetNcycles());
    int signZ = sgn(helix->momentum->GetZ());
    double tShift = helix->tShift;
    double slopeAbs = helix->slopeAbs;
    
    for(double n=0;n<nCycles;n+=1){
      if(n>0
         || (signZ > 0 && t1 > tShift)
         || (signZ < 0 && t1 < tShift)
         ){
        z1 = origin->GetZ() + slopeAbs*(t1 + signZ*n*2*TMath::Pi());
        
        points.push_back(Point(x1, y1, z1));
      }
      if(n>0
         || (signZ > 0 && t2 > tShift)
         || (signZ < 0 && t2 < tShift)
         ){
        z2 = origin->GetZ() + slopeAbs*(t2 + signZ*n*2*TMath::Pi());
        points.push_back(Point(x2, y2, z2));
      }
    }
    
  }
  return points;
}

void HelixProcessor::CalculateNregularPoints(unique_ptr<Helix> &helix, int limit)
{
  vector<vector<Point>> pointsByLine = pointsProcessor->SplitPointsIntoLines(helix->points,
                                                                             config->GetLinesToleranceForRegularity());
  vector<double> possibleDistances;
  set<double> possibleDistancesSet;
  helix->nRegularPoints = 0;
  
  int nPointsForDistance;
  double zRegularityTolerance = config->GetZregularityTolerance();
  bool found;
  double testingDistance;
  
  for(auto line : pointsByLine){
    int iPoint;
    for(iPoint=0; iPoint < (int)line.size()-1; iPoint++){
      testingDistance = pointsProcessor->distance(line[iPoint], line[iPoint+1]);
      found = false;
      for(double dd : possibleDistances){
        if(fabs(testingDistance-dd) < zRegularityTolerance){found = true;break;}
      }
      if(found) continue;
      
      possibleDistances.push_back(testingDistance);
      nPointsForDistance = 0;
      
      for(auto line2 : pointsByLine){
        for(int i=0;i<(int)line2.size();i++){
          if(std::abs(pointsProcessor->distance(line2[0], line2[i])-i*testingDistance) < zRegularityTolerance)  nPointsForDistance++;
        }
      }
      if(nPointsForDistance > helix->nRegularPoints){
        helix->nRegularPoints = nPointsForDistance;
        if(helix->nRegularPoints > limit) return;
      }
    }
  }
}
