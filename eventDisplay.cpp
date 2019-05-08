#include "Helpers.hpp"
#include "Event.hpp"
#include "EventSet.hpp"
#include "Fitter.hpp"
#include "Display.hpp"
#include "ConfigManager.hpp"
#include "HelixProcessor.hpp"
#include "TrackProcessor.hpp"

uint searchRun = 1;
uint searchLumi = 1;
unsigned long long searchEvent = 2662;

string configPath = "configs/eventDisplay.md";
string cutLevel = "after_L0/";//after_L1/";

xtracks::EDataType dataType = xtracks::kSignal;
int setIter = kWino_M_300_cTau_10;
int iEvent = 2;

// 6  (q+, vz+, pz-) OK
// 10 (q+, vz+, pz-) OK - RECO
// 11 (q+, vz-, pz-) OK [strong breaking]
// 13 (q-, vz-, pz-) OK
// 14 (q-, vz+, pz+) OK
// 15 (q+, vz+, pz+) OK
// 18 (q-, vz+, pz+), (q+, vz+, pz+) OK [nice two helices]
// 19 (q-, vz-, pz-), (q+, vz-, pz-) OK [crazy stuff, probably 3rd soft particle there...]
// 20 (q+, vz+, pz+) OK [very high p_z]
// 21 (q-, vz+, pz-), (q+, vz+, pz+) OK [two big helices, one strongly breaking]
// 23 (q+, vz+, pz+) OK

bool injectPion = false;
bool fitHelix = true;

Display *display;
shared_ptr<EventSet> events;

map<string,any> filteredPointsOptions = {
  {"title", "Filtered Points"},
	{"binsMin" , 0},
	{"binsMax" , 100000},
  {"markerStyle", 20},
  {"markerSize", 1.0},
  {"color", kYellow}
};

void DrawHitsOrClusters(const shared_ptr<Event> event, int pointsType)
{
  
  vector<shared_ptr<Point>> hitsOrClusters;
  
  map<string,any> drawingOptions = {
    {"markerStyle", (pointsType==2) ? 20 : 22},
    {"markerSize", (pointsType==2) ? 1.0 : 2.0},
  };
  string typeName;
  
  if(pointsType == 0){
    if(!config.drawPionSimHits) return;
    
    hitsOrClusters = event->GetPionSimHits();
    drawingOptions["color"] = kCyan;
    typeName = "Pions hits ";
  }
  else if(pointsType == 1){
    if(!config.drawCharginoSimHits) return;
    
    hitsOrClusters = event->GetCharginoSimHits();
    drawingOptions["color"] = kMagenta;
    typeName = "Charginos hits ";
  }
  else if(pointsType == 2){
    if(!config.drawTrackerClusters) return;
    
    hitsOrClusters = event->GetTrackerClusters();
    drawingOptions["color"] = kYellow;
    typeName = "Tracker clusters ";
  }
  else if(pointsType == 3){
    if(!config.drawPionClusters) return;
    
    hitsOrClusters = event->GetPionClusters();
    drawingOptions["color"] = kBlue;
    typeName = "Pion clusters ";
  }
  
  map<string, vector<shared_ptr<Point>>> hitsOrClustersBySubDet;
  
  for(auto &[iter, name] : subDetMap){
    hitsOrClustersBySubDet[name] = vector<shared_ptr<Point>>();
  }
  
  for(auto &hit : hitsOrClusters){
    hitsOrClustersBySubDet[hit->GetSubDetName()].push_back(hit);
  }
  
  for(auto &[name, hitsVector] : hitsOrClustersBySubDet){
    if(hitsVector.size() == 0) continue;
    
    drawingOptions["title"] = (typeName+name).c_str();
    display->DrawSimplePoints(hitsVector, drawingOptions);
  }
}

shared_ptr<Event> GetEvent()
{
  EventSet events;
	//  events->LoadEventsFromFiles("/");
	events.LoadEventsFromFiles(cutLevel);
	
//  auto event = events.GetEvent(dataType, searchRun, searchLumi, searchEvent);
	auto event = events.At(dataType, setIter, iEvent);
	
	if(!event){
		cout<<"eventDisplay -- event not found"<<endl;
		exit(0);
	}
	
  JetCut jetCut;
	
	jetCut.SetPt(range<double>(30.0, inf));
	jetCut.SetChargedHadronEnergyFraction(range<double>(0.01,0.99));
	jetCut.SetNeutralHadronEnergyFraction(range<double>(0.01,0.99));
	
  eventProcessor.ApplyJetCut(event, jetCut);
	
  event->LoadAdditionalInfo();
  
	return event;
}

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  // create event display
  config = ConfigManager(configPath);
  display = new Display();
	
	auto event = GetEvent();
	shared_ptr<Track> track = event->GetTrack(0);
	
	const map<string,any> dedxOptions = {
		{"title", "dE/dx clusters"},
		{"binsMin" , 1},
		{"binsMax" , 10},
		{"nBins" , 5},
		{"markerStyle", 20},
		{"markerSize", 2.0}
	};
	
	display->DrawEvent(event, dedxOptions);
  event->Print();

  // -----------------------------------------------------------------------------------------------------
  // Helix fitting part
  // -----------------------------------------------------------------------------------------------------
  
  cout<<"Preparing hits, track and pion's helix"<<endl;
  
  DrawHitsOrClusters(event, 0); // pion sim hits
  DrawHitsOrClusters(event, 1); // chargino sim hits
  DrawHitsOrClusters(event, 2); // tracker clusters
  DrawHitsOrClusters(event, 3); // pion rec clusters
  
  
	const map<string,any> trueHelixOptions = {
		{"title", "True helix"},
		{"markerStyle", 20},
		{"markerSize", 0.2},
		{"color", kGreen}
	};
	
  auto allSimplePoints = event->GetTrackerClusters();
  
	if(injectPion){
		// Draw decay point to make sure that it's correctly located
    vector<shared_ptr<Point>> decayPoint;
    decayPoint.push_back(make_shared<Point>(track->GetDecayPoint().GetX(),
                                            track->GetDecayPoint().GetY(),
                                            track->GetDecayPoint().GetZ()));
		
		const map<string,any> decayPointOptions = {
			{"title", "Decay Point"},
			{"markerStyle", 22},
			{"markerSize", 2.0},
			{"color", kGreen}
		};
		
		display->DrawSimplePoints(decayPoint, decayPointOptions);
		
		// Draw true pion helix
    Helix pionHelix;
    helixProcessor.GetRandomPionHelix(track, pionHelix);
		display->DrawHelix(pionHelix,trueHelixOptions);
		cout<<"\n\nInjected pion helix:"<<endl;
		pionHelix.Print();
		
		// Calculate and draw points along the helix that hit the silicon
    auto pionPoints = helixProcessor.GetPointsHittingSilicon(pionHelix);
		for(auto &p : pionPoints){p->SetIsPionHit(true);}
		pionHelix.SetPoints(pionPoints);
		
		const map<string,any> pionPointsOptions = {
			{"title", "Pion points"},
			{"markerStyle", 21},
			{"markerSize", 1.6},
			{"color", kMagenta}
		};
		
		display->DrawSimplePoints(pionPoints, pionPointsOptions);
		
		// inject hits from pion into all points in the tracker
		allSimplePoints.insert(allSimplePoints.end(), pionPoints.begin(), pionPoints.end());
	}
	else{
    vector<Helix> truePionHelices;
    event->GetGenPionHelices(truePionHelices);
		
		for(auto &helix : truePionHelices){
			display->DrawHelix(helix,trueHelixOptions);
			helix.SetPoints(allSimplePoints);
			auto helixPoints = helix.GetPoints();
      filteredPointsOptions["title"] = "true helix points";
      filteredPointsOptions["color"] = kRed;
//      display->DrawSimplePoints(helixPoints, filteredPointsOptions);
			
			cout<<"\n\nTrue pion helix:"<<endl;
			helix.Print();
		}
	}
	
  for(auto p = allSimplePoints.begin(); p != allSimplePoints.end();){
    shared_ptr<Point> point = *p;
    if(point->GetSubDetName() == "TID" || point->GetSubDetName() == "TEC") p = allSimplePoints.erase(p);
    else p++;
  }
  
	if(fitHelix){
		cout<<"Fitting best helix"<<endl;
		auto fitter = make_unique<Fitter>();
		
    auto pionClusters = event->GetPionClusters();

    for(auto &p : allSimplePoints){
      int layer = -1;
      double minDist = inf;
      double pointR = sqrt(pow(p->GetX(), 2) + pow(p->GetY(), 2));
      
      for(int iLayer=0; iLayer<nLayers; iLayer++){
        double pointLayerDist = fabs(layerR[iLayer] - pointR);
        
        if(pointLayerDist < minDist){
          minDist = pointLayerDist;
          layer = iLayer;
        }
      }
      
      if(layer == track->GetNtrackerLayers()){
        Point trackPoint(layerR[layer] * cos(track->GetPhi())    + 10*event->GetVertex()->GetX(),
                         layerR[layer] * sin(track->GetPhi())    + 10*event->GetVertex()->GetY(),
                         layerR[layer] / tan(track->GetTheta())  + 10*event->GetVertex()->GetZ());

        if(pointsProcessor.distance(make_shared<Point>(trackPoint), p) < 100){
          pionClusters.push_back(p);
        }
      }
    }
    
    map<string,any> pionClustersOptions = {
      {"title", "Pion clusters"},
      {"binsMin" , 0},
      {"binsMax" , 100000},
      {"markerStyle", 20},
      {"markerSize", 2.0},
      {"color", kCyan}
    };
    
    display->DrawSimplePoints(pionClusters, pionClustersOptions);
    
    vector<int> rndIndices = { 2940,4839,3233,69,5098,1973,719,638,4618,2246,126,3029,98,3303,4436,2404,4013,384,4561,2588,3258,1339,1592,695,4015,318,2920,2655,2705,1463,5287,1830,1147,5261,2726,4184,436,854,5396,211,4325,4000,4872,4807,4204,1450,1034,735,2142,2498,5423,3070,1111,199,2843,1530,4205,5310,631,4892,474,411,2355,5057,3771,4725,5337,4038,3395,1624,3685,3659,2261,2270,4529,2739,3971,1858,3937,3497,4720,4575,1781,763,737,5419,2890,3097,1637,714,4655,3581,4442,4552,312,2205,410,601,5348,2104,2458,1404,5396,1459,744,1290,3196,38,5239,2848,2715,5173,4176,1086,1555,831,912,3518,2957,186,99,1220,3446,2735,5374,1836,4227,1006,2791,2660,4222,3686,1799,3741,1296,3466,1856,3639,4096,5390,5311,3037,1865,4700,2696,1145,1374,1155,2424,1614,843,4943,95,4081,1554,4391,4958,1674,1381,219,403,160,2476,3637,4077,497,3388,682,2879,4538,2592,3959,3376,2297,1280,5212,5304,3122,2664,3658,1364,2794,0,2684,5044,1734,4964,1812,770,2990,1622,2976,3722,2272,3261,1057,1416,3289,3997,1583,5025,5231,1172,2986,2538,2084,67,4056,4337,42,3655,3865,3003,1475,2524,3031,299,5220,4351,2237,775,1422,4931,3198,2929,4183,1777,844,425,253,2626,5425,873,3713,2089,4335,3638,2847,1453,5163,4673,183,5116,2784,2685,2072,5419,4474,4327,3116,1313,4211,1027,3441,4459,5202,5350,3885,4072,2916,2180,1647,4626,2628,1702,3266,3670,2692,215,5237,2478,4239,2473,3951,3513,4262,4687,1475,261,1969,4535,4590,2928,3812,4782,2278,2180,610,5406,2468,2277,2949,337,1531,3959,2400,713,4004,1219,4594,1636,5013,4508,2673,3005,1851,105,2636,1262,5180,3531,2779,3970,3393,1349,2443,5084,2150,5427,3550,3142,2141,2319,4541,1716,307,3384,1001,930,5149,2904,3928,83,2907,1344,1720,401,3683,3954,1765,2383,4254,3117,1453,1783,1570,3084,533,2361,2008,3540,1415,4235,4795,4871,999,1613,3664,5054,2645,2029,1055,4925,3603,496,1088,823,1099,4936,1246,5256,2324,1057,112,5144,291,629,1931,4497,1174,3978,2219,2980,4636,3739,3250,3138,4586,1181,2106,787,5197,1662,2759,3158,653,723,3617,4213,372,5065,3210,818,305,1385,4660,910,3498,380,1293,1654,1601,4155,1335,4124,570,1449,4514,2932,4400,4840,3593,1816,1646,1939,2940,1934,2755,2714,1311,4466,3698,1455,3044,2631,2858,4560,694,3023,2178,2355,2936,2693,2645,3821,595,2398,2795,4528,4223,2722,3910,1632,1209,5104,1836,2026,490,4671,1762,757,2763,1324,3340,1645,1653,4401,1044,4078,1676,3401,5216,231,2208,1671,3039,5433,1349,3013,3471,3446,2027,2500,1297,569,2618,1970,3518,3426,985,4876,5159,4641,5374,5181,5390,196,1096,530,1646,280,3780,903,5207,3196,5129,82,2391,129,2492,1386,2565,5294,3654,2298,310,391,2833,3707,2911,707,4646,446,1952,2320,3758,1179,1664,3736,746,5306,2803,2857,1340,443,2698,1029,117,886,5169,2455,4506,3632,1117,3857,3623,615,1274,2085,902,4407,1573,2245,879,5175,381,5256,4307,515,2237,2674,3182,2415,2978,3583,1848,1762,849,2804,65,2790,12,963,117,4788,334,1631,5095,2146,3754,4356,3402,1966,2924,4039,4153,3735,2053,2159,1468,3712,583,2628,186,4848,225,641,477,1672,3338,4403,2250,3421,744,5239,5288,3875,1711,5241,1387,3715,130,2045,2518,2312,5184,152,5109,3365,3208,5306,4002,4318,1030,3003,4385,1255,1993,4782,841,1873,5117,4515,3537,713,4152,1014,4447,3541,4107,2683,5026,5048,1176,774,4192,3721,3950,4093,1289,4940,1073,5389,2367,791,4032,4603,4120,4970,5313,5310,1748,165,1295,3037,3471,256,3854,2316,4735,1636,1746,95,980,1415,4841,583,4603,3243,861,458,3078,4351,3713,3591,4896,406,2163,3684,1885,1225,238,1925,4199,3790,2386,3492,462,3771,4287,5332,1353,1615,3984,1933,3579,2273,4410,3092,4185,3691,5352,5099,3011,4382,1047,1885,1850,2314,3798,3164,2701,4180,3884,2593,5233,2061,4426,1401,2635,1335,414,5375,3958,1162,1411,1741,4587,2405,3336,432,3334,455,2914,4823,3081,3471,4685,3948,115,56,1511,4661,40,2016,5060,2493,1251,790,160,1345,3033,4385,5289,2654,4959,4083,2746,4448,4,1225,3589,75,2240,2005,896,3513,2459,608,2946,2994,916,876,5416,2939,3023,3684,3304,2486,4309,2072,2285,2145,4861,1350,503,5345,4357,1255,456,3573,1449,3716,4745,1290,5114,2678,2824,564,4536,1252,2867,3681,1268,4007,2690,4593,4957,1473,1040,1266,775,3800,2272,5377,4416,1745,5158,575,4075,3767,215,476,1225,189,3920,3062,3512,3743,35,1363,5300,3990,1485,1902,5182,1157,831,2166,1599,629,2037,5038,4228,2168,3225,2654,646,3335,1544,3645,4761,4703,5134,3278,2079,3328,4422,2933,124,3103,26,3081,1907,632,1993,4182,851,4095,2110,2637,744,1802,5341,1708,1639,2376,4298,3112,611,3181,797,1220,815,3749,3167,142,4750,2438,4013,3135,2556,2831,5137,522,5211,1563,1975,4160,3246,4972,4637,156,3732,2769,4242,3844,472,3058,2557,491,3165,4861,3659,922,446,1931,1974,473,3428,3108,1450,4893,5275,806,4539,124,2874,1879,2657,4533,4139,2165,1076,3744,2639,3895,2952,1940,1548,3051,2187,699,5291,2456,1225,363,4699,4271,4786,5279,3784,3687,3420,5413,3187,251,2706,4114,1115,851,4972,96,1347,1699,3241,4883,928,3691,4740,2,3690,1462,4734,1266,1574,1642,1112,4898,1672,2500,1942,3085,3494,440,2899,4634,3164,3727,3261, };
    
    // Turn this on to inject some noise
    cout<<endl;
    for(int i=0;i<1000;i++){
      int r = rndIndices[i];
//      int r = RandInt(0, (int)allSimplePoints.size()-1);
      cout<<r<<",";
      pionClusters.insert(pionClusters.end(),allSimplePoints[r]);
    }
    cout<<endl;
    
    display->DrawSimplePoints(pionClusters, filteredPointsOptions);
    
    vector<Helix> fittedHelices = fitter->FitHelices(pionClusters, *track, *event->GetVertex());
    		
    map<string,any> bestHelixOptions = {
      {"title", "Best helix"},
      {"markerStyle", 20},
      {"markerSize", 0.2},
      {"color", kRed}
    };
    
    map<string,any> helixVertexOptions = {
      {"title", "Helix vertex"},
      {"binsMin" , 0},
      {"binsMax" , 100000},
      {"markerStyle", 20},
      {"markerSize", 2.0},
      {"color", kYellow}
    };
    
    for(int iHelix=0; iHelix<fittedHelices.size(); iHelix++){

      bestHelixOptions["title"] = ("Helix "+to_string(fittedHelices[iHelix].uniqueID)).c_str();
      display->DrawShrinkingHelix(fittedHelices[iHelix], bestHelixOptions);
      fittedHelices[iHelix].Print();
      
      helixVertexOptions["markerStyle"] = 20;
      vector<shared_ptr<Point>> helixVertex = fittedHelices[iHelix].GetPoints();
      display->DrawSimplePoints(helixVertex, helixVertexOptions);
    }
  }
  
     
  gEve->Redraw3D(true);
  theApp.Run();
  return 0;
}
