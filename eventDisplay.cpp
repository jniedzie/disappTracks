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
string cutLevel = "after_L1/4layers/";//after_L1/";

xtracks::EDataType dataType = xtracks::kSignal;
int setIter = kWino_M_300_cTau_10;
int iEvent = 31;

// 13 - should be possible, two helices
// 14 - nice helix, but probably missing first strip hit
// 26 - beautiful example, must be possible, but track misreconstructed
// 27 - two helices, one looks fine
// 31 - two helices, one has many hits and looks good (OK)

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
//  events.LoadEventsFromFiles(cutLevel);
  events.LoadEventFromFiles(dataType, setIter, iEvent, cutLevel);
	
//  auto event = events.GetEvent(dataType, searchRun, searchLumi, searchEvent);
	auto event = events.At(dataType, setIter, 0);
	
	if(!event){
		cout<<"eventDisplay -- event not found"<<endl;
		exit(0);
	}
	
  JetCut jetCut;
	
	jetCut.SetPt(range<double>(30.0, inf));
	jetCut.SetChargedHadronEnergyFraction(range<double>(0.01,0.99));
	jetCut.SetNeutralHadronEnergyFraction(range<double>(0.01,0.99));
	
  eventProcessor.ApplyJetCut(event, jetCut);
	
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
    vector<Helix> truePionHelices = event->GetGenPionHelices();
		
		for(auto &helix : truePionHelices){
			display->DrawHelix(helix,trueHelixOptions);
//      helix.SetPoints(allSimplePoints);
//      auto helixPoints = helix.GetPoints();
//      filteredPointsOptions["title"] = "true helix points";
//      filteredPointsOptions["color"] = kRed;
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
/*
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
  */
    map<string,any> pionClustersOptions = {
      {"title", "All clusters"},
      {"binsMin" , 0},
      {"binsMax" , 100000},
      {"markerStyle", 20},
      {"markerSize", 2.0},
      {"color", kCyan}
    };

    
    vector<int> rndIndices = { 3271,817,358,2770,1294,1231,3937,142,2281,3131,2754,2122,3986,2757,1196,2668,4029,4059,3577,721,3372,1175,195,4454,3275,4331,3537,3874,4337,2242,4000,4255,3457,1815,756,3945,1811,99,3342,3262,1656,3830,4226,2327,4022,1383,2522,1732,666,4096,2755,1687,173,836,1120,2726,1501,2287,3546,1022,246,1277,2140,496,2439,2084,1239,4281,3106,2178,954,2741,2503,3305,2184,173,1526,842,4355,2473,2065,81,754,1395,3177,4101,1984,4451,3208,1849,3062,1768,3635,8,1380,2587,4570,3470,2771,2545,702,1562,4124,767,2667,880,4376,3570,1472,3745,2311,1177,4186,3065,3551,3389,3497,2245,1135,1881,4128,2852,4479,3376,4327,4037,381,2530,426,980,3137,1325,2753,843,2897,226,1087,318,736,2766,3391,3392,1506,4243,4066,2960,104,3787,778,1709,638,2322,2447,3998,2567,137,4454,1109,2282,4427,2127,2183,3029,4080,2555,2632,3123,3697,2739,4141,2520,1128,378,1300,913,3254,2821,1529,965,2232,2503,1747,3481,1515,4007,2563,4490,1353,2540,49,1872,1616,403,1319,2302,2319,1731,152,1248,1297,1912,201,1183,3038,1195,2866,761,1926,496,27,3271,3122,490,4407,585,926,205,1729,4518,3763,1509,2466,3152,1375,3848,1651,624,938,186,3241,3416,826,2272,2707,552,2311,2835,621,4348,2757,3846,3765,1887,4115,4155,854,604,1405,345,2591,1073,4335,2301,875,3478,1776,1460,1217,3984,2458,2052,3625,4155,1236,3188,1306,3498,3347,314,2819,3745,1228,4145,4033,1117,3589,4097,2527,113,132,2811,3655,317,1710,3805,1383,4458,2635,3344,2378,3804,773,1214,4386,2487,2264,3759,1192,220,3439,2338,2860,3916,1727,1402,3619,1371,666,800,699,2659,3824,2713,4418,442,3668,236,4384,467,668,3667,2334,2697,1586,1972,2468,3745,4299,3373,3737,1337,3009,2852,3153,704,251,1525,1171,3769,52,74,1533,2139,1169,4116,2750,3365,2182,3497,2903,3471,2929,2455,2840,1383,1843,1749,1880,1461,681,3816,2110,2145,3724,1576,42,3379,62,1667,2194,1675,3289,4098,3316,110,3570,2073,3591,1968,175,4288,1846,792,1588,2852,1218,2134,576,4351,2148,862,1057,4240,1007,3885,2605,1257,3121,378,2530,238,525,3364,2478,2026,3533,2109,1345,981,883,1102,2704,3885,457,2020,4039,3065,920,2546,2922,3195,320,3001,4547,3540,1786,1728,1553,1434,856,3079,1737,2456,452,1318,1599,3089,1805,3532,1459,3003,2682,3171,526,19,1795,2701,2105,3652,1961,217,2220,1283,98,838,2297,3264,1254,13,416,2324,1131,3710,705,1237,4274,824,3642,3773,2540,2551,3,2882,3873,4416,2131,1612,2221,779,3505,3458,1538,1584,2463,1539,1567,20,4060,4084,3748,3577,1866,438,3655,4474,2348,3517,294,656,3300,208,1239,3353,1129,3897,538,1285,3845,2838,555,2220,336,264,4064,2078,3767,2241,1785,210,194,56,4026,2171,1149,1350,3654,221,4102,4080,4177,2846,2841,2690,4571,3833,738,265,674,3139,354,2917,4108,3783,4025,2842,1015,563,3691,1024,4102,38,1989,2053,1257,3345,887,396,2813,4161,1055,2610,1589,1716,1334,833,2592,3322,1448,3322,2176,3583,2540,1057,1986,2233,4053,3398,4238,2803,4394,3543,3511,1768,3745,3863,235,2370,1493,4077,6,411,3904,1525,194,1774,3255,4328,598,2641,4482,3804,1937,3794,741,242,2845,2373,984,4077,4346,3453,2446,243,1102,1628,1238,4186,2122,3477,758,2549,1605,2110,512,988,3263,1448,4193,4215,1905,4444,4530,4254,599,1572,3219,4296,2753,4027,4539,4306,739,4364,2655,1627,1850,3041,3397,3969,2197,1391,1921,3413,4507,979,2867,3654,2715,95,2633,1293,2004,1962,3575,3867,2056,1236,2720,1745,3682,1993,693,3481,3020,366,1014,3612,1578,3920,3874,500,3373,1823,3065,2692,1185,2321,166,4404,3985,1743,4204,356,4142,1491,4534,4530,2676,2540,2831,1161,1090,55,1551,4563,314,2696,2244,3409,1743,1033,308,1859,1339,2864,2783,2031,4048,3166,3816,4118,3720,1642,837,1121,132,2566,2398,3213,2743,196,4216,2714,1892,2070,2236,3954,3383,2541,3150,1241,2774,1538,4080,3274,530,2113,3882,4005,1667,407,22,2407,1396,1195,1688,2925,4230,3928,257,334,1232,386,901,3174,1954,879,3100,1981,3738,3397,928,4550,1112,725,4261,4265,88,1100,1629,485,105,4256,1274,3126,2463,428,4301,3028,2039,942,4527,145,2566,13,757,304,2304,1827,4334,273,1755,4372,3909,2112,1654,4255,676,1280,2376,2627,3322,2349,1936,1589,2037,2313,653,2142,1157,1637,2347,364,1900,2710,2461,1784,128,949,4482,329,1086,2724,1208,2934,2669,3412,1523,351,3140,2349,4512,2616,46,1158,3581,2185,3134,4181,1424,1017,1418,171,1381,470,1980,563,872,4537,3959,628,1645,155,2609,2550,73,1496,716,2205,413,1489,798,4074,3909,3640,4188,2576,4328,2926,3593,1376,4524,235,1646,1334,2821,4390,4366,1250,2464,4550,2733,1722,4019,1788,167,3618,337,1854,189,1147,3364,3961,803,34,1590,3750,1541,2956,2232,1922,554,264,859,3050,1372,3964,3160,2850,2892,2425,3364,2330,175,1491,4027,3205,3503,4040,1977,720,681,292,553,460,1076,1393,3421,4075,2368,1590,751,3636,3320,2929,1860,487,3847,4162,2386,1702,662,4512,4473,2573,1294,143,1412,78,3484,1970,4103,2272,2225,2606,59,2378,3780,40,28,4075,272,4429,4088,1870,3409,2354,1308,335,2408,2540,485,1243,2241,1399,234,918,414,890,4364,2140,519,1252,317,1894,2601,2345,828,1875,1099,2207,4467,2402,3486,1175,4547,3199,4117,417,789,1734,4305,2197,1017,3387,4373,2422,4551,971,476,901,4364,2675,419,2879,343,395,753,778,1453,2381,625,2761,3907,3110,246,30,2321,85,160,230,3485,2188,4220,3090,439,4277,1858,1881,3083,4359,3711,1368,3101,4064,304,224,3761,266,4177,3781,2254,407,1279,1110,4017,2817,4369,3554,965,3878,1668,4287,2089,4091,2131,1296,32,802,4050,1928,4237,4309,746,683,20,3113,2018,2145,2506,3274,711,3872,3769,1913,1356,1683,3938,2193,4424,1792,2240,3124,651,2642,2764,678,2891,735,44,4035,3499,3928,3586,2718,4039,625,2155,4364,4163,471,294,641,1314,1217,1310,4112,2404,1573,2223,1430,1659,2925,3177,2427,1730,4440,961,3788,76,2689,386,3527,3241,3947,81,4021,3804,2705,4065,1549,1730,2725,153,2563,1316,2102,2922,3649,779,4337,2191,2345,4133,4431,2168,2815,2263,1198,1231,1146,2462,2551,4094,1906,4352,1596,990,2413,1497,1955,256,659,4364,2985,2157,3122,4126,4189,3608,3066,2945,2680,2041,2281,179,2295,2366,417,3524,2352,3276,900,3428,3392,3421,1790,2455,1980,3033,52,4208,1553,2399,2278,3830,1173,933,4527,4515,4451,4303,882,259,4297,595,4239,2043,320,2563,4560,2488,283,2563,1442,2604,2272,3243,1083,1337,1779,1759,4491,2446,218,1557,1236,2139,4457,1222,3717,2178,1143,4435,4442,1206,2290,3581,183,4073,1134,3635,2873,7,1083,262,110,20,4137,3372,1543,2416,854,2642,3627,658,405,3933,1658,1663,1821,346,4249,971,1222,1772,456,4394,1929,913,3577,3647,3279,903,4457,1855,1536,4327,1901,4377,3359,3979,1856,3865,4250,875,3867,2861,3118,1323,3009,4556,1185,4142,3242,758,3488,497,1768,530,2538,463,728,3911,154,2995,1404,3489,21,1061,3222,1928,776,235,4528,4122,438,1138,3154,4325,496,2229,99,1254,4340,3688,4079,120,2143,907,1083,4460,3663,3546,4304,4264,4146,14,250,384,338,403,856,2690,3543,2623,3749,4169,4342,84,285,3405,2740,3146,1893,523,498,2539,1707,778,3361,3521,2070,2149,666,2609,2717,1244,2858,2304,1151,2715,4307,3927,2659,1704,4253,2944,272,2758,349,1384,604,3489,3584,4461,1596,228,4273,3170,1879,809,1628,1642,2403,1425,2980,4438,1876,2932,4472,1581,559,3544,993,4135,4559,1184,106,831,1158,1186,2866,1624,1559,876,2599,1461,3680,837,584,1790,1722,661,2742,2198,1138,1823,4159,1428,3701,4307,4257,3587,1483,3808,4053,1572,4132,1179,2869,3647,374,4135,2524,3962,2321,4480,3857,426,2993,46,2778,614,2139,793,1003,3478,910,2248,4214,176,1693,1817,66,3773,1927,1943,896,2357,1579,2222,3067,3019,2548,1514,4437,281,3045,2483,3419,2217,4550,3486,2419,2678,631,285,639,1424,406,2223,4009,4476,4102,3089,1188,3134,549,1815,1330,0,3563,2574,830,697,2248,3564,507,1649,2332,3253,1901,117,512,1598,415,757,892,1222,3501,1453,130,1511,1456,4443,1601,245,1049,3059,136,1804,2077,496,2629,2648,1432,3026,1542,1117,1147,3499,36,2814,3618,3364,4181,4274,462,4103,2792,4268,87,2614,2894,430,3863,1113,359,141,2117,4416,3802,3794,3004,2965,4285,151,4029,556,2869,3033,1215,4311,2721,2045,4101,3964,655,661,2945,2176,1687,901,3184,4049,2451,2839,3391,1816,3396,3049,1923,2916,812,1930,801,976,439,3202,1117,2446,4431,735,4078,3498,4241,43,3513,615,230,4290,3555,2585,3375,4409,1644,4427,4131,1815,790,332,2989,2907,4536,597,2836,1337,1224,1648,3084,487,383,3268,1506,2544,991,2836,2838,1148,1544,4240,3125,3228,4307,4056,1052,2762,1303,4247,3254,715,2034,2525,2240,4556,4229,1210,1843,3209,1525,2191,2514,51,3927,911,1872,1775,4015,200,2477,832,711,4111,3148,2857,380,3872,1715,2748,2229,546,727,1914,2058,1935,4109,401,3658,1777,3681,4357,3082,2897,552,2726,1560,2197,3908,473,2656,103,1780,1855,1773,3322,2391,1011,2989,1737,3846,4251,3986,3924,2878,742,3075,717,2625,2797,1226,3968,4210,3024,2378,2487,1884,2769,822,1691,4420,3653,4341,2455,2885,2350,2027,2045,2592,1519,1083,2235,3008,2028,2255,3345,3208,739,738,3054,2027,162,1323,3887,401,2791,2448,373,2611,1157,2497,1190,4110,3971,3248,4318,4149,376,2948,977,292,1098,4454,3624,2291,4265,1057,2240,4526,1720,4053,1116,2919,3971,1186,2168,1909,3403,2428,1589,1184,3188,2545,4200,1574,561,1582,1718,1937,4534,4191,3014,2112,1304,4035,1175,1853,606,3212,3006,265,2647,3849,4180,3829,3416,1278,1717,4176,2190,78,3540,4141,3384,36,3403,1224,1201,4111,1607,11,4502,3420,74,4002,4057,1103,3171,1585,5,4172,108,2361,4126,400,576,2251,4327,2777,1952,1022,3538,3093,2044,599,3880,1068,3223,1108,3218,3500,386,2476,805,2738,1474,748,3544,2154,4006,2261,383,4181,635,4461,4444,151,3876,437,2535,3874,2553,4306,1406,317,4516,1885,805,3743,4493,1243,2104,3867,3280,814,3900,3669,3532,135,2505,3307,1786,4224,3785,110,3061,2603,3519,700,3466,3445,1463,1197,2077,10,4240,895,1102,2811,176,2951,1115,258,308,3592,467,1447,4483,671,3675,3725,3329,1070,1230,2903,3215,223,3802,3704,1331,1086,3204,1191,1806,2962,680,1168,1324,512,2621,827,879,689,966,1323,3990,3748,2865,4505,3318,2588,1960,1868,3103,1681,3239,1145,1588,4403,2428,3243,1350,1892,3409,967,3124,4260,452,3517,1122 };
    
    // Turn this on to inject some noise
    cout<<endl;
    for(int i=0; i<0; i++){
      int r = rndIndices[i];
//      int r = RandInt(0, (int)allSimplePoints.size()-1);
//      cout<<r<<",";
      pionClusters.insert(pionClusters.end(),allSimplePoints[r]);
    }
    cout<<endl;
    
    
    
    display->DrawSimplePoints(pionClusters, pionClustersOptions);
    
    auto start = now();
//    vector<Helix> fittedHelices = fitter->FitHelices(allSimplePoints, *track, *event->GetVertex());
    vector<Helix> fittedHelices = fitter->FitHelices(pionClusters, *track, *event->GetVertex());
    auto end = now();
    
    cout<<"Fitting time: "<<duration(start, end)<<endl;
    
    
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
