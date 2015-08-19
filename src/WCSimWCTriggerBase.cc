#include "WCSimWCTriggerBase.hh"
#include "WCSimWCPMT.hh"
#include "WCSimWCDigi.hh"
#include "WCSimWCHit.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4DigiManager.hh"
#include "G4ios.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#include "WCSimDetectorConstruction.hh"
#include "WCSimPmtInfo.hh"
#include "WCSimDarkRateMessenger.hh"

#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include <vector>
// for memset
#include <cstring>
#include <iostream>

#ifndef WCSIMWCTRIGGERBASE_VERBOSE
//#define WCSIMWCTRIGGERBASE_VERBOSE
#endif

const double WCSimWCTriggerBase::offset = 950.0 ; // ns. apply offset to the digit time
const double WCSimWCTriggerBase::eventgateup = 950.0 ; // ns. save eventgateup ns after the trigger time
const double WCSimWCTriggerBase::eventgatedown = -400.0 ; // ns. save eventgateup ns before the trigger time
const double WCSimWCTriggerBase::LongTime = 100000.0 ; // ns = 0.1ms. event time


WCSimWCTriggerBase::WCSimWCTriggerBase(G4String name,
				       WCSimDetectorConstruction* myDetector,
				       WCSimWCDAQMessenger* myMessenger)
  :G4VDigitizerModule(name)
{
  G4String colName = "WCDigitizedCollection";
  this->myDetector = myDetector;
  collectionName.push_back(colName);
  DigiHitMap.clear();
  TriggerTimes.clear();
  TriggerTypes.clear(); 
  TriggerInfos.clear(); 

  if(myMessenger) {
    DAQMessenger = myMessenger;
    DAQMessenger->TellMeAboutTheTrigger(this);
    DAQMessenger->TellTrigger();
  }

  digitizeCalled = false;
}
  
WCSimWCTriggerBase::~WCSimWCTriggerBase(){
}

void WCSimWCTriggerBase::AdjustNHitsThresholdForNoise()
{
  int npmts = this->myDetector->GetTotalNumPmts();
  double trigger_window_seconds = nhitsWindow * 1E-9;
  double dark_rate_Hz = PMTDarkRate * 1000;
  double average_occupancy = dark_rate_Hz * trigger_window_seconds * npmts;
  
  G4cout << "Average number of PMTs in detector active in a " << nhitsWindow
	 << "ns window with a dark noise rate of " << PMTDarkRate
	 << "kHz is " << average_occupancy
	 << " (" << npmts << " total PMTs)"
	 << G4endl
	 << "Updating the NHits threshold, from " << nhitsThreshold
	 << " to " << nhitsThreshold + round(average_occupancy) << G4endl;
  nhitsThreshold += round(average_occupancy);
}

void WCSimWCTriggerBase::Digitize()
{
  if(nhitsAdjustForNoise && !digitizeCalled) {
    AdjustNHitsThresholdForNoise();
    digitizeCalled = true;
  }

  //Input is collection of all digitized hits that passed the threshold
  //Output is all digitized hits which pass the trigger
  
  DigiHitMap.clear();
  TriggerTimes.clear();
  TriggerTypes.clear(); 
  TriggerInfos.clear(); 

  //This is the output digit collection
  DigitsCollection = new WCSimWCDigitsCollection ("/WCSim/glassFaceWCPMT",collectionName[0]);

  G4DigiManager* DigiMan = G4DigiManager::GetDMpointer();

  // Get the Digitized hits collection ID
  G4int WCDCID = DigiMan->GetDigiCollectionID("WCDigitizedStoreCollection");
  // Get the PMT Digits Collection
  WCSimWCDigitsCollection* WCDCPMT = 
    (WCSimWCDigitsCollection*)(DigiMan->GetDigiCollection(WCDCID));

  // Do the work  
  if (WCDCPMT) {
    DoTheWork(WCDCPMT);
  }
  
  StoreDigiCollection(DigitsCollection);
}

void WCSimWCTriggerBase::AlgNHits(WCSimWCDigitsCollection* WCDCPMT, bool remove_hits, bool test)
{

  //if test is true, we run the algorithm with 1/2 the threshold, and kTriggerNHitsTest
  //for testing multiple trigger algorithms at once
  int this_nhitsThreshold = nhitsThreshold;
  TriggerType_t this_triggerType = kTriggerNHits;
  if(test) {
    this_nhitsThreshold /= 2;
    this_triggerType = kTriggerNHitsTest;
  }

  //Now we will try to find triggers
  //loop over PMTs, and Digits in each PMT.  If nhits > Threshhold in a time window, then we have a trigger

  int ntrig = 0;
  int window_start_time = 0;
  int window_end_time   = WCSimWCTriggerBase::LongTime - nhitsWindow;
  int window_step_size  = 5; //step the search window along this amount if no trigger is found
  float lasthit;
  std::vector<int> digit_times;
  bool first_loop = true;

  G4cout << "WCSimWCTriggerBase::AlgNHits. Number of entries in input digit collection: " << WCDCPMT->entries() << G4endl;
#ifdef WCSIMWCTRIGGERBASE_VERBOSE
  int temp_total_pe = 0;
  for (G4int i = 0 ; i < WCDCPMT->entries() ; i++) {
    temp_total_pe += (*WCDCPMT)[i]->GetTotalPe();
  }
  G4cout << "WCSimWCTriggerBase::AlgNHits. " << temp_total_pe << " total p.e. input" << G4endl;
#endif

  // the upper time limit is set to the final possible full trigger window
  while(window_start_time <= window_end_time) {
    int n_digits = 0;
    float triggertime; //save each digit time, because the trigger time is the time of the first hit above threshold
    bool triggerfound = false;
    digit_times.clear();
    
    //Loop over each PMT
    for (G4int i = 0 ; i < WCDCPMT->entries() ; i++) {
      //int tube=(*WCDCPMT)[i]->GetTubeID();
      //Loop over each Digit in this PMT
      for ( G4int ip = 0 ; ip < (*WCDCPMT)[i]->GetTotalPe() ; ip++) {
	int digit_time = (*WCDCPMT)[i]->GetTime(ip);
	//hit in trigger window?
	if(digit_time >= window_start_time && digit_time <= (window_start_time + nhitsWindow)) {
	  n_digits++;
	  digit_times.push_back(digit_time);
	}
	//G4cout << digit_time << G4endl;
	//get the time of the last hit (to make the loop shorter)
	if(first_loop && (digit_time > lasthit))
	  lasthit = digit_time;
      }//loop over Digits
    }//loop over PMTs

    //if over threshold, issue trigger
    if(n_digits > this_nhitsThreshold) {
      ntrig++;
      //The trigger time is the time of the first hit above threshold
      std::sort(digit_times.begin(), digit_times.end());
      triggertime = digit_times[this_nhitsThreshold];
      triggertime -= (int)triggertime % 5;
      TriggerTimes.push_back(triggertime);
      TriggerTypes.push_back(this_triggerType);
      TriggerInfos.push_back(std::vector<Float_t>(1, n_digits));
      triggerfound = true;
    }

#ifdef WCSIMWCTRIGGERBASE_VERBOSE
    if(n_digits)
      G4cout << n_digits << " digits found in 200nsec trigger window ["
	     << window_start_time << ", " << window_start_time + nhitsWindow
	     << "]. Threshold is: " << this_nhitsThreshold << G4endl;
#endif

    //move onto the next go through the timing loop
    if(triggerfound) {
      window_start_time = triggertime + WCSimWCTriggerBase::eventgateup;
    }//triggerfound
    else {
      window_start_time += window_step_size;
    }

    //shorten the loop using the time of the last hit
    if(first_loop) {
#ifdef WCSIMWCTRIGGERBASE_VERBOSE
      G4cout << "Last hit found to be at " << lasthit
	     << ". Changing window_end_time from " << window_end_time
	     << " to " << lasthit - (nhitsWindow - 10)
	     << G4endl;
#endif
      window_end_time = lasthit - (nhitsWindow - 10);
      first_loop = false;
    }
  }
  
  //call FillDigitsCollection() if at least one trigger was issued
  G4cout << "Found " << ntrig << " NHit triggers" << G4endl;
  FillDigitsCollection(WCDCPMT, remove_hits, this_triggerType);
}

void WCSimWCTriggerBase::AlgNHitsThenITC(WCSimWCDigitsCollection* WCDCPMT, bool remove_hits)
{

  //Now we will try to find triggers
  //loop over PMTs, and Digits in each PMT. 
  // If nhits > Threshhold in a time window, then we have a trigger
  // If this cut fails, attempt an ITC ratio trigger

  int ntrig = 0;
  int window_start_time = 0;
  int window_end_time   = WCSimWCTriggerBase::LongTime - nhitsWindow;
  int window_step_size  = 5; //step the search window along this amount if no trigger is found
  float lasthit;
  std::vector<int> digit_times;
  std::vector<int> digit_times_itc_small, digit_times_itc_large;
  bool first_loop = true;

  G4cout << "WCSimWCTriggerBase::AlgNHits. Number of entries in input digit collection: " << WCDCPMT->entries() << G4endl;
#ifdef WCSIMWCTRIGGERBASE_VERBOSE
  int temp_total_pe = 0;
  for (G4int i = 0 ; i < WCDCPMT->entries() ; i++) {
    temp_total_pe += (*WCDCPMT)[i]->GetTotalPe();
  }
  G4cout << "WCSimWCTriggerBase::AlgNHits. " << temp_total_pe << " total p.e. input" << G4endl;
#endif

  // the upper time limit is set to the final possible full trigger window
  while(window_start_time <= window_end_time) {
    int n_digits = 0;
    int n_digits_itc_small = 0, n_digits_itc_large = 0;
    float triggertime; //save each digit time, because the trigger time is the time of the first hit above threshold
    bool triggerfound = false;
    digit_times.clear();
    
    //Loop over each PMT & count NDigits in window [window_start_time, window_start_time + nhitsWindow]
    //Also count in two extra windows for the ITC cut
    // [window_start_time, window_start_time + itcSmallWindow]  
    // [window_start_time - itcLargeWindowLow, window_start_time + itcLargeWindowHigh]  
    for (G4int i = 0 ; i < WCDCPMT->entries() ; i++) {
      //int tube=(*WCDCPMT)[i]->GetTubeID();
      //Loop over each Digit in this PMT
      for ( G4int ip = 0 ; ip < (*WCDCPMT)[i]->GetTotalPe() ; ip++) {
	int digit_time = (*WCDCPMT)[i]->GetTime(ip);
	//hit in trigger window?
	if(digit_time >= window_start_time && digit_time <= (window_start_time + nhitsWindow)) {
	  n_digits++;
	  digit_times.push_back(digit_time);
	}
        //hit in the small ITC window?
        if((digit_time >= window_start_time) && (digit_time <= (window_start_time + itcSmallWindow))) {
          n_digits_itc_small++;
          digit_times_itc_small.push_back(digit_time);
        }
        //hit in the large ITC window?
        if((digit_time >= (window_start_time - itcLargeWindowLow)) && digit_time <= (window_start_time + itcLargeWindowHigh)) {
          n_digits_itc_large++;
          digit_times_itc_large.push_back(digit_time);
        }
	//G4cout << digit_time << G4endl;
	//get the time of the last hit (to make the loop shorter)
	if(first_loop && (digit_time > lasthit))
	  lasthit = digit_time;
      }//loop over Digits
    }//loop over PMTs

    //if over threshold, issue trigger
    if(n_digits > nhitsThreshold) {
      ntrig++;
      //The trigger time is the time of the first hit above threshold
      std::sort(digit_times.begin(), digit_times.end());
      triggertime = digit_times[nhitsThreshold];
      triggertime -= (int)triggertime % 5;
      //save the trigger information
      TriggerTimes.push_back(triggertime);
      TriggerTypes.push_back(kTriggerNHits);
      TriggerInfos.push_back(std::vector<Float_t>(1, n_digits));
      triggerfound = true;
#ifdef WCSIMWCTRIGGERBASE_VERBOSE
      G4cout << EnumAsString(kTriggerNHits) << " trigger passed with time " << triggertime << G4endl;
#endif
    }//NDigits trigger passed

    //The simple NHits trigger hasn't been passed. See if the ITC ratio trigger can be passed
    double itc_ratio = (double)n_digits_itc_small / (double)n_digits_itc_large;
    if(!triggerfound && itc_ratio > itcRatioThreshold) {
      ntrig++;
      //The trigger time is the time of the first hit TEMPORARY
      std::sort(digit_times.begin(), digit_times.end());
      triggertime = digit_times[0];
      triggertime -= (int)triggertime % 5;
      std::vector<Float_t> triggerinfo;
      triggerinfo.push_back(itc_ratio);
      triggerinfo.push_back(n_digits_itc_small);
      triggerinfo.push_back(n_digits_itc_large);
      //save the trigger information
      TriggerTimes.push_back(triggertime);
      TriggerTypes.push_back(kTriggerITCRatio);
      TriggerInfos.push_back(triggerinfo);
      triggerfound = true;
#ifdef WCSIMWCTRIGGERBASE_VERBOSE
      G4cout << EnumAsString(kTriggerITCRatio) << " trigger passed with time " << triggertime << G4endl;
#endif
    }//ITC trigger passed
    

#ifdef WCSIMWCTRIGGERBASE_VERBOSE
    if(n_digits)
      G4cout << n_digits << " digits found in 200nsec trigger window ["
	     << window_start_time << ", " << window_start_time + nhitsWindow
	     << "]. Threshold is: " << nhitsThreshold << G4endl;
#endif

    //move onto the next go through the timing loop
    if(triggerfound) {
      window_start_time = triggertime + WCSimWCTriggerBase::eventgateup;
    }//triggerfound
    else {
      window_start_time += window_step_size;
    }

    //shorten the loop using the time of the last hit
    if(first_loop) {
#ifdef WCSIMWCTRIGGERBASE_VERBOSE
      G4cout << "Last hit found to be at " << lasthit
	     << ". Changing window_end_time from " << window_end_time
	     << " to " << lasthit - (nhitsWindow - 10)
	     << G4endl;
#endif
      window_end_time = lasthit - (nhitsWindow - 10);
      first_loop = false;
    }
  }
  
  //call FillDigitsCollection() if at least one trigger was issued
  G4cout << "Found " << ntrig << " NHitThenITC triggers" << G4endl;
  FillDigitsCollection(WCDCPMT, remove_hits, kTriggerUndefined);
}

void WCSimWCTriggerBase::FillDigitsCollection(WCSimWCDigitsCollection* WCDCPMT, bool remove_hits, TriggerType_t save_triggerType)
{

  G4String WCIDCollectionName = myDetector->GetIDCollectionName();
  G4float timingConstant = 0.0;
  WCSimPMTObject * PMT = myDetector->GetPMTPointer(WCIDCollectionName); //for hit time smearing


  //saveFailuresMode = 0 - save only triggered events
  //saveFailuresMode = 1 - save both triggered & not triggered events
  //saveFailuresMode = 2 - save only not triggered events
  if(TriggerTimes.size()) {
    if(saveFailuresMode == 2)
      return;
  }
  else {
    if(saveFailuresMode == 0)
      return;
    TriggerTypes.push_back(kTriggerFailure);
    TriggerTimes.push_back(saveFailuresTime);
    TriggerInfos.push_back(std::vector<Float_t>(1, -1));
    save_triggerType = kTriggerFailure;
  }

  //  WCSimPMTObject * PMT = myDetector->GetPMTPointer("glassFaceWCPMT"); //for hit time smearing
  //>>>>>>> 588f6a1a57690c567128cc870173732325d6317d

  //Loop over trigger times
  for(unsigned int itrigger = 0; itrigger < TriggerTimes.size(); itrigger++) {
    TriggerType_t triggertype = TriggerTypes[itrigger];
    //check if we've already saved this trigger
    if(triggertype != save_triggerType && save_triggerType != kTriggerUndefined)
      continue;
    float         triggertime = TriggerTimes[itrigger];
    std::vector<Float_t> triggerinfo = TriggerInfos[itrigger];
    float lowerbound = triggertime + WCSimWCTriggerBase::eventgatedown;
    float upperbound = triggertime + WCSimWCTriggerBase::eventgateup;

#ifdef WCSIMWCTRIGGERBASE_VERBOSE
    G4cout << "Saving trigger " << itrigger << " of type " << WCSimEnumerations::EnumAsString(triggertype)
	   << " in time range [" << lowerbound << ", " << upperbound << "]"
	   << " with trigger time " << triggertime
	   << " and additional trigger info";
    for(std::vector<Float_t>::iterator it = triggerinfo.begin(); it != triggerinfo.end(); ++it)
      G4cout << " " << *it;
    G4cout << G4endl;
#endif

    //loop over PMTs
    for (G4int i = 0; i < WCDCPMT->entries(); i++) {
      int tube=(*WCDCPMT)[i]->GetTubeID();
      //loop over digits
      for ( G4int ip = 0; ip < (*WCDCPMT)[i]->GetTotalPe(); ip++){
	int digit_time  = (*WCDCPMT)[i]->GetTime(ip);
	if(digit_time >= lowerbound && digit_time <= upperbound) {
	  //hit in event window
	  //add it to DigitsCollection
	  float peSmeared = (*WCDCPMT)[i]->GetPe(ip);
	  float Q = (peSmeared > 0.5) ? peSmeared : 0.5;
	  G4double digihittime = -triggertime
	    + WCSimWCTriggerBase::offset
	    + digit_time
	    + PMT->HitTimeSmearing(Q);
	  if(digihittime < 0)
	    continue;

	  //int parentID    = (*WCDCPMT)[i]->GetParentID(ip);
	  std::vector< std::pair<int,int> > digitized_composition = (*WCDCPMT)[i]->GetDigiCompositionInfo();
	  std::vector< std::pair<int,int> > triggered_composition;
	  for(std::vector< std::pair<int,int> >::iterator it = digitized_composition.begin(); it != digitized_composition.end(); ++it) {
	    if((*it).first == ip) {
	      triggered_composition.push_back(std::make_pair(itrigger, (*it).second));
	    }
	    else if ((*it).first > ip)
	      break;
	  }//loop over digitized_composition
	  //add hit
	  if ( DigiHitMap[tube] == 0) {
	    WCSimWCDigi* Digi = new WCSimWCDigi();
	    Digi->SetTubeID(tube);
	    //Digi->AddParentID(parentID);
	    Digi->AddGate  (itrigger,triggertime);
	    Digi->SetTime  (itrigger,digihittime);
	    Digi->SetPe    (itrigger,peSmeared);
	    Digi->AddPe    (digihittime);
	    Digi->AddDigiCompositionInfo(triggered_composition);
	    DigiHitMap[tube] = DigitsCollection->insert(Digi);
	  }
	  else {
	    //(*DigitsCollection)[DigiHitMap[tube]-1]->AddParentID(parentID);
	    (*DigitsCollection)[DigiHitMap[tube]-1]->AddGate(itrigger, triggertime);
	    (*DigitsCollection)[DigiHitMap[tube]-1]->SetTime(itrigger, digihittime);
	    (*DigitsCollection)[DigiHitMap[tube]-1]->SetPe  (itrigger, peSmeared);
	    (*DigitsCollection)[DigiHitMap[tube]-1]->AddPe  (digihittime);
	    (*DigitsCollection)[DigiHitMap[tube]-1]->AddDigiCompositionInfo(triggered_composition);
	  }
	  if(remove_hits)
	    (*WCDCPMT)[i]->RemoveDigitizedGate(ip);
	}//digits within trigger window
      }//loop over Digits
    }//loop over PMTs
  }//loop over Triggers
  G4cout << "WCSimWCTriggerBase::FillDigitsCollection. Number of entries in output digit collection: " << DigitsCollection->entries() << G4endl;

}


void WCSimWCTriggerBase::AlgNHitsThenSubNHits(WCSimWCDigitsCollection* WCDCPMT, bool remove_hits)
{
 

  nhitsVTXmap.clear();
  nhitsmap.clear();
  vtxVector.clear();
 
 
  // Get the info for pmt positions                                                                                                                                              
  std::vector<WCSimPmtInfo*> *pmts = myDetector->Get_Pmts();

  //Now we will try to find triggers
  //loop over PMTs, and Digits in each PMT. 
  // If nhits > Threshhold in a time window, then we have a trigger
  // If this cut fails, attempt an ITC ratio trigger

  int ntrig = 0;
  int window_start_time = 0;
  int window_end_time   = WCSimWCTriggerBase::LongTime - nhitsWindow;
  int window_step_size  = 5; //step the search window along this amount if no trigger is found
  float lasthit;
  std::vector<int> digit_times;
  //   std::vector<int> digit_times_vtx[100][100][100];
  
  //int digit_times_vtx[100][100][100][2000];
  std::vector<int> digit_times_itc_small, digit_times_itc_large;
  bool first_loop = true;

  G4cout << "WCSimWCTriggerBase::AlgNHits. Number of entries in input digit collection: " << WCDCPMT->entries() << G4endl;
#ifdef WCSIMWCTRIGGERBASE_VERBOSE
  int temp_total_pe = 0;
  for (G4int i = 0 ; i < WCDCPMT->entries() ; i++) {
    temp_total_pe += (*WCDCPMT)[i]->GetTotalPe();
  }
  G4cout << "WCSimWCTriggerBase::AlgNHits. " << temp_total_pe << " total p.e. input" << G4endl;
#endif

  //Generate verticies
  //skip for now to make the algorithm

  //generate vertices automatically.
  //superK
  float radius,height;
  height = 36000./2.; //mm middle of tank == 0 so half the height
  radius = 33681.5/2.; //mm
  //generate a box of verticies.  Save only those within the radius given.
  //5m = 5000mm spacing
  //
  int n_vtx = 0;
  for(int i=-1*height;i <= height;i=i+5000) {
    for(int j=-1*radius;j<=radius;j=j+5000) {
      for(int k=-1*radius;k<=radius;k=k+5000) {
	if(sqrt(j*j+k*k) > radius)
	  continue;
	vtxVector.push_back(G4ThreeVector(j*1.,k*1.,i*1.));
	n_vtx++;
      }
    }
  }
 
//   vtxVector.push_back(G4ThreeVector(0.,0.,0.));
//   vtxVector.push_back(G4ThreeVector(0.,0.,5000.));
//   vtxVector.push_back(G4ThreeVector(0.,0.,-5000.));
//   vtxVector.push_back(G4ThreeVector(0.,5000.,0.));
//   vtxVector.push_back(G4ThreeVector(0.,-5000.,0.));
//   vtxVector.push_back(G4ThreeVector(5000.,0.,0.));
//   vtxVector.push_back(G4ThreeVector(-5000.,0.,0.));
//   vtxVector.push_back(G4ThreeVector(0.,0.,10000.));
//   vtxVector.push_back(G4ThreeVector(0.,0.,-10000.));
//   vtxVector.push_back(G4ThreeVector(0.,10000.,0.));
//   vtxVector.push_back(G4ThreeVector(0.,-10000.,0.));
//   vtxVector.push_back(G4ThreeVector(10000.,0.,0.));
//   vtxVector.push_back(G4ThreeVector(-10000.,0.,0.));
   nhitsmap.resize(int(floor(WCSimWCTriggerBase::LongTime/5)));
   nhitsVTXmap.resize(n_vtx);
   for(int i = 0;i<n_vtx;i++) 
     nhitsVTXmap.at(i).resize(int(floor(WCSimWCTriggerBase::LongTime/5)));

//    nhitsVTXmap.at(0).resize(int(floor(WCSimWCTriggerBase::LongTime/5)));
//    nhitsVTXmap.at(1).resize(int(floor(WCSimWCTriggerBase::LongTime/5)));
//    nhitsVTXmap.at(2).resize(int(floor(WCSimWCTriggerBase::LongTime/5)));
//    nhitsVTXmap.at(3).resize(int(floor(WCSimWCTriggerBase::LongTime/5)));
//    nhitsVTXmap.at(4).resize(int(floor(WCSimWCTriggerBase::LongTime/5)));
//    nhitsVTXmap.at(5).resize(int(floor(WCSimWCTriggerBase::LongTime/5)));
//    nhitsVTXmap.at(6).resize(int(floor(WCSimWCTriggerBase::LongTime/5)));
//    nhitsVTXmap.at(7).resize(int(floor(WCSimWCTriggerBase::LongTime/5)));
//    nhitsVTXmap.at(8).resize(int(floor(WCSimWCTriggerBase::LongTime/5)));
//    nhitsVTXmap.at(9).resize(int(floor(WCSimWCTriggerBase::LongTime/5)));
//    nhitsVTXmap.at(10).resize(int(floor(WCSimWCTriggerBase::LongTime/5)));
//    nhitsVTXmap.at(11).resize(int(floor(WCSimWCTriggerBase::LongTime/5)));
//    nhitsVTXmap.at(12).resize(int(floor(WCSimWCTriggerBase::LongTime/5)));
   Float_t hit_pos[3];
  
  float triggertime; //save each digit time, because the trigger time is the time of the first hit above threshold
  bool triggerfound = false;
 
  //internal to this routine we want to save trigger windows from the std Nhits trigger so to search for triggers 
  //using test vertices outside this window.

  //loop over all hits to fill maps in 5ns time windows.

   digit_times.clear();
   
   
  //speed of light in water
  double cLightH20 = CLHEP::c_light/1.3330; 
  
    for (G4int i = 0 ; i < WCDCPMT->entries() ; i++) {
      int tube=(*WCDCPMT)[i]->GetTubeID();
      //get PMT position
      WCSimPmtInfo* pmtinfo = (WCSimPmtInfo*)pmts->at( tube - 1 );
      //hit positions in cm, conveet to mm as the CLHEP untis are in mm and ns
      hit_pos[0] = 10*pmtinfo->Get_transx();
      hit_pos[1] = 10*pmtinfo->Get_transy();
      hit_pos[2] = 10*pmtinfo->Get_transz();
      
      //loop over vertices, calculate distance from vertex to PMT.
      for(int j=0;j<vtxVector.size();j++) {
	
	//		std::cout<<tube<<"HIT "<<j<<" "<<hit_pos[0]<<" "<<hit_pos[1]<<" "<<hit_pos[2]<<" "<<sqrt(hit_pos[0]*hit_pos[0]+hit_pos[1]*hit_pos[1])<<" "<<
	//vtxVector.at(j).x()<<" "<<vtxVector.at(j).y()<<" "<<vtxVector.at(j).z()<<" "<<CLHEP::c_light<<"\n";
	Double_t diff_x = hit_pos[0]-vtxVector.at(j).x();                                                                                                                  
	Double_t diff_y = hit_pos[1]-vtxVector.at(j).y();                                                                                                               
	Double_t diff_z = hit_pos[2]-vtxVector.at(j).z();
	//distance in mm
	Double_t dist_tube = sqrt(diff_x*diff_x+diff_y*diff_y+diff_z*diff_z);
	//time for light to travel this distance in ns
	Double_t tofPMT = dist_tube/cLightH20;
	//Loop over each Digit in this PMT                                                                                                         
	for ( G4int ip = 0 ; ip < (*WCDCPMT)[i]->GetTotalPe() ; ip++) {
	  int digit_time = (*WCDCPMT)[i]->GetTime(ip);
	  if(digit_time >= WCSimWCTriggerBase::LongTime)
	    continue;
	  int digit_time_5ns = digit_time/5.;
	  //only fill the non-corrected map once, the first time
	  if(j==0) {
	    nhitsmap[int(floor(digit_time_5ns))]++;
	  }
	  //std::cout<<"time tof "<<int(floor(digit_time_5ns-tofPMT/5.))<<"\n";
	  //time is shifted by 250ns in order to ensure non-zero times (correcting for speed of light)
	  if(int(floor(digit_time_5ns-tofPMT/5.+75)) < 0)
	    continue;
	  //	  std::cout<<"VTX "<<int(floor(digit_time_5ns-tofPMT/5.+75))<<" "<<digit_time_5ns<<" "<<tofPMT/5.<<"\n";
	  nhitsVTXmap.at(j)[int(floor(digit_time_5ns-tofPMT/5.+75))]++;
	 
	}    
      }
    }





  
// //     //Now time to find triggers in these maps
// //     //windowVTX
// //     //find normal triggers
      G4int acc = 0; // accumulated # hits within time window...
// //     std::map< G4int, G4int>::iterator _mGateKeeper, _mNextGate;
      G4float RealOffset;
// //     for( _mGateKeeper = nhitsmap.begin() ; _mGateKeeper != nhitsmap.end() ; _mGateKeeper++)
// //     {
// // 	acc = 0;
// // 	_mNextGate = _mGateKeeper;
// //         RealOffset = 0.0; 				// will need to add the offset later
// // 									// 40 means + 200ns
// // 									// so check 39 bins ahead in the histogram..
// // 	bool triggered = false;
// // 	while ( _mNextGate != nhitsmap.lower_bound( _mGateKeeper->first + 39)
// // 		&& _mNextGate->first <= _mGateKeeper->first + 39 		// but not more than 200ns away though!
// // 	      )
// // 	{
// // 	  acc += _mNextGate->second;
// //           if (!triggered &&  acc > nhitsThreshold)
// // 	  {
// // 	    //RealOffset = _mGateKeeper->first*5.0;
// // 	    RealOffset = _mNextGate->first*5.0;
// // 	    	    TriggerTimes.push_back(RealOffset);
// // 	    TriggerTypes.push_back(kTriggerNHitsSKDETSIM);
// // 	    //std::cerr << "found a trigger..." << RealOffset/5.0  <<"\n";
// // 	    _mGateKeeper = nhitsmap.lower_bound( _mNextGate->first + G4int(WCSimWCTriggerBase::eventgateup )/5. );
// // 	    std::cerr.flush();
// // 	    triggered = true;
// //           }
// // 	  _mNextGate++;							// look at the next time bin with hits
// // 	}//while
// // 	if ( acc > nhitsThreshold)
// //           TriggerInfos.push_back(std::vector<Float_t>(1, acc));
// //     }

    //Now we want to use the tof corrected times to try to find other triggers
    //instead of a 200ns window we will use a smaller one 30ns.
    //This will have to be done in a slightly different manner looping over all 5ns time windows from
    //0 to Longtime.  This is due to wanting to check the same time window for each test vertex
    
    int count = 0;
    std::vector<std::pair<int,int> > PossibleTrigger;
    std::vector<int> PossibleTriggerCount;
    //loop from time = 0 to 20ns before the Longtime value
    while(count <= WCSimWCTriggerBase::LongTime/5 - 5) {
      int maxnpmt=0;
      int vtxindex=-1;
      for(int j=0;j<vtxVector.size();j++) {
	acc = 0;
	int lowerbound = count;
	int upperbound = count + 4; //ie. 20ns
	//sum the number of hit PMTs in this time window
	int npmt = 0;
	for(int k=lowerbound;k<upperbound;k++) {
	  npmt = npmt + nhitsVTXmap.at(j).at(k);	  
	}
       
	if(npmt > maxnpmt) {
	  maxnpmt = npmt;
	  vtxindex = j;
	}      
      }
      //      std::cout<<"MAX "<<vtxindex<<" "<<count <<" "<<count+3<<" "<<maxnpmt<<"\n";
      if(maxnpmt > nhitsThreshold) {//trigger
	std::cout<<"MAX "<<vtxindex<<" "<<count <<" "<<count+3<<" "<<maxnpmt<<" vtx x "<<vtxVector.at(vtxindex).x()<<
	  " y "<<vtxVector.at(vtxindex).y()<<" z "<<vtxVector.at(vtxindex).z()<<"\n";
	count = count + 4;
	PossibleTrigger.push_back(std::make_pair(vtxindex,count*5));
	PossibleTriggerCount.push_back(maxnpmt);
      }
      else {
	count++;
      }
    }
    //sort the possible triggers by time  
    std::sort(PossibleTrigger.begin(), PossibleTrigger.end(), sort_pair_second<int, int>());
    
    //loop over triggers from earliest to latest
    std::vector<std::pair<int,int> > TriggerPairsCorT;
    std::vector<std::pair<int,int> > TriggerNormT;
    int time_start=-9999;

    for(int i=0;i<PossibleTrigger.size();i++) {
      //once a trigger is found, we must jump to 950ns in the future before searching for the next
      if(PossibleTrigger.at(i).second > time_start) {
	TriggerPairsCorT.push_back(std::make_pair(PossibleTrigger.at(i).first,PossibleTrigger.at(i).second));
	time_start = PossibleTrigger.at(i).second + WCSimWCTriggerBase::eventgateup;
    
	//set trigger details to be passed onto fill routine
	//trigger time in non-corrected times is timeTrig.
	int triggertime = PossibleTrigger.at(i).second; //correct for the 200ns offset added earlier
	triggertime -= (int)triggertime % 5;
	std::vector<Float_t> triggerinfo;
	triggerinfo.push_back(vtxVector.at((Float_t)PossibleTrigger.at(i).first).x());
	triggerinfo.push_back(vtxVector.at((Float_t)PossibleTrigger.at(i).first).y());
	triggerinfo.push_back(vtxVector.at((Float_t)PossibleTrigger.at(i).first).z());
	TriggerTypes.push_back(kTriggerSubNHits);
	TriggerInfos.push_back(triggerinfo);
	TriggerTimes.push_back(triggertime);
      }
    }

    FillDigitsCollection(WCDCPMT, remove_hits, kTriggerSubNHits);
}
