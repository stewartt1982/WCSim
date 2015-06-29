#ifndef WCSimWCTriggerNHitsThenSubNHits_h
#define WCSimWCTriggerNHitsThenSubNHits_h 1

#include "WCSimWCTriggerBase.hh"

class WCSimWCTriggerNHitsThenSubNHits : public WCSimWCTriggerBase
{
public:

  //not recommended to override these methods
  WCSimWCTriggerNHitsThenSubNHits(G4String name, WCSimDetectorConstruction*, WCSimWCDAQMessenger*);
  ~WCSimWCTriggerNHitsThenSubNHits();
  
private:
  void DoTheWork(WCSimWCDigitsCollection* WCDCPMT);

};

#endif








