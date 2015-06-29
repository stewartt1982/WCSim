#include "WCSimDarkRateMessenger.hh"
#include "WCSimWCDigitizer.hh"
#include "WCSimWCAddDarkNoise.hh"

#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"

WCSimDarkRateMessenger::WCSimDarkRateMessenger(WCSimWCDigitizer* darkratepoint)
  : WCSimDigitize(darkratepoint)
{
  Initalise();

  //inform the DarkRateMessenger which constructor was called
  constructor = 0;
}

WCSimDarkRateMessenger::WCSimDarkRateMessenger(WCSimWCAddDarkNoise* darkratepoint)
  : WCSimAddDarkNoise(darkratepoint)
{
  Initalise();
  
  //inform the DarkRateMessenger which constructor was called   
  constructor = 1;
}

void WCSimDarkRateMessenger::Initalise()
{
  WCSimDir = new G4UIdirectory("/DarkRate/");
  WCSimDir->SetGuidance("Commands to change the dark noise frequency of the simulation");
  
  SetFrequency = new G4UIcmdWithADoubleAndUnit("/DarkRate/SetDarkRate",this);
  SetFrequency->SetGuidance("Commands to change the dark noise frequency of the simulation");
  SetFrequency->SetParameterName("DarkRate",true);
  SetFrequency->SetDefaultValue(0);
  //kilohertz is 10e-6
  SetFrequency->SetUnitCategory("Frequency");
  SetFrequency->SetDefaultUnit("kHz");
  SetFrequency->SetUnitCandidates("Hz kHz MHz GHz");

  SetConversionRate = new G4UIcmdWithADouble("/DarkRate/SetConvert",this);
  SetConversionRate->SetGuidance("Caribrate the frequency of dark noise after digitization");
  SetConversionRate->SetParameterName("DigiCorr",true);
  SetConversionRate->SetDefaultValue(1);


  //Mode 0 - Add dark rate in window defined by /DarkRate/SetDarkLow and /DarkRate/SetDarkHigh
  //If not set default is 0 and 100000ns
  //Mode 1 - Add dark rate to a window of size /DarkRate/SetDarkWindow around each hit
  SetDarkMode = new G4UIcmdWithAnInteger("/DarkRate/SetDarkMode",this);
  SetDarkMode->SetGuidance("Set the mode of adding dark noise to the event");
  SetDarkMode->SetParameterName("DarkMode",true);
  SetDarkMode->SetDefaultValue(0);
  
  SetDarkLow = new G4UIcmdWithADoubleAndUnit("/DarkRate/SetDarkLow",this);
  SetDarkLow->SetGuidance("Set the lower GEANT time limit to add dark noise");
  SetDarkLow->SetParameterName("DarkLow",true);
  SetDarkLow->SetDefaultValue(0.);
  SetDarkLow->SetUnitCategory("Time");
  SetDarkLow->SetDefaultUnit("ns");

  SetDarkHigh = new G4UIcmdWithADoubleAndUnit("/DarkRate/SetDarkHigh",this);
  SetDarkHigh->SetGuidance("Set the upper GEANT time limit to add dark noise");
  SetDarkHigh->SetParameterName("DarkHigh",true);
  SetDarkHigh->SetDefaultValue(100000.);
  SetDarkHigh->SetUnitCategory("Time");
  SetDarkHigh->SetDefaultUnit("ns");

  SetDarkWindow = new G4UIcmdWithADoubleAndUnit("/DarkRate/SetDarkWindow",this);
  SetDarkWindow->SetGuidance("Set the window width to add dark noise");
  SetDarkWindow->SetParameterName("DarkLow",true);
  SetDarkWindow->SetDefaultValue(5000.);
  SetDarkWindow->SetUnitCategory("Time");
  SetDarkWindow->SetDefaultUnit("ns");
}

WCSimDarkRateMessenger::~WCSimDarkRateMessenger()
{

  delete SetFrequency;
  delete SetConversionRate;
  delete SetDarkMode;
  delete SetDarkLow;
  delete SetDarkHigh;
  delete SetDarkWindow;
  delete WCSimDir;
}

void WCSimDarkRateMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{

  if(command == SetFrequency){
    // Since kHz is 10e-3 for this class we must multiply by a 10e6 factor
    // to make default units in kHz
    double const conversion_to_kHz = 1000000;
    if(constructor == 0)
      WCSimDigitize->SetDarkRate(conversion_to_kHz * SetFrequency->GetNewDoubleValue(newValue));
    else if(constructor == 1)
      WCSimAddDarkNoise->SetDarkRate(conversion_to_kHz * SetFrequency->GetNewDoubleValue(newValue));
    G4cout << "Setting Dark Rate " << conversion_to_kHz * SetFrequency->GetNewDoubleValue(newValue) << G4endl;
  }

  if(command == SetConversionRate){
    if(constructor == 0)
      WCSimDigitize->SetConversion(SetConversionRate->GetNewDoubleValue(newValue));
    else if(constructor == 1)
      WCSimAddDarkNoise->SetConversion(SetConversionRate->GetNewDoubleValue(newValue));
    G4cout << "Setting conversion value " << SetConversionRate->GetNewDoubleValue(newValue) << G4endl;
  }

  if(command == SetDarkMode){
    if(constructor == 0){
      //      WCSimDigitize->SetDarkMode(SetDarkMode->GetNewDoubleValue(newValue));
      }
    else if(constructor == 1) {
      WCSimAddDarkNoise->SetDarkMode(SetDarkMode->GetNewIntValue(newValue));
      G4cout << "Setting DarkMode value " << SetDarkMode->GetNewIntValue(newValue) << G4endl;
    }
  }

  if(command == SetDarkLow){
    if(constructor == 0){
      //      WCSimDigitize->SetDarkLow(SetDarkLow->GetNewDoubleValue(newValue));
    }   
    else if(constructor == 1) {
      WCSimAddDarkNoise->SetDarkLow(SetDarkLow->GetNewDoubleValue(newValue));
      G4cout << "Setting DarkLow value " << SetDarkLow->GetNewDoubleValue(newValue) << G4endl;
    }
  }

  if(command == SetDarkHigh){
    if(constructor == 0) {
      //      WCSimDigitize->SetDarkHigh(SetDarkHigh->GetNewDoubleValue(newValue));
    }   
    else if(constructor == 1) {
      WCSimAddDarkNoise->SetDarkHigh(SetDarkHigh->GetNewDoubleValue(newValue));
      G4cout << "Setting DarkHigh value " << SetDarkHigh->GetNewDoubleValue(newValue) << G4endl;
    }
  }

  if(command == SetDarkWindow){
    if(constructor == 0){
      //      WCSimDigitize->SetDarkWidth(SetDarkWidth->GetNewDoubleValue(newValue));
      }
    else if(constructor == 1) {
      WCSimAddDarkNoise->SetDarkWindow(SetDarkWindow->GetNewDoubleValue(newValue));
      G4cout << "Setting DarkWidth value " << SetDarkWindow->GetNewDoubleValue(newValue) << G4endl;
    }
  }

}
