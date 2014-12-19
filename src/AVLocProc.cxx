//
// Basic event processor for AV locations
// S.J.M.Peeters@sussex.ac.uk, June 2014
//
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <RAT/Log.hh>

#include "include/AVLocProc.h"
#include "include/AVLocTools.h"

using namespace std;

// process event and print something - just a first try
bool ProcessEventBasic(RAT::DS::Root * rDS)
{
  for( int iEV = 0; iEV < rDS->GetEVCount(); ++iEV) {
    RAT::DS::EV* rEV = rDS->GetEV(iEV);
    for( int ipmt = 0; ipmt < rEV->GetPMTCalCount(); ++ipmt) {
      int pmtID = rEV->GetPMTCal(ipmt)->GetID();
      Double_t PMTTime = rEV->GetPMTCal(ipmt)->GetTime();
      cerr << pmtID << ": " << PMTTime << endl;
    }
  }
  return true;
}

// process event and fill the avproc ntuple
bool ProcessEventNtuple(RAT::DS::Root * rDS, TNtuple * ntuple, 
			LEDInfo & led_info, PMTInfo & pmt_info)
{
  double GTTriggerDelay = RAT::DB::Get()->GetLink("DAQ")->GetD("gtriggerdelay");	
  for( int iEV = 0; iEV < rDS->GetEVCount(); ++iEV) {
    RAT::DS::EV* rEV = rDS->GetEV(iEV);
    for( int ipmt = 0; ipmt < rEV->GetPMTCalCount(); ++ipmt) {
      int PMTID = rEV->GetPMTCal(ipmt)->GetID();
      TVector3 dist(pmt_info.x_pos[PMTID],pmt_info.y_pos[PMTID],pmt_info.z_pos[PMTID]);
      dist -= led_info.position;
      double EVoffset = 500 - GTTriggerDelay - (rEV->GetGTrigTime());
      Double_t PMTTime = rEV->GetPMTCal(ipmt)->GetTime()-EVoffset;
      ntuple->Fill(led_info.nr,led_info.sub,PMTID,PMTTime,dist.Mag());
    }
  }
  return true;
}
