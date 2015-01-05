//
// Basic event processor for AV locations
// S.J.M.Peeters@sussex.ac.uk, June 2014
//
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <RAT/Log.hh>
#include <RAT/DS/PMT.hh>
#include "include/AVLocProc.h"
#include "include/AVLocTools.h"

using namespace std;

// process event and print something - just a first try
bool ProcessEventBasic(RAT::DS::Entry * rDS)
{
  for( int iEV = 0; iEV < rDS->GetEVCount(); ++iEV) {
    RAT::DS::EV&  rEV = rDS->GetEV(iEV);
    RAT::DS::CalPMTs& pmtList  = rEV.GetCalPMTs();
    for( int ipmt = 0; ipmt < pmtList.GetCount(); ++ipmt) {
      int pmtID = pmtList.GetPMT(ipmt).GetID();
      Double_t PMTTime = pmtList.GetPMT(ipmt).GetTime();
      cerr << pmtID << ": " << PMTTime << endl;
    }
  }
  return true;
}

// process event and fill the avproc ntuple
bool ProcessEventNtuple(RAT::DS::Entry * rDS, TNtuple * ntuple, 
			LEDInfo & led_info, PMTInfo & pmt_info)
{
  double GTTriggerDelay = RAT::DB::Get()->GetLink("DAQ")->GetD("gtriggerdelay");	
  for( int iEV = 0; iEV < rDS->GetEVCount(); ++iEV) {
    RAT::DS::EV& rEV = rDS->GetEV(iEV);
    RAT::DS::CalPMTs& pmtList  = rEV.GetCalPMTs();
    for( int ipmt = 0; ipmt < pmtList.GetCount(); ++ipmt) {
      int PMTID = pmtList.GetPMT(ipmt).GetID();
      TVector3 dist(pmt_info.x_pos[PMTID],pmt_info.y_pos[PMTID],pmt_info.z_pos[PMTID]);
      dist -= led_info.position;
      double EVoffset = 500 - GTTriggerDelay - (MHzTicks2Seconds(rEV.GetClockCount50(),50)*1000000000);
      Double_t PMTTime = pmtList.GetPMT(ipmt).GetTime()-EVoffset;
      ntuple->Fill(led_info.nr,led_info.sub,PMTID,PMTTime,dist.Mag());
    }
  }
  return true;
}

//Taken from DataQualityProc.cc
double MHzTicks2Seconds( unsigned long int ticks, int frequency ) 
{
    double seconds = 1.0 / static_cast<double>( frequency );
    seconds *= 1.0e-6;
    seconds *= static_cast<double>( ticks );
    return seconds;
}

