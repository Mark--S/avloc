//
// Basic event processor for AV locations
// S.J.M.Peeters@sussex.ac.uk, June 2014
//
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <RAT/Log.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/MCPMT.hh>
#include <RAT/DS/MCPE.hh>
#include "include/AVLocProc.h"
#include "include/AVLocTools.h"
#include <cstdio>

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
    for( int ipmt = 0; ipmt < pmtList.GetCount(); ipmt++) {
      int PMTID = pmtList.GetPMT(ipmt).GetID();
      TVector3 dist(pmt_info.x_pos[PMTID],pmt_info.y_pos[PMTID],pmt_info.z_pos[PMTID]);
      //printf("pmt position %f %f %f\n",dist.X(),dist.Y(),dist.Z());
      //printf("led X:%f ledy:%f ledz:%f\n",led_info.position.X(),led_info.position.Y(),led_info.position.Z());
      dist -= led_info.position;
      //printf("difference %f %f %f\n",dist.X(),dist.Y(),dist.Z());
      //printf("Trigger Delay %f\n",GTTriggerDelay);
      double EVoffset = 500 - GTTriggerDelay- 100; 
      EVoffset = 0;
      //cout << "HIT TAC: "<<pmtList.GetPMT(ipmt).GetTime()<<endl;
      //printf("Universal time %f Universal Time Days %u Universal Time Seconds %u  Clock Ticks:%llu EVOffset %f\n",rEV.GetUniversalTime().GetNanoSeconds(),rEV.GetUniversalTime().GetDays(),rEV.GetUniversalTime().GetSeconds(),rEV.GetClockCount50(),EVoffset);
      Double_t PMTTime = pmtList.GetPMT(ipmt).GetTime()-EVoffset;
      //cout << "PMT TIME: "<<PMTTime<<endl;
      if(PMTTime>-400){
          ntuple->Fill(led_info.nr,led_info.sub,PMTID,PMTTime,dist.Mag());
      }
    }
  }
  return true;
}

// process event and fill the avproc ntuple
bool ProcessEventNtupleMC(RAT::DS::Entry * rDS, TNtuple * ntuple, 
			LEDInfo & led_info, PMTInfo & pmt_info)
{
  double GTTriggerDelay = RAT::DB::Get()->GetLink("DAQ")->GetD("gtriggerdelay");	
  int detectorEventCount = rDS->GetEVCount(); 
  RAT::DS::MC mc = rDS->GetMC();
  if(detectorEventCount > 1){
      cout << "More than one detector event in tellie flash cutting" << endl;
      return true;
  };
  for( int iEV = 0; iEV < detectorEventCount; ++iEV) {
    RAT::DS::EV& rEV = rDS->GetEV(iEV);
    RAT::DS::CalPMTs& pmtList  = rEV.GetCalPMTs();
    for( int ipmt = 0; ipmt < pmtList.GetCount(); ++ipmt) {
      int PMTID = pmtList.GetPMT(ipmt).GetID();
      TVector3 dist(pmt_info.x_pos[PMTID],pmt_info.y_pos[PMTID],pmt_info.z_pos[PMTID]);
      //printf("pmt position %f %f %f\n",dist.X(),dist.Y(),dist.Z());
      //printf("led X:%f ledy:%f ledz:%f\n",led_info.position.X(),led_info.position.Y(),led_info.position.Z());
      dist -= led_info.position;
      //printf("difference %f %f %f\n",dist.X(),dist.Y(),dist.Z());
      RAT::DS::MCEV  mcEV = rDS->GetMCEV(iEV);
      double gtTime = mcEV.GetGTTime();
      //printf("Trigger Delay %f\n",GTTriggerDelay);
      double EVoffset = 500 - GTTriggerDelay - gtTime; 
      //cout << "EVoffset: "<<EVoffset<<endl;
      //printf("Universal time %f Universal Time Days %u Universal Time Seconds %u  Clock Ticks:%llu EVOffset %f\n",rEV.GetUniversalTime().GetNanoSeconds(),rEV.GetUniversalTime().GetDays(),rEV.GetUniversalTime().GetSeconds(),rEV.GetClockCount50(),EVoffset);
      Double_t PMTTime = pmtList.GetPMT(ipmt).GetTime()-EVoffset;
      //RAT::DS::MCPMT mcPMT = GetMCPMT(PMTID,mc);
      //int numPE = mcPMT.GetMCPECount();
      //if(numPE>1) continue;
      ntuple->Fill(led_info.nr,led_info.sub,PMTID,PMTTime,dist.Mag());
    }
  }
  return true;
}
//Method to get the corresponding MCPMT from the corresponding event
RAT::DS::MCPMT GetMCPMT(int pmtId, RAT::DS::MC& mc){
    int pmtCount = mc.GetMCPMTCount();
    for(int i=0; i<pmtCount; i++){
       RAT::DS::MCPMT testPMT = mc.GetMCPMT(i);
       if(testPMT.GetID()==pmtId) return testPMT;
    }
}

