//
// Basic event processor for AV locations
// S.J.M.Peeters@sussex.ac.uk, June 2014
//
#ifndef __AVLOCPROC_H__
#define __AVLOCPROC_H__

#include <TNtuple.h>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/MC.hh>

#include "include/AVLocTools.h"

bool ProcessEventBasic(RAT::DS::Entry * rDS);
bool ProcessEventNtuple(RAT::DS::Entry * rDS, TNtuple * ntuple, 
			LEDInfo & led_info, PMTInfo & pmt_info);

bool ProcessEventNtupleMC(RAT::DS::Entry * rDS, TNtuple * ntuple, 
			LEDInfo & led_info, PMTInfo & pmt_info);

RAT::DS::MCPMT GetMCPMT(int pmtId, RAT::DS::MC& mc);
#endif
