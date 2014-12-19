//
// Basic event processor for AV locations
// S.J.M.Peeters@sussex.ac.uk, June 2014
//
#ifndef __AVLOCPROC_H__
#define __AVLOCPROC_H__

#include <TNtuple.h>
#include <RAT/DS/Root.hh>

#include "include/AVLocTools.h"

bool ProcessEventBasic(RAT::DS::Root * rDS);
bool ProcessEventNtuple(RAT::DS::Root * rDS, TNtuple * ntuple, 
			LEDInfo & led_info, PMTInfo & pmt_info);

#endif
