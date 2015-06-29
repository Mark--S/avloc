//
// Library for AV Location software - toolkit
// S.J.M.Peeters@sussex.ac.uk, June 2014
//
#ifndef __AVLOCTOOLS_H__
#define __AVLOCTOOLS_H__
using namespace std;
#include <string>
#include <TH1D.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TVector3.h>
#include <RAT/DB.hh>
#include <RAT/DS/Entry.hh>
#include <TGraph.h>
#include <RAT/DS/Run.hh>

// function to load a root file
void LoadRootFile(string filename, TTree **tree, RAT::DS::Entry **rDS, RAT::DS::Run **rRun);

// function to load the RAT database
void LoadDataBase(string logname);

// Storage for LED information 
struct LEDInfo {
  string   name;
  int      nr;    // fibre number
  int      sub;   // 0 for A fibre, 1 for B fibre
  TVector3 position;
  TH1D *   spectrum;
};

struct PMTInfo {  
  vector<double> x_pos;
  vector<double> y_pos;
  vector<double> z_pos;
  vector<double> x_dir;
  vector<double> y_dir;
  vector<double> z_dir;
};

// Get PMT position - in RAT, this is the centre of the bucket, ie 5.67 cm in front of the PMT
// (source: SNO NIM)
PMTInfo GetPMTpositions(void);

// function to extract relevant LED info (uses part of filename and assumes db loaded)
LEDInfo GetLEDInfoFromFileName(string filename);
LEDInfo GetLEDInfoFromFibreNr(int fibre_nr, int subnr); // subnr: 0 for A fibre 1 for B fibre
LEDInfo GetLEDInfoFromFibreName(string fibre_name);

// function to gain access to a (re)writable ntuple for summary data for avloc
TNtuple * GetNtuple(TFile ** fpointer, TString filename);

struct PhysicsNr {
  double value;
  double error;
};

// function to convert LED spectrum convoluted with the group velocity in water into an average and an error
PhysicsNr GroupVelocity(string fibre_name);

// Calculate time of flight in ns
PhysicsNr TimeOfFlight(TVector3 inject, TVector3 detect, PhysicsNr n_h2o, double offset);

#endif
