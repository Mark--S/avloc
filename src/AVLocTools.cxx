//
// Library for AV Location software - toolkit
// S.J.M.Peeters@sussex.ac.uk, June 2014
//
#include <assert.h>
#include <iostream>

#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <RAT/Log.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DU/GroupVelocity.hh>
#include <RAT/DB.hh>

#include "include/AVLocTools.h"
// function to load a root file
void LoadRootFile(string filename, TTree **tree, RAT::DS::Entry **rDS, RAT::DS::Run **rRun)
{
  TFile *file = new TFile(filename.data());
  (*tree) = (TTree*)file->Get( "T" );
  TTree *runTree = (TTree*)file->Get("runT");
  assert(runTree);
  *rDS = new RAT::DS::Entry();
  (*tree)->SetBranchAddress( "ds", &(*rDS) );
  assert(rDS);
  *rRun = new RAT::DS::Run();
  assert(rRun);
  runTree->SetBranchAddress( "run", &(*rRun) );
  runTree->GetEntry();
}

// function to load database
void LoadDataBase(string logname)
{
  RAT::Log::Init(logname.data());
  RAT::DB * db = RAT::DB::Get();
  assert(db);	
  char* glg4data = getenv("GLG4DATA");
  if (glg4data == static_cast<char*>(NULL)) {
    cerr << "ratzdab::ratdb: Environment variable $GLG4DATA must be set" << endl;
    assert(glg4data);
  }
  string data = string(glg4data);
  cout << "LoadDataBase: loading defaults..." << endl;
  db->LoadDefaults();
}

PMTInfo GetPMTpositions(void) {
  cout << "Loading PMT positions" << endl;
  const double offset = 56.7; // difference front PMT and bucket in mm
  RAT::DB* db = RAT::DB::Get();
  assert(db);
  char* ratroot = getenv("RATROOT");
  if (ratroot == static_cast<char*>(NULL)) {
      cerr << "Environment variable $RATROOT must be set" << endl;
      assert(ratroot);
  }
  string rat     = string(ratroot);
    string pmtfile = rat;
  pmtfile += "/data/pmt/snoman.ratdb";
  db->LoadFile(pmtfile);
  RAT::DBLinkPtr pmtInfo = db->GetLink("PMTINFO");
  assert(pmtInfo);
  PMTInfo pmt_info;
  cout << "Getting Arrays" << endl;
  pmt_info.x_pos = pmtInfo->GetDArray("x");
  pmt_info.y_pos = pmtInfo->GetDArray("y");
  pmt_info.z_pos = pmtInfo->GetDArray("z");
  cout << "Got Position arrays" << endl;
  vector<double> xDir = pmtInfo->GetDArray("v");
  vector<double> yDir = pmtInfo->GetDArray("u");
  vector<double> zDir = pmtInfo->GetDArray("w");
  cout << "Obtained ARRAYS" << endl;
  for ( unsigned int i = 0 ; i < xDir.size() ; ++i ) {
    TVector3 pos(pmt_info.x_pos[i],pmt_info.y_pos[i],pmt_info.z_pos[i]);
    TVector3 dir(xDir[i],yDir[i],zDir[i]);
    cerr << pos.Mag() << " -> "; 
    pos -= dir.Unit()*offset;
    pmt_info.x_pos[i] = pos.X();
    pmt_info.y_pos[i] = pos.Y();
    pmt_info.z_pos[i] = pos.Z();
    cerr << pos.Mag() << endl;
  }
  return pmt_info;
} 

// function to extract relevant LED info (uses part of filename and assumes db loaded)
LEDInfo GetLEDInfoFromFileName(string filename)
{ 
  string shortname  = filename.substr(filename.find_last_of('/')+1);
  std::cout << "shortname " << shortname << std::endl;
  string fibre_name = shortname.substr(6, 6);
  std::cout << "fibrename " << fibre_name << std::endl;
  return GetLEDInfoFromFibreName(fibre_name);
}


LEDInfo GetLEDInfoFromFibreNr(int fibre_nr, int subnr) 
{
  string fibre_name = "FT";
  char   nr[8];
  sprintf(nr,"%03i",fibre_nr);
  fibre_name += nr;
  assert (subnr == 0 || subnr == 1);
  subnr == 0 ? fibre_name += 'A' : fibre_name += 'B';
  return GetLEDInfoFromFibreName(fibre_name);
}

LEDInfo GetLEDInfoFromFibreName(string fibre_name) 
{
  RAT::DB * db = RAT::DB::Get();
  assert(db);
  LEDInfo led_info;
  led_info.name = fibre_name;
  
  // get number
  //std::cout << led_info.name << std::endl;
  string nr = led_info.name.substr(2,3);
  led_info.nr = atoi(nr.data());
  char letter = led_info.name[5];
  //std::cout << letter << std::endl;
  if      ( letter == 'A' ) { led_info.sub = 0; }
  else if ( letter == 'B' ) { led_info.sub = 1; }
  else { cerr << "Unknown sub fibre: " << letter << " in " << led_info.name << endl; }

  // get position  
  RAT::DBLinkPtr led_db = db->GetLink("FIBRE",led_info.name.data());
  assert(led_db);
  led_info.position.SetX(led_db->GetD("x"));
  led_info.position.SetY(led_db->GetD("y"));
  led_info.position.SetZ(led_db->GetD("z"));

  // get wavelength spectrum
  vector<Float_t> amp = db->GetLink("ELLIEWAVE","TELLIE503")->GetFArrayFromD("dist_wl_intensity");
  vector<Float_t> wl  = db->GetLink("ELLIEWAVE","TELLIE503")->GetFArrayFromD("dist_wl");
  assert(wl.size() == amp.size() );
  char title[128];
  sprintf(title,"%s (%.2f,%.2f,%.2f)",led_info.name.data(),
	  led_info.position.X(),led_info.position.Y(),led_info.position.Z());
  int nbins = wl.size();
  float min = wl.front();
  float max = wl.back();
  assert ( max > min );
  assert ( nbins != 0 );
  double bw  = (max-min)/(float)nbins;
  try{
      delete gROOT->FindObject("hLEDData");
      led_info.spectrum = new TH1D("hLEDData",title,nbins,min-0.5*bw,max+0.5*bw);
  }
  catch(int e){
      led_info.spectrum = new TH1D("hLEDData",title,nbins,min-0.5*bw,max+0.5*bw);
  }
  led_info.spectrum->SetXTitle("wavelength (nm)");
  for (unsigned int i = 0 ; i < wl.size() ; ++i) {
    led_info.spectrum->SetBinContent(i+1,amp[i]);
  }

  // report 
  /*
  cout << "LED: " << led_info.name;
  cout << " (" << led_info.nr << " - " << led_info.sub  << ") @ (";
  cout << led_info.position.X() << ",";
  cout << led_info.position.Y() << ",";
  cout << led_info.position.Z() << ")";
  cout << "; wavelength is (" << led_info.spectrum->GetMean() << " +/- ";
  cout << led_info.spectrum->GetRMS() << " nm)" << endl;
  */

  return led_info;
}

// function to gain access to a (re)writable ntuple for summary data for avloc
// content:
// 0 - fibre_nr  : fibre number
// 1 - fibre_sub : sub-fibre ID - 0 is fibre A, 1 is fibre B
// 2 - lcn       : logical channel number for the PMT
// 3 - time      : hit time (ns)
// 4 - dist      : distance from fibre (mm)
TNtuple * GetNtuple(TFile ** fpointer,TString outputFile)
{
  TNtuple * ntuple = NULL;
  TString filename = "summary_ntuple.root";
  filename+=outputFile;
  if ( gSystem->FindFile("./",filename) == 0 ) {
    cout << "CREATING NEW FILE " << filename << endl;
    (*fpointer) = new TFile(filename,"NEW");
    ntuple = new TNtuple("avloctuple","avloctuple","fibre_nr:fibre_sub:lcn:time:dist");
  }
  else {
    cout << "REUSING OLD FILE " << filename << endl;
    (*fpointer) = new TFile(filename,"UPDATE");
    ntuple = (TNtuple*)(*fpointer)->Get("avloctuple");
    assert(ntuple);
  }
  return ntuple;
}

// helper function to extract a y value from a graph, given x, 
// using linear interpolation between points
double get_value(TGraph * g, double x)
{
  int np = g->GetN();
  double * vx = g->GetX();
  double * vy = g->GetY();
  if ( x < vx[0] || x > vx[np-1] ) {
    cerr << "AVLocTools::get_value : x value out of range :" 
    << x << " vs [" << vx[0] << "," << vx[np-1] << "]" << endl;
    
    return 0.;
  }
  double value = 0.;
  for (int i = 0 ; i < np ; ++i ) {
    if ( vx[i] > x ) {
      value = vy[i-1] + (vx[i]-x)*(vy[i]-vy[i+1])/(vx[i]-vx[i+1]);
      break;
    }
  }
  return value;
}

// function to convert LED spectrum convoluted with the group velocity in water 
// into an average and an error
// for the weighted average calculation in one loop, see:
// http://en.wikipedia.org/wiki/Mean_square_weighted_deviation
PhysicsNr GroupVelocity(string fibre_name) 
{
  LEDInfo   led_info = GetLEDInfoFromFibreName(fibre_name);
  int nbins  = led_info.spectrum->GetNbinsX();
  double sumw  = 0.; // sum_i (w_i)
  double sum   = 0.; // sum_i (w_i * x_i)
  double sumsq = 0.; // sum_i (w_i * x_i^2)
  //Unsure of of 10^-4 scaling factor  required for calc by distance method as 400nm->3.103125 * 1e-6 units of energy
  double hc = 0.197*6.28318530718*10e-4;
  //const double f = h/e*c*1E-6*1E9; //  eV*nm = 1E-6 MeV * 1E9 nm
  //const RAT::DU::GroupVelocity& vel =  RAT::DU::Utility::Get()->GetGroupVelocity();
  for ( int i=1 ; i <= nbins ; ++i ){
    double wi = led_info.spectrum->GetBinContent(i);
    double energy = hc/led_info.spectrum->GetBinCenter(i);
    //printf("wavlength:%f energy:%f\n", led_info.spectrum->GetBinCenter(i),energy);
    //cout << "energy: "<< energy << endl;
    const double time = RAT::DU::Utility::Get()->GetGroupVelocity().CalcByDistance(0.0,0.0,1.0,energy);
    //cout.precision(15);
    //cout << "time: " << time << endl;
    double xi = 1.0/time;
    sumw  += wi;
    sum   += wi*xi;
    sumsq += wi*xi*xi;
  }
  PhysicsNr vgroup;
  if ( sumw == 0. ) {
    cerr << "AVLocTools::GroupVelocity : sum of weights is zero" << endl;
    vgroup.value = 0.;
    vgroup.error = 0.;
  }
  else {
    vgroup.value = sum/sumw;
    vgroup.error = (sumsq*sumw-sum*sum)/(sumw*sumw);
  }
  return vgroup;
}

PhysicsNr TimeOfFlight(TVector3 inject, TVector3 detect, PhysicsNr n_h2o, double offset) {
  PhysicsNr tof;
  assert(n_h2o.value);
  double vv = 299792458000./n_h2o.value;
  double ev = vv*(n_h2o.error/n_h2o.value);
  TVector3 reflect = (inject+detect).Unit() * 6050.;
  tof.value = ((inject - reflect*offset).Mag() + (detect - reflect*offset).Mag())/vv*1E9; // to ns
  tof.error = ev/vv*tof.value;
  return tof;	
}



