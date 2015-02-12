// 
// Process TELLIE data for AV reflection fit
// S.J.M.Peeters@sussex.ac.uk, June 2014
//
#include <iostream>
#include <string>
#include <cmath>

#include <TBenchmark.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TFile.h>
#include <TTree.h>

#include <RAT/DB.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DU/GroupVelocity.hh>
#include <RAT/DU/Utility.hh>

#include "include/AVLocTools.h"
#include "include/AVLocProc.h"

using namespace std;

int main(int argc,char **argv)
{
  TBenchmark benchmark;
  benchmark.Start("MAKENTUPLE");
  string filename;
  if ( argc != 2 ) {
    cerr << "Usage: " << argv[0] << " <filename>" << endl;
    return 1;
  }
  else{
    filename = argv[1];
  }
  RAT::DS::Entry * rDS  = NULL;
  RAT::DS::Run  * rRun = NULL;
  TTree         * tree = NULL;
  LoadRootFile(filename,&tree,&rDS,&rRun);
  assert(rDS);
  assert(rRun);
  assert(tree);  
  LoadDataBase("make_ntuple.log");
  LEDInfo led_info = GetLEDInfoFromFileName(filename);
  //cout << "FILENAME IS: "<<filename<<endl;
  unsigned first = filename.find("/");
  unsigned last = filename.find(".");
  string filenameNew = filename.substr(first+1,last-first-1);
  TFile * ntuple_file = NULL;
  TNtuple * ntuple = GetNtuple(&ntuple_file,filenameNew);
  assert(ntuple_file);
  assert(ntuple);

  // get dispersion graph and change x values from MeV to nm
  // y values are presumable the speed in mm/ns
  const double c = 299792458;
  const double h = 6.62606957e-34;
  const double e = 1.602176565e-19;
  const double f = h/e*c*1E-6*1E9; //  eV*nm = 1E-6 MeV * 1E9 nm
  /*TGraph* AVGraph = rRun->GetGroupVelocityTime()->GetWaterGraph();
  double * x = AVGraph->GetX();
  double * y = AVGraph->GetY();
  int np = AVGraph->GetN();
  TGraph * water_group_velocity = new TGraph();
  int cntr = 0;
  for (int i = np-1; i > 0 ; --i ) {
    double lambda = f/x[i];
    water_group_velocity->SetPoint(cntr++,lambda,y[i]);
  }
  water_group_velocity->SetName("water_group_velocity");
  water_group_velocity->Write();
  */
    // Load PMT information
  char* ratroot = getenv("RATROOT");
  if (ratroot == static_cast<char*>(NULL)) {
    cerr << "Environment variable $RATROOT must be set" << endl;
    assert(ratroot);
  }
  string rat     = string(ratroot);
  string pmtfile = rat;
  pmtfile += "/data/pmt/snoman.ratdb";
  RAT::DB * db = RAT::DB::Get();
  assert(db);
  db->LoadFile(pmtfile);
  RAT::DBLinkPtr pmtInfo = db->GetLink("PMTINFO");
  assert(pmtInfo);
  PMTInfo pmt_info;
  pmt_info.x_pos = pmtInfo->GetDArray("x");
  pmt_info.y_pos = pmtInfo->GetDArray("y");
  pmt_info.z_pos = pmtInfo->GetDArray("z");
  string geofile = rat;
  geofile += "/data/geo/snoplus_water.geo";
  db->Load(geofile);
  //RAT::DB::Get()->LoadDefaults();
  db->Load(pmtfile);
  RAT::DU::Utility::Get()->BeginOfRun();
  
  for( int iEvent = 0; iEvent < tree->GetEntries() ; ++iEvent) {
    tree->GetEntry(iEvent);
    if ( !ProcessEventNtuple(rDS,ntuple,led_info,pmt_info) ) return 1; 
  }

  ntuple->Write();
  ntuple_file->Close();
  
  PhysicsNr vgroup = GroupVelocity(led_info.name );
  cout << "vgroup = " << vgroup.value << " +/- " << vgroup.error << " ns/mm (" << (vgroup.error/vgroup.value)*1000. << " permille error)" << endl;
  double test_time = 5000./vgroup.value; // time for 5 meters
  cout << "For 5 m: " << test_time << " +/- " << (vgroup.error/vgroup.value)*test_time << " ns" << endl;
  double n_ref = (c*1000./1E9)/vgroup.value;
  cout << "The effective refractive index is therefore: " << n_ref << " +/-" << (vgroup.error/vgroup.value)*n_ref << " (compare to 1.33)" << endl;
  
  benchmark.Stop("MAKENTUPLE");
  benchmark.Show("MAKENTUPLE");
  return 0;
}

