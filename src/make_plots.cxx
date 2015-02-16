// 
// Plot TELLIE ntuple for AV reflection fit
// S.J.M.Peeters@sussex.ac.uk, June 2014
//
#include <iostream>
#include <string>

#include <TBenchmark.h>
#include <TFile.h>
#include <TNtuple.h>

#include <RAT/DB.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/Run.hh>
#include <RAT/DU/Utility.hh>

#include "include/AVLocTools.h"
#include "include/AVLocProc.h"
#include "include/AVLocPlot.h"

using namespace std;

int main(int argc,char **argv)
{
  //TBenchmark benchmark;
  //benchmark.Start("MAKEPLOTS");
  string ntuple_filename;
  string plot_filename;
  double distance;
  int fibre_nr;
  int sub_nr;
  if ( argc != 6 ) {
    cerr << "Usage: " << argv[0] << " <ntuple filename> <output filename for plots> <distance cut (mm)> <fibre nr> <sub_nr>" << endl;
    return 1;
  }
  else {
    ntuple_filename = argv[1];
    plot_filename   = argv[2];
    distance        = atof(argv[3]);
    fibre_nr        = atoi(argv[4]);
    sub_nr          = atoi(argv[5]);
  }
  LoadDataBase("make_plots.log");
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
  
  TFile * ntuple_file = new TFile(ntuple_filename.data(),"READ");
  if ( !ntuple_file->IsOpen() ) {
    cerr << "Could not open file " << ntuple_filename << endl;
    return 0;
  }
  TNtuple * ntuple = (TNtuple*)ntuple_file->Get("avloctuple");
  assert(ntuple);
  
  TFile * plot_file = new TFile(plot_filename.data(),"CREATE");
  if ( !plot_file->IsOpen() ) {
    cerr << "Could not open file " << plot_filename << endl;
    return 0;
  }
  TH2D * hflatmap = flatmap_ntuple(ntuple,distance,fibre_nr,sub_nr,0.,500.,1);
  //time_histograms(ntuple,distance,fibre_nr,sub_nr);
  plot_offset(ntuple,distance,fibre_nr,sub_nr);
  plotAverageHitOffset(ntuple,distance);
  cout << "Made Histograms"<<endl;
  plot_file->Close();
  cout << "Closed the file"<<endl;
  
  //benchmark.Stop("MAKEPLOTS");
  //benchmark.Show("MAKEPLOTS");
  return 0;
}

