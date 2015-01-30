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

#include "include/AVLocTools.h"
#include "include/AVLocProc.h"
#include "include/AVLocPlot.h"

using namespace std;

int main(int argc,char **argv)
{
  TBenchmark benchmark;
  benchmark.Start("MAKEPLOTS");
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
  hflatmap->Write();
  plot_file->Close();
  
  benchmark.Stop("MAKEPLOTS");
  benchmark.Show("MAKEPLOTS");
  return 0;
}

