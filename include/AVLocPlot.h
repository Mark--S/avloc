//
// Plotting tools for AV location project
// S.J.M.Peeters@sussex.ac.uk, June 2014
//
#ifndef __AVLOCPLOT_H__
#define __AVLOCPLOT_H__

#include <cmath>
#include <string>

#include <TNtuple.h>
#include <TH2.h>
#include <TVector2.h>
#include <TVector3.h>

// Tools for plotting on a SNO+ flat map
// Originaly from Ken Clark, via James Waterfield
// From http://www.rwgrayprojects.com/rbfnotes/polyhed/PolyhedraData/Icosahedralsahedron/Icosahedralsahedron.pdf

const double a = 1.0 / 5.5;
const double b = a * sqrt( 3.0 ) / 2.0;

const TVector2 A12a = TVector2( a / 2.0, 0.0 );
const TVector2 A12b = TVector2( 3.0 * a / 2.0, 0.0 );
const TVector2 A12c = TVector2( 5.0 * a / 2.0, 0.0 );
const TVector2 A12d = TVector2( 7.0 *a / 2.0, 0.0 );
const TVector2 A12e = TVector2( 9.0 * a / 2.0, 0.0 );
const TVector2 A2a = TVector2( 0.0, b );
const TVector2 A2b = TVector2( 5.0 * a, b );
const TVector2 A17a = TVector2( a / 2.0 , 2.0 * b );
const TVector2 A17b = TVector2( 11.0 * a / 2.0 , 2.0 * b );
const TVector2 A51a = TVector2( a, 3.0 * b );
const TVector2 A51b = TVector2( 2.0 * a, 3.0 * b );
const TVector2 A51c = TVector2( 3.0 * a, 3.0 * b );
const TVector2 A51d = TVector2( 4.0 * a, 3.0 * b );
const TVector2 A51e = TVector2( 5.0 * a, 3.0 * b );
const TVector2 A27 = TVector2( 4.0 * a, b );
const TVector2 A46 = TVector2( 3.0 * a, b );
const TVector2 A31 = TVector2( 2.0 * a, b );
const TVector2 A6 = TVector2( a, b );
const TVector2 A37 = TVector2( 9.0 * a / 2.0 , 2.0 * b );
const TVector2 A33 = TVector2( 3.0 * a / 2.0 , 2.0 * b );
const TVector2 A58 = TVector2( 5.0 * a / 2.0 , 2.0 * b );
const TVector2 A54 = TVector2( 7.0 * a / 2.0 , 2.0 * b );
 
TVector2 TransformCoord( const TVector3& V1, const TVector3& V2, const TVector3& V3, const TVector2& A1, const TVector2& A2, const TVector2& A3,const TVector3& P );
TVector2 IcosProject( TVector3 pmtPos );

// use ntuple to plot flat map
TH2D * flatmap_ntuple(TNtuple * ntuple, double distance = 100000., int fibre_nr = 44, int sub_nr = 0, double time_min = 0., double time_max = 500. , bool in = true);

// use ntuple to plot the time histograms
void time_histograms(TNtuple * ntuple, double distance = 100000., int fibre_nr = 44, int sub_nr = 0);

// plot histogram of PMT distance vs recorded offset: 
// (mean time observed - mean time expected)*vgroup
// distance: max PMT distance to be considered (mm)
// fibre_nr: fibre nr to analyse 
void plot_offset(TNtuple * ntuple, double distance = 5000., int fibre_nr = 44, int sub_nr = 0);
#endif
