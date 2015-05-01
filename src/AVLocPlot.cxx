//
// Plotting tools for AV location project
// S.J.M.Peeters@sussex.ac.uk, June 2014
//
#include "include/AVLocPlot.h"
#include "include/AVLocTools.h"

#include <RAT/DS/Entry.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DB.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DU/GroupVelocity.hh>
#include <RAT/DU/LightPathCalculator.hh>

#include <TF1.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1D.h>
#include <TCanvas.h>

#include <map>

using namespace std;

TVector2 TransformCoord( const TVector3& V1, const TVector3& V2, const TVector3& V3, const TVector2& A1, const TVector2& A2, const TVector2& A3,const TVector3& P ) {
    TVector3 xV = V2 - V1;
    TVector3 yV = ( ( V3 - V1 ) + ( V3 - V2 ) ) * 0.5;
    TVector3 zV = xV.Cross( yV ).Unit();

    double planeD = V1.Dot( zV );

    double t = planeD / P.Dot( zV );

    TVector3 localP = t*P - V1;

    TVector2 xA = A2 - A1;
    TVector2 yA = ( ( A3 - A1 ) +( A3 - A2 ) ) * 0.5;

    double convUnits = xA.Mod() / xV.Mag();

    TVector2 result;
    result = localP.Dot( xV.Unit() ) * xA.Unit() * convUnits;
    result += localP.Dot( yV.Unit() ) * yA.Unit() * convUnits + A1;
    return result;
}

TVector2 IcosProject( TVector3 pmtPos ){
    TVector3 pointOnSphere( pmtPos.X(), pmtPos.Y(), pmtPos.Z() );
    pointOnSphere = pointOnSphere.Unit();
    pointOnSphere.RotateX( -45.0 );
    // From http://www.rwgrayprojects.com/rbfnotes/polyhed/PolyhedraData/Icosahedralsahedron/Icosahedralsahedron.pdf                                                                                                                                                                                               
    const double t = ( 1.0 + sqrt( 5.0 ) ) / 2.0;
    const TVector3 V2 = TVector3( t * t, 0.0, t * t * t ).Unit();
    const TVector3 V6 = TVector3( -t * t, 0.0, t * t * t ).Unit();
    const TVector3 V12 = TVector3( 0.0, t * t * t, t * t ).Unit();
    const TVector3 V17 = TVector3( 0.0, -t * t * t, t * t ).Unit();
    const TVector3 V27 = TVector3( t * t * t, t * t, 0.0 ).Unit();
    const TVector3 V31 = TVector3( -t * t * t, t * t, 0.0 ).Unit();
    const TVector3 V33 = TVector3( -t * t * t, -t * t, 0.0 ).Unit();
    const TVector3 V37 = TVector3( t * t * t, -t * t, 0.0 ).Unit();
    const TVector3 V46 = TVector3( 0.0, t * t * t, -t * t ).Unit();
    const TVector3 V51 = TVector3( 0.0, -t * t * t, -t * t ).Unit();
    const TVector3 V54 = TVector3( t * t, 0.0, -t * t * t ).Unit();
    const TVector3 V58 = TVector3( -t * t, 0.0, -t * t * t ).Unit();
    // Faces {{ 2, 6, 17}, { 2, 12, 6}, { 2, 17, 37}, { 2, 37, 27}, { 2, 27, 12}, {37, 54, 27},                                                                                                                                                                                                                    
    // {27, 54, 46}, {27, 46, 12}, {12, 46, 31}, {12, 31, 6}, { 6, 31, 33}, { 6, 33, 17},                                                                                                                                                                                                                          
    // {17, 33, 51}, {17, 51, 37}, {37, 51, 54}, {58, 54, 51}, {58, 46, 54}, {58, 31, 46},                                                                                                                                                                                                                         
    // {58, 33, 31}, {58, 51, 33}}                                                                                                                                                                                                                                                                                 
    vector<TVector3> IcosahedralCentres;
    IcosahedralCentres.push_back( ( V2 + V6 + V17 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V2 + V12 + V6 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V2 + V17 + V37 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V2 + V37 + V27 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V2 + V27 + V12 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V37 + V54 + V27 ) * ( 1.0 / 3.0 ) );

    IcosahedralCentres.push_back( ( V27 + V54 + V46 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V27 + V46 + V12 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V12 + V46 + V31 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V12 + V31 + V6 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V6 + V31 + V33 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V6 + V33 + V17 ) * ( 1.0 / 3.0 ) );

    IcosahedralCentres.push_back( ( V17 + V33 + V51 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V17 + V51 + V37 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V37 + V51 + V54 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V58 + V54 + V51 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V58 + V46 + V54 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V58 + V31 + V46 ) * ( 1.0 / 3.0 ) );

    IcosahedralCentres.push_back( ( V58 + V33 + V31 ) * ( 1.0 / 3.0 ) );
    IcosahedralCentres.push_back( ( V58 + V51 + V33 ) * ( 1.0 / 3.0 ) );

    vector<double> distFromCentre;
    unsigned int uLoop;
    for( uLoop = 0; uLoop < IcosahedralCentres.size(); uLoop++ ){
        distFromCentre.push_back( ( IcosahedralCentres[uLoop] - pointOnSphere ).Mag() );
    }
    const int face = min_element( distFromCentre.begin(), distFromCentre.end() ) - distFromCentre.begin() + 1;
    TVector2 resultPosition;
    switch(face)
    {
        case 1://{ 2, 6, 17}
            resultPosition = TransformCoord( V2, V6, V17, A2a, A6, A17a, pointOnSphere );
            break;
        case 2://{ 2, 12, 6}                        
            resultPosition = TransformCoord( V2, V12, V6, A2a, A12a, A6, pointOnSphere );
            break;
        case 3://{ 2, 17, 37}
            resultPosition = TransformCoord( V2, V17, V37, A2b, A17b, A37, pointOnSphere );
            break;
        case 4://{ 2, 37, 27}
            resultPosition = TransformCoord( V2, V37, V27, A2b, A37, A27, pointOnSphere );
            break;
        case 5://{ 2, 27, 12}
            resultPosition = TransformCoord( V2, V27, V12, A2b, A27, A12e, pointOnSphere );
            break;
        case 6://{37, 54, 27}
            resultPosition = TransformCoord( V37, V54, V27, A37, A54, A27, pointOnSphere );
            break;
        case 7://{27, 54, 46}
            resultPosition = TransformCoord( V27, V54, V46, A27, A54, A46, pointOnSphere );
            break;
        case 8://{27, 46, 12}
            resultPosition = TransformCoord( V27, V46, V12, A27, A46, A12d, pointOnSphere );
            break;
        case 9://{12, 46, 31}
            resultPosition = TransformCoord( V12, V46, V31, A12c, A46, A31, pointOnSphere );
            break;
        case 10://{12, 31, 6}
            resultPosition = TransformCoord( V12, V31, V6, A12b, A31, A6, pointOnSphere );
            break;
        case 11://{ 6, 31, 33}
            resultPosition = TransformCoord( V6, V31, V33, A6, A31, A33, pointOnSphere );
            break;
        case 12://{ 6, 33, 17}
            resultPosition = TransformCoord( V6, V33, V17, A6, A33, A17a, pointOnSphere );
            break;
        case 13://{17, 33, 51}
            resultPosition = TransformCoord( V17, V33, V51, A17a, A33, A51a, pointOnSphere );
            break;
        case 14://{17, 51, 37}
            resultPosition = TransformCoord( V17, V51, V37, A17b, A51e, A37, pointOnSphere );
            break;
        case 15://{37, 51, 54}
            resultPosition = TransformCoord( V37, V51, V54, A37, A51d, A54, pointOnSphere );
            break;
        case 16://{58, 54, 51}
            resultPosition = TransformCoord( V58, V54, V51, A58, A54, A51c, pointOnSphere );
            break;
        case 17://{58, 46, 54}
            resultPosition = TransformCoord( V58, V46, V54, A58, A46, A54, pointOnSphere );
            break;
        case 18://{58, 31, 46}
            resultPosition = TransformCoord( V58, V31, V46, A58, A31, A46, pointOnSphere );
            break;
        case 19://{58, 33, 31}
            resultPosition = TransformCoord( V58, V33, V31, A58, A33, A31, pointOnSphere );
            break;
        case 20://{58, 51, 33}
            resultPosition = TransformCoord( V58, V51, V33, A58, A51b, A33, pointOnSphere );
            break;
    }
    return TVector2( resultPosition.X(), 2.0 * resultPosition.Y() );
}


// use ntuple to plot flat map
TH2D * flatmap_ntuple(TNtuple * ntuple, double distance, int fibre_nr, int sub_nr, double time_min, double time_max, bool in)
{
    // PMT info
    RAT::Log::Init("/dev/null");
    RAT::DB* db = RAT::DB::Get();
    assert(db);
    char* ratroot = getenv("RATROOT");
    if (ratroot == static_cast<char*>(NULL)) {
        cerr << "Environment variable $RATROOT must be set" << endl;
        assert(ratroot);
    }
    string rat     = string(ratroot);
    string pmtfile = rat;
    pmtfile += "/data/pmt/airfill2.ratdb";
    db->LoadFile(pmtfile);
    RAT::DBLinkPtr pmtInfo = db->GetLink("PMTINFO");
    assert(pmtInfo);
    vector<double> xPos = pmtInfo->GetDArray("x");
    vector<double> yPos = pmtInfo->GetDArray("y");
    vector<double> zPos = pmtInfo->GetDArray("z");

    // Loop over ntuple
    Double_t pmt_hits[10000] = {0};

    unsigned int nev = ntuple->GetEntries();
    for (unsigned int i=0 ; i < nev ; ++i) {
        ntuple->GetEntry(i);
        int  fibre  = (int)ntuple->GetArgs()[0];
        int    sub  = (int)ntuple->GetArgs()[1];
        int    lcn  = (int)ntuple->GetArgs()[2];
        double time = (double)ntuple->GetArgs()[3];
        double dis  = (double)ntuple->GetArgs()[4]; 
        bool in_distance;
        if ( in ) {
            dis < distance ? in_distance = true : in_distance = false;
        }
        else {
            dis > distance ? in_distance = true : in_distance = false;

        }
        if (in_distance && fibre == fibre_nr && sub == sub_nr && time >= time_min && time < time_max) {
            pmt_hits[lcn] += 1;
        }
    }

    // Make plot
    const int xbins = 300;
    const int ybins = 300;
    TH2D *hflatmap = new TH2D("hflatmap","SNO+ flatmap",xbins, 0 , 1, ybins, 0 , 1);
    for (unsigned int i = 0; i < 10000; i++){
        TVector3 pmtPos(xPos[i],yPos[i],zPos[i]);
        TVector2 icosProj = IcosProject(pmtPos);
        int xbin = int((1-icosProj.X())*xbins);
        int ybin = int((1-icosProj.Y())*ybins);
        int bin = hflatmap->GetBin(xbin,ybin);
        if(pmt_hits[i]>0){
            printf("Hits on PMT %d : %d\n",i,pmt_hits[i]);
            hflatmap->SetBinContent(bin,pmt_hits[i]);
        }
    }
    return hflatmap;
}


// use ntuple to plot the time histograms
void time_histograms(TNtuple * ntuple, double distance, int fibre_nr, int sub_nr)
{
    // Loop over ntuple
    TH1I * histo_map[10000];
    for (unsigned int i = 0 ; i < 10000 ; ++i ) histo_map[i] = NULL;
    unsigned int nev = ntuple->GetEntries();
    cout << "entries: " << nev << endl;
    for (unsigned int i = 0 ; i < nev ; ++i) {
        ntuple->GetEntry(i);
        double dist = (double)ntuple->GetArgs()[4]; 
        if ( dist < distance ) {
            int    lcn  = (int)ntuple->GetArgs()[2];
            double time = (double)ntuple->GetArgs()[3]; 
            int  fibre  = (int)ntuple->GetArgs()[0];
            int  sub    = (int)ntuple->GetArgs()[1];
            if ( histo_map[lcn] == NULL ) {
                char name[128];
                sprintf(name,"pmt%i",lcn);
                histo_map[lcn] = new TH1I(name,name,51,-0.5,50.5);
                histo_map[lcn]->SetXTitle("time (ns)");
            }
            //cout << "Filling Histogram "<<endl;
            //This is the line causing the bug time offset not like the old stuff
            if ( fibre == fibre_nr && sub == sub_nr && time > 0. && time < 50. ) {
                histo_map[lcn]->Fill(time);
            }
        }
    }

    // save histograms to file (needs to open!)
    TH1D * time_summary = new TH1D("time_summary","average hit time for each PMT",
            10001,-0.5,10000.5);
    TH1D * time_histo = new TH1D("time_histo","time distribution for reflections",
            501,-0.5,50.5);
    time_summary->SetXTitle("LCN");
    time_summary->SetYTitle("hit_time (ns)");
    for (unsigned int i = 0 ; i < 10000 ; ++i ) {
        if (histo_map[i] != NULL ) {
            // if at least 30 entries, calculate mean and rms
            //cout << "histo map " << i << " entries "<<histo_map[i]->GetEntries()<<endl;
            if (histo_map[i]->GetEntries() > 30 ) {
                histo_map[i]->Fit("gaus");
                histo_map[i]->Write();
                TF1 * f = histo_map[i]->GetFunction("gaus");
                assert(f);
                double mu = f->GetParameter(1);
                double si = f->GetParError(1);
                time_summary->SetBinContent(i+1,mu);
                time_summary->SetBinError  (i+1,si);
                //cout << "Filling histogram" << endl;
                time_histo->Fill(mu,1./(si*si));
            }	
        }
    }
    time_summary->Write();
    time_histo->Fit("gaus");
    time_histo->Write();
}

void plot_offset(TNtuple * ntuple, double distance, int fibre_nr, int sub_nr)
{
    int nBins = 200;
    LEDInfo   led    = GetLEDInfoFromFibreNr(fibre_nr, sub_nr);
    PMTInfo pmt_info = GetPMTpositions();
    RAT::DU::GroupVelocity  gv = RAT::DU::Utility::Get()->GetGroupVelocity();
    RAT::DU::LightPathCalculator lp = RAT::DU::Utility::Get()->GetLightPathCalculator();
    lp.SetELLIEReflect(true);
    lp.SetAVOffset(0);
    // effective refractive index:
    // need to get this from the database but is in data now ... hardcoded, i.e. improve!!
    cout << "Set up Light Path Calculator"<<endl;
    PhysicsNr n_h2o; 
    n_h2o.value = 1.3637; 
    n_h2o.error = 0.0021;
    // Loop over ntuple
    TH1I * histo_map[10000];
    TH1I * histo_mapPE[10000];
    TH1I * histo_mapNotOffset[10000];
    TH1D * distanceInAV[10000];
    TH1D * distanceInWater[10000];
    TH1D * distanceInScint[10000];
    for (unsigned int i = 0 ; i < 10000 ; ++i ) histo_map[i] = NULL;
    for (unsigned int i = 0 ; i < 10000 ; ++i ) histo_mapPE[i] = NULL;
    for (unsigned int i = 0 ; i < 10000 ; ++i ) histo_mapNotOffset[i] = NULL;
    for (unsigned int i = 0 ; i < 10000 ; ++i ) distanceInAV[i] = NULL;
    for (unsigned int i = 0 ; i < 10000 ; ++i ) distanceInWater[i] = NULL;
    for (unsigned int i = 0 ; i < 10000 ; ++i ) distanceInScint[i] = NULL;
    unsigned int nev = ntuple->GetEntries();
    for (unsigned int i = 0 ; i < nev ; ++i) {
        ntuple->GetEntry(i);
        double dist = (double)ntuple->GetArgs()[4]; 
        if ( dist < distance ) {
            //cout << "GOt Entry "<<i<<endl;
            int    lcn    = (int)ntuple->GetArgs()[2];
            double time   = (double)ntuple->GetArgs()[3]; 
            int    fibre  = (int)ntuple->GetArgs()[0];
            double peTime = (double)ntuple->GetArgs()[5];
            //cout << "Got Args"<<endl;
            if ( histo_map[lcn] == NULL ){
                char name[128];
                char namePE[128];
                char nameNotOffset[128];
                char nameAV[128];
                char nameWater[128];
                char nameScint[128];
                sprintf(name,"pmt%i",lcn);
                sprintf(nameNotOffset,"pmt no offset%i",lcn);
                sprintf(namePE,"pmtPE%i",lcn);
                sprintf(nameAV,"dist in AV %i",lcn);
                sprintf(nameWater,"dist in Water %i",lcn);
                sprintf(nameScint,"dist in Scint  %i",lcn);
                histo_map[lcn] = new TH1I(name,name,51,-25.5,25.5);
                histo_mapPE[lcn] = new TH1I(namePE,namePE,51,-25.5,25.5);
                histo_mapNotOffset[lcn] = new TH1I(nameNotOffset,nameNotOffset,51,0,50);
                distanceInAV[lcn] = new TH1D(nameAV,nameAV,50,-20,200);
                distanceInWater[lcn] = new TH1D(nameWater,nameWater,500,20,8000);
                distanceInScint[lcn] = new TH1D(nameScint,nameScint,50,-20,200);
                histo_map[lcn]->SetXTitle("time (ns)");
                histo_mapPE[lcn]->SetXTitle("time (ns)");
                histo_mapNotOffset[lcn]->SetXTitle("time (ns)");
                distanceInAV[lcn]->SetXTitle("distance (mm)");
                distanceInScint[lcn]->SetXTitle("distance (mm)");
                distanceInWater[lcn]->SetXTitle("distance (mm)");
            }
           // cout << "Set UP histo map" <<endl;
            if ( fibre == fibre_nr && time > 0. && time < 50. ) {
                TVector3 PMT_pos(pmt_info.x_pos[lcn],pmt_info.y_pos[lcn],pmt_info.z_pos[lcn]);
                //PhysicsNr tof = TimeOfFlight(led.position, PMT_pos, n_h2o, 1.);
                double localityVal = 1.0;
                double energy = lp.WavelengthToEnergy(506.787);
                lp.CalcByPosition(led.position, PMT_pos, energy, localityVal);
                double distInWater = lp.GetDistInWater();
                double distInScint = lp.GetDistInInnerAV();
                double distInAV = lp.GetDistInAV();
                double timeOfFlight = gv.CalcByDistance(distInScint,distInAV,distInWater,energy);
                histo_mapPE[lcn]->Fill(time-peTime);
                histo_map[lcn]->Fill(time-timeOfFlight);
                //cout << time << endl;
                histo_mapNotOffset[lcn]->Fill(time);
                distanceInAV[lcn]->Fill(distInAV);
                distanceInWater[lcn]->Fill(distInWater);
                distanceInScint[lcn]->Fill(distInScint);
            }
        }
    }

    // save histograms to file (needs to open!)
    TH1D * time_summary = new TH1D("time_summary","average hit time for each PMT",
            10001,-0.5,10000.5);
    TH1D * time_summary_offset  = new TH1D("time_summary_fucntionOfDistance","average hit time offset with Distance",
            nBins,0.0,distance);
    TH1D * time_summary_Distance  = new TH1D("time_summary_DistanceNoOffset","average hit time with Distance",
            nBins,0.0,distance);
    char title[128];
    sprintf(title,"time distribution for reflections, fibre %i-%i",fibre_nr,sub_nr);
    TH1D * time_histo = new TH1D("time_histo",title, nBins+1,-10.05,10.05);
    TH1D * time_histo_PE = new TH1D("time_histo_PE",title, nBins+1,-10.05,10.05);
    time_summary->SetXTitle("LCN");
    time_summary->SetYTitle("hit_time (ns)");
    time_summary_offset->SetXTitle("Distance (mm)");
    time_summary_offset->SetYTitle("Hit Time (ns)");
    time_summary_Distance->SetXTitle("Distance (mm)");
    time_summary_Distance->SetYTitle("Hit Time (ns)");
    for (unsigned int i = 0 ; i < 10000 ; ++i ) {
        if (histo_map[i] != NULL ) {
            // if at least 30 entries, calculate mean and rms
            TVector3 PMT_pos(pmt_info.x_pos[i],pmt_info.y_pos[i],pmt_info.z_pos[i]);
            double dist = (PMT_pos-led.position).Mag();
            printf("Hits on PMT timing %d : %d\n",i,histo_map[i]->GetEntries());
            if (histo_map[i]->GetEntries() > 30 && dist<distance) {
                histo_map[i]->Fit("gaus");
                histo_mapPE[i]->Fit("gaus");
                histo_map[i]->Write();
                histo_mapNotOffset[i]->Fit("gaus");
                TF1 * f = histo_map[i]->GetFunction("gaus");
                TF1 * fPE = histo_mapPE[i]->GetFunction("gaus");
                TF1 * fNotOffset = histo_mapNotOffset[i]->GetFunction("gaus");
                assert(f);
                double mu = f->GetParameter(1);
                double si = f->GetParError(1);
                double muPE = fPE->GetParameter(1);
                double siPE = fPE->GetParError(1);
                double muNotOffset = fNotOffset->GetParameter(1);
                double siNotOffset = fNotOffset->GetParError(1);
                //converting distance to bins 2500mm dist cut and 100 bins
                int binNumber = (int )(nBins*dist/distance);
                printf("PMT Number: %d Hit offset: %f distance:%f  bin Number %d\n",i,mu,dist,binNumber);
                cout << "mu: "<<muNotOffset<<" si: "<<siNotOffset<<endl;
                time_summary_offset->SetBinContent(binNumber,mu);
                time_summary_offset->SetBinError(binNumber,si);
                time_summary->SetBinContent(i+1,mu);
                time_summary->SetBinError  (i+1,si);
                time_summary_Distance->SetBinContent(binNumber,muNotOffset);
                time_summary_Distance->SetBinError(binNumber,siNotOffset);
                time_histo->Fill(mu,1./(si*si));
                time_histo_PE->Fill(muPE,1./(siPE*siPE));
            }	
        }
    }
    cout << "Writing out offset histograms" << endl;
    time_summary->Write();
    time_summary_offset->Write();
    time_summary_Distance->Write();
    time_histo->Fit("gaus");
    time_histo_PE->Fit("gaus");
    time_histo->SetXTitle("ns");
    time_histo_PE->SetXTitle("ns");
    time_histo->Write();
    time_histo_PE->Write();

}

//Method to iterate over all the fibres in the NTuple get hit histograms for distance bins and fit these to guassians drawing a histogram with mean offset and error
void plotAverageHitOffset(TNtuple * ntuple, double distance){
    RAT::DU::Utility::Get()->BeginOfRun();
    unsigned int nBins = 100;
    TH1D * distanceMap[nBins];
    PMTInfo pmt_info = GetPMTpositions();
    //cout << "Getting Group Velocity"<<endl;
    RAT::DU::GroupVelocity  gv = RAT::DU::Utility::Get()->GetGroupVelocity();
    //cout << "Getting LightPath Calculator"<<endl;
    RAT::DU::LightPathCalculator lp = RAT::DU::Utility::Get()->GetLightPathCalculator();
    lp.SetELLIEReflect(true);
    cout << "Set up light path calculator"<<endl;
    // effective refractive index:
    // need to get this from the database but is in data now ... hardcoded, i.e. improve!!
    TH1D * time_histo = new TH1D("time_histo_AllPMTS","time_histo_AllPMTS", 101,-10.05,10.05);
    time_histo->SetXTitle("Offset (ns)");
    // Loop over ntuple
    cout << "Setting distance maps up" << endl;
    for (unsigned int i = 0 ; i < nBins ; ++i ) distanceMap[i] = 0;
    unsigned int nev = ntuple->GetEntries();
    //cout << nev << endl;
    for (unsigned int i = 0 ; i < nev ; ++i){
        //if (i%1000 == 0 ) cout << "Got entry: "<<i<<endl;
        ntuple->GetEntry(i);
        double dist = (double)ntuple->GetArgs()[4]; 
        int    lcn    = (int)ntuple->GetArgs()[2];
        double time   = (double)ntuple->GetArgs()[3]; 
        if(time > 0. && time < 50. && dist < distance  ){
            int    fibre  = (int)ntuple->GetArgs()[0];
            int sub = (int) ntuple->GetArgs()[1];
            //Getting bin number from distance 
            int binNum = (int)((dist/distance)*nBins);
            if ( distanceMap[binNum] == NULL ){
                //cout << "Setting up histo"<<endl;
                char name[128];
                char nameAV[128];
                char nameWater[128];
                char nameScint[128];
                sprintf(name,"binNum %d",binNum);
                distanceMap[binNum] = new TH1D(name,name,50,-10,10);
                distanceMap[binNum]->SetXTitle("offset (ns)");
            }
            //cout << "Getting pmt Position"<<endl;
            TVector3 PMT_pos(pmt_info.x_pos[lcn],pmt_info.y_pos[lcn],pmt_info.z_pos[lcn]);
            //cout << "Getting fibbre information" <<endl;
            LEDInfo   led    = GetLEDInfoFromFibreNr(fibre, sub);
            //PhysicsNr tof = TimeOfFlight(led.position, PMT_pos, n_h2o, 1.);
            TVector3 hypPMTPos(0,0,0);
            double lambda = 508;
            double localityVal = 1.0;
            double energy = 0.00000243658;
            //cout << "Calculating by distance"<<endl;
            lp.CalcByPosition(led.position, PMT_pos, energy, localityVal);
            double distInWater = lp.GetDistInWater();
            double distInScint = lp.GetDistInInnerAV();
            double distInAV = lp.GetDistInAV();
            //cout << "Calculating Time of Flight bin Num: "<<binNum<<endl;
            double timeOfFlight = gv.CalcByDistance(distInScint,distInAV,distInWater,energy);
            distanceMap[binNum]->Fill(time-timeOfFlight);
        }
    }
    cout << "Finished getting entries now doing average histogram"<<endl;
    TH1D * time_summary_offset_Average = new TH1D("hitTimeAsFunctionOfDistanceforAllFibres","average hit time offset with Distance all Fibres",nBins,0.0,distance);
    time_summary_offset_Average->SetXTitle("Distance from fibre (mm)");
    time_summary_offset_Average->SetYTitle("Hit Time (ns)");
    cout << "Set up histogram"<<endl;
    //time_summary_offset_Average->Print("ALL");
    for (unsigned int i = 0 ; i < nBins ; ++i ) {
        if (distanceMap[i] != NULL ) {
            // if at least 30 entries, calculate mean and rms
            if (distanceMap[i]->GetEntries() > 30) {
                cout << "Fitting histogram"<<endl;
                distanceMap[i]->Fit("gaus");
                TF1 * f = distanceMap[i]->GetFunction("gaus");
                assert(f);
                double mu = f->GetParameter(1);
                double si = f->GetParameter(2);
                cout << "Setting bin error and values "<<i<<"   mu   "<<mu<<" si  "<<si<<endl;
                time_summary_offset_Average->SetBinContent(i+1,mu);
                time_summary_offset_Average->SetBinError(i+1,si);
                time_histo->Fill(mu,1.0/(si*si));

            }	
        }
    }
    //time_summary_offset_Average->Print("ALL");
    cout << "Writing out histogram"<<endl;
    //time_summary_offset_Average->SetDirectory(gDirectory->pwd());
    time_histo->Fit("gaus");
    time_histo->Write();
    time_summary_offset_Average->Write();
}





void plotAVFlightDifference(TNtuple * ntuple1, TNtuple * ntuple2 ,double distance, int fibre_nr, int sub_nr)
{
    int nBins = 200;
    LEDInfo   led    = GetLEDInfoFromFibreNr(fibre_nr, sub_nr);
    PMTInfo pmt_info = GetPMTpositions();
    // effective refractive index:
    // need to get this from the database but is in data now ... hardcoded, i.e. improve!!
    PhysicsNr n_h2o; 
    n_h2o.value = 1.3637; 
    n_h2o.error = 0.0021;
    // Loop over ntuple
    //Stores Difference of Ntuples
    TH1I * histo_map[10000];
    //Stores the values of hit times from both ntuples
    TH1I * histo_map1[10000];
    TH1I * histo_map2[10000];
    for (unsigned int i = 0 ; i < 10000 ; ++i ) histo_map1[i] = NULL;
    for (unsigned int i = 0 ; i < 10000 ; ++i ) histo_map2[i] = NULL;
    unsigned int nev = ntuple1->GetEntries();
    unsigned int nev2 = ntuple2->GetEntries();
    //Cutting off excess entries
    if(nev2<nev) nev=nev2;
    for (unsigned int i = 0 ; i < nev ; ++i) {
        ntuple1->GetEntry(i);
        double dist = (double)ntuple1->GetArgs()[4]; 
        double dist2 = (double)ntuple2->GetArgs()[4]; 
        if ( dist < distance ) {
            //cout << "GOt Entry "<<i<<endl;
            int    lcn    = (int)ntuple1->GetArgs()[2];
            double time   = (double)ntuple1->GetArgs()[3]; 
            int    fibre  = (int)ntuple1->GetArgs()[0];
            double peTime = (double)ntuple1->GetArgs()[5];
            //cout << "Got Args"<<endl;
            if ( histo_map1[lcn] == NULL ){
                char name[128];
                sprintf(name,"Newpmt%i",lcn);
                histo_map1[lcn] = new TH1I(name,name,51,0,50);
                histo_map1[lcn]->SetXTitle("time (ns)");
            }
           // cout << "Set UP histo map" <<endl;
            if ( fibre == fibre_nr && time > 0. && time < 50. ) {
                histo_map1[lcn]->Fill(time);
            }
        }
    }

    nev = ntuple2->GetEntries();
    for (unsigned int i = 0 ; i < nev ; ++i) {
        ntuple2->GetEntry(i);
        double dist = (double)ntuple2->GetArgs()[4]; 
        if ( dist < distance ) {
            //cout << "GOt Entry "<<i<<endl;
            int    lcn    = (int)ntuple2->GetArgs()[2];
            double time   = (double)ntuple2->GetArgs()[3]; 
            int    fibre  = (int)ntuple2->GetArgs()[0];
            double peTime = (double)ntuple2->GetArgs()[5];
            //cout << "Got Args"<<endl;
            if ( histo_map2[lcn] == NULL ){
                char name[128];
                sprintf(name,"2ndpmt%i",lcn);
                histo_map2[lcn] = new TH1I(name,name,51,0,50);
                histo_map2[lcn]->SetXTitle("time (ns)");
            }
           // cout << "Set UP histo map" <<endl;
            if ( fibre == fibre_nr && time > 0. && time < 50. ) {
                histo_map2[lcn]->Fill(time);
            }
        }
    }
    cout << "Got the ntuples"<<endl;
    // save histograms to file (needs to open!)
    char title[128];
    TH1D * time_summary_offset  = new TH1D("time_summaryDifferenceBetween unshifted and shifted","average hit time offset with Distance",
            nBins,0.0,distance);
    time_summary_offset->SetXTitle("Distance (mm)");
    time_summary_offset->SetYTitle("Hit Time No Offset - Hit Time offset (ns)");
    //Subtracting the histograms
    for(unsigned int i=0; i<10000; i++){ 
        if(histo_map1[i]==NULL || histo_map2[i]==NULL){
            histo_map1[i]=NULL;
            histo_map2[i]=NULL;
            continue;
        }
        histo_map1[i]->Add(histo_map2[i],-1);
    }
    cout << "Added the PMT Histograms"<<endl;
    for (unsigned int i = 0 ; i < 10000 ; ++i ) {
        if (histo_map[i] != NULL ) {
            cout << "Not NULL"<<endl;
            // if at least 30 entries, calculate mean and rms
            if (histo_map1[i]->GetEntries() > 30) {
                cout << "Fitting histogram for pmt; "<<i << endl;
                histo_map1[i]->Fit("gaus");
                cout << "Fitted histogram for pmt; "<<i << endl;
                TF1 * f = histo_map1[i]->GetFunction("gaus");
                assert(f);
                double mu = f->GetParameter(1);
                double si = f->GetParameter(2);
                TVector3 PMT_pos(pmt_info.x_pos[i],pmt_info.y_pos[i],pmt_info.z_pos[i]);
                double dist = (PMT_pos-led.position).Mag();
                //converting distance to bins 2500mm dist cut and 100 bins
                int binNumber = (int )(nBins*dist/distance);
                printf("bin Number %d\n",binNumber);
                time_summary_offset->SetBinContent(binNumber,mu);
                time_summary_offset->SetBinError(binNumber,si);
            }	
        }
    }
    cout << "Writing out offset histograms now" << endl;
    time_summary_offset->Write();
}

