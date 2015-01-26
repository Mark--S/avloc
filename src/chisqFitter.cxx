#include <include/avlocTools.h>
#include <RAT/DB.hh>
#include <TFile.h>
#include <TNTuple.h>
#include <TVector3.h>
#include <TMinuit.h>
#include <RAT/DU/Utility.hh>
#include <RAT/DU/GroupVelocity.hh>
#include <RAT/DB.hh>
#include <iostream>

using namespace std;
int fibre_nr;
int sub_nr;
//Distance cut for the PMTS
double distCut = 1000.;
//PMT and LED info being loaded in from db
LEDInfo led;     
PMTInfo pmts;  
TNtuple * ntuple;
TVector3  LEDPos;
TVector3 * PMTPos;
double vg;
double error = (1.5*1.5 + 2.0*2.0);

//Minuit unction to minimise
void funcn(Int_t & npar, Double_t * deriv, Double_t& f, Double_t * par, Int_t flag){
    int events = ntuple->GetEntries();
    double chisq=0;
    for(int i=0; i<events; i++){
        ntuple->GetEntry(i);
        int  fibre  = (int)ntuple->GetArgs()[0];
        int    sub  = (int)ntuple->GetArgs()[1];
        int    lcn  = (int)ntuple->GetArgs()[2];
        double time = (double)ntuple->GetArgs()[3];
        double dis  = (double)ntuple->GetArgs()[4]; 
        if(dis<distCut){
            //cout << "Getting pmt vector for: "<< lcn << " x: "<<  pmts.x_pos[lcn]<< "  y: " << pmts.y_pos[lcn]<< "  z: " << pmts.z_pos[lcn]<<endl;
            PMTPos =  new TVector3(pmts.x_pos[lcn],pmts.y_pos[lcn],pmts.z_pos[lcn]);
            //cout << "Got vector" << endl;
            led = GetLEDInfoFromFibreNr(fibre,sub);
            LEDPos = led.position;
            TVector3 norm = (*PMTPos+LEDPos);
            norm.SetMag(1.0);
            double calcTime = (par[0]*norm-LEDPos).Mag()-(*PMTPos-par[0]*norm).Mag();
            calcTime/=vg;
            chisq+= (time-calcTime)*(time-calcTime)/error;
            delete PMTPos;
        }
    }
    f = chisq;
}

int main(){
    //Obtaining the fibres and sub fibre we want to fire from
    LoadDataBase("fitter.log");
    pmts = GetPMTpositions();
    cout << "Obtained PMT positions "<< pmts.x_pos[1000] << endl;
    //Loading up root file
    TFile * input = TFile::Open("summary_ntuple.root");
    ntuple = (TNtuple*)input->Get("avloctuple");
    const double time = RAT::DU::Utility::Get()->GetGroupVelocity().CalcByDistance(0.0,0.0,1.0);
    vg = 1.0/time;
    TMinuit min(1);
    min.SetFCN(funcn);
    min.SetErrorDef(0.5);
    min.SetPrintLevel(1); 
    min.DefineParameter(0,"Radius of AV",9000,50,8500,9500);
    min.Migrad();
    return 0;
}


