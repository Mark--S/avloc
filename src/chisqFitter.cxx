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
#include <vector>

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
TVector3  PMTPos;
double vg;
double error = (1.5*1.5 + 2.0*2.0);
//Number of times pmt is hit
double * numHits;
//Average hit time for each PMT 
double * hitTimes;
//Errors on each hit hime
double * hitErrors;
//Number of pmts
int numPMTS;
//fibreNumbers being fired from
vector<double> fibreNumbers;
//fibreNumber iterator
int fibreNum;
double trialFunction(int,int ,double);
void timeCuts(int);


//Minuit unction to minimise
void funcn(Int_t & npar, Double_t * deriv, Double_t& f, Double_t * par, Int_t flag){
    double chisq=0;
    for(int i=0; i<numPMTS; i++){
        double trial = trialFunction(fibreNum,i,par[0]);
        //cout << "hit Times Avg: "<<hitTimes[i]<<" hitErrors "<<hitErrors[i]<<" PMT Number "<<i<<endl;
        chisq+=((trial-hitTimes[i])/hitErrors[i])*((trial-hitTimes[i])/hitErrors[i]);
    }
    cout << "Chisq: "<<chisq<<endl;
    f = chisq;
}

int main(){
    //Obtaining the fibres and sub fibre we want to fire from
    LoadDataBase("fitter.log");    
    RAT::DB* db = RAT::DB::Get();
    db->Load( "/Users/markstringer/Documents/PHD/rat/snoing/install/rat-5.0.0/data/geo/snoplus_water.geo" );
    db->Load("/Users/markstringer/Documents/PHD/rat/snoing/install/rat-5.0.0/data/pmt/snoman.ratdb");
    std::cout << "DOING BEGIN OF RUN§" << std::endl;
    RAT::DU::Utility::Get()->BeginOfRun();
    std::cout << "DONe BEGIN OF RUN§" << std::endl;
    pmts = GetPMTpositions();
    numPMTS = pmts.x_pos.size();
    cout << "Obtained PMT positions "<< pmts.x_pos[1000] << endl;
    //Loading up root file
    TFile * input = TFile::Open("summary_ntupleOffsetWorking.root");
    ntuple = (TNtuple*)input->Get("avloctuple");
    //Obtaining the number of fibres
    cout << "Getting fibre Numbers" << endl;
    TH1D * fibreHisto = new TH1D("fibreNr","fibreNr",100,0,100);
    int events = ntuple->GetEntries();
    for(int i=0; i<events; i++){
        ntuple->GetEntry(i);
        int  fibre  = (int)ntuple->GetArgs()[0];
        fibreHisto->Fill(fibre);
    }

    for(int i=1; i<=fibreHisto->GetXaxis()->GetNbins();i++){
        if(fibreHisto->GetBinContent(i)!=0){
            fibreNumbers.push_back(i);
        }
    }
    cout << "Obtained fibre numbers" << endl;

    //THIS RETURNS NAN NEED TO FIX POSSIBLY OTHER BUGS IN CODE WHERE THIS USED AS WELL
    double time = RAT::DU::Utility::Get()->GetGroupVelocity().CalcByDistance(0.0,0.0,1.0);
    cout << "Time for light to traval 1 mm "<<time<<endl;
    vg = 1.0/time;
    cout << "Group Velocity: " << vg << endl;
    //Allocating hit time means, errors, and num Hits;
    cout<<"Starting Fits"<< endl;
    for(int i=0; i<(int)fibreNumbers.size(); i++){
        fibreNum = fibreNumbers[i];
        numHits = new double[numPMTS];
        hitTimes = new double[numPMTS];
        hitErrors = new double[numPMTS];
        timeCuts(fibreNum);
        TMinuit min(1);
        min.SetFCN(funcn);
        min.SetErrorDef(0.5);
        min.SetPrintLevel(1);
        string fibreRad = "Radius of AV using fibre "+std::to_string(fibreNumbers[i]);
        min.DefineParameter(0,fibreRad.c_str(),6000,50,5500,6500);
        min.Migrad();
        delete[] numHits;
        delete[] hitTimes;
        delete[] hitErrors;
    }
    return 0;
}

//Method to perform the time and distance cuts and fill up the hitTimes and hitError arrays with the data to be fitted to
void timeCuts(int fibreNumber){
    //Calculating time cut limits Using fibre FT003A and PMT LCN 2755
    double upperTime = trialFunction(4,3944,5500);
    double lowerTime = trialFunction(4,3944,6500);
    int events = ntuple->GetEntries();
    for(int i=0; i<events; i++){
        ntuple->GetEntry(i);
        int  fibre  = (int)ntuple->GetArgs()[0];
        int    sub  = (int)ntuple->GetArgs()[1];
        int    lcn  = (int)ntuple->GetArgs()[2];
        double time = (double)ntuple->GetArgs()[3];
        double dis  = (double)ntuple->GetArgs()[4]; 
        //Fibre number cut
        cout << "Upper Time "<<upperTime<<" Lower Time "<<lowerTime<<" Actual Time"<<time<<endl;
        if(fibre == fibreNumber){
            //Time Cuts
            if(time>lowerTime && time<upperTime){
                cout << "Got past cuts"<<endl;
                numHits[lcn]++;
                hitTimes[lcn]+=time;
            }
        }
    }
    //Iterating over the hitTimes array and dividing by numhits to get average
    for(int i=0; i<numPMTS; i++){
        hitTimes[i]/=numHits[i];
        //cout << "Hit Times: "<<hitTimes[i]<<" num Hits "<<numHits[i]<<endl;
    }
    //Now iterating over again to get the errors 
    for(int i=0; i<events; i++){
        ntuple->GetEntry(i);
        int  fibre  = (int)ntuple->GetArgs()[0];
        int    sub  = (int)ntuple->GetArgs()[1];
        int    lcn  = (int)ntuple->GetArgs()[2];
        double time = (double)ntuple->GetArgs()[3];
        double dis  = (double)ntuple->GetArgs()[4]; 
        //Fibre Number cut
        if(fibre==fibreNumber){
            //Time Cuts
            if(time>lowerTime && time< upperTime){
                hitErrors[lcn]+= (time-hitTimes[lcn])*(time-hitTimes[lcn]);
            }
        }
    }
    //Dividing the hitErrors array by N-1 and sqrting to get errors
    for(int i=0; i<numPMTS; i++){
        hitErrors[i]/=i-1;
        hitErrors[i] = sqrt(hitErrors[i]);
    };

};

double trialFunction(int fibreNumber, int LCN,double RAV){
    led = GetLEDInfoFromFibreNr(fibreNumber,0);
    PMTPos =  TVector3(pmts.x_pos[LCN],pmts.y_pos[LCN],pmts.z_pos[LCN]);
    LEDPos = led.position;
    TVector3 norm = (PMTPos+LEDPos);
    norm.SetMag(1.0);
    double calcTime = ((RAV*norm)-LEDPos).Mag()+(PMTPos-(RAV*norm)).Mag();
    //cout <<" Distance Travelled "<<calcTime<<endl;
    calcTime = calcTime/vg;
    //cout << "calcTime: "<<calcTime<< " Fibre number: "<<fibreNumber<<" PMT LCN: "<<LCN<<endl;
    return calcTime;
};








