#include <include/AVLocTools.h>
#include <RAT/DB.hh>
#include <TFile.h>
//#include <TNTuple.h>
#include <TVector3.h>
#include <TMinuit.h>
#include <RAT/DU/Utility.hh>
#include <RAT/DU/GroupVelocity.hh>
#include <RAT/DB.hh>
#include <iostream>
#include <vector>
#include <sstream>
#include <RAT/DU/GroupVelocity.hh>
#include <RAT/DU/LightPathCalculator.hh>
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
RAT::DU::LightPathCalculator lp;
RAT::DU::GroupVelocity gv;
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
double trialFunction(int,int,double);
void timeCuts(int);


//Minuit unction to minimise
void funcn(Int_t & npar, Double_t * deriv, Double_t& f, Double_t * par, Int_t flag){
    double chisq=0;
    for(int i=0; i<numPMTS; i++){
        if(numHits[i]==0){
            continue;
        }

        double trial = trialFunction(fibreNum,i,par[0]);
        //cout << "hit Times Avg: "<<hitTimes[i]<<" hitErrors "<<hitErrors[i]<<" PMT Number "<<i<<endl;
        chisq+=((trial-hitTimes[i])/hitErrors[i])*((trial-hitTimes[i])/hitErrors[i]);
//        cout << "hitimes average: "<< hitTimes[i] << "Hit times errors: " << hitErrors[i] << endl;
    }
    cout << "Chisq: "<<chisq<<endl;
    f = chisq;
}

int main(){
    stringstream ss;
    //Obtaining the fibres and sub fibre we want to fire from
    LoadDataBase("fitter.log");    
    RAT::DB* db = RAT::DB::Get();
    char* ratroot = getenv("RATROOT");
    if (ratroot == static_cast<char*>(NULL)) {
        cerr << "Environment variable $RATROOT must be set" << endl;
        assert(ratroot);
    }
    string rat     = string(ratroot);
    string pmtfile = rat;
    pmtfile += "/data/pmt/snoman.ratdb";
    string geofile = rat;
    geofile += "/data/geo/snoplus_water.geo";
    db->Load(pmtfile);
    db->Load(geofile);
    std::cout << "DOING BEGIN OF RUN§" << std::endl;
    RAT::DU::Utility::Get()->BeginOfRun();
    std::cout << "DONe BEGIN OF RUN§" << std::endl;
    pmts = GetPMTpositions();
    numPMTS = pmts.x_pos.size();
    gv = RAT::DU::Utility::Get()->GetGroupVelocity();
    lp = RAT::DU::Utility::Get()->GetLightPathCalculator();
    lp.SetELLIEReflect(true);
    cout << "Obtained PMT positions "<< pmts.x_pos[1000] << endl;
    //Loading up root file
    TFile * input = TFile::Open("totalNtupleBackend.root");
    ntuple = (TNtuple*)input->Get("avloctuple");
    //Obtaining the number of fibres
    cout << "Getting fibre Numbers" << endl;
    TH1I * fibreHisto = new TH1I("fibreNr","fibreNr",100,0,100);
    int events = ntuple->GetEntries();
    for(int i=0; i<events; i++){
        ntuple->GetEntry(i);
        int  fibre  = (int)ntuple->GetArgs()[0];
        fibreHisto->Fill(fibre-1);
    }

    for(int i=1; i<=fibreHisto->GetXaxis()->GetNbins();i++){
        if(fibreHisto->GetBinContent(i)!=0){
            fibreNumbers.push_back(i);
        }
    }
    cout << "Obtained fibre numbers" << endl;
    for(unsigned int i=0; i<fibreNumbers.size(); i++){
        cout << "Fibres: "<<fibreNumbers[i]<<endl;
    }
    //THIS RETURNS NAN NEED TO FIX POSSIBLY OTHER BUGS IN CODE WHERE THIS USED AS WELL
    double time = RAT::DU::Utility::Get()->GetGroupVelocity().CalcByDistance(0.0,0.0,1.0);
    cout << "Time for light to traval 1 mm "<<time<<endl;
    vg = 1.0/time;
    cout << "Group Velocity: " << vg << endl;
    //Allocating hit time means, errors, and num Hits;
    cout<<"Starting Fits"<< endl;
    for(unsigned int i=0; i<fibreNumbers.size(); i++){
        fibreNum = fibreNumbers[i];
        numHits = (double *) malloc(sizeof(double)*numPMTS);
        hitTimes = (double *) malloc(sizeof(double)*numPMTS);
        hitErrors = (double *) malloc(sizeof(double)*numPMTS);
        for(int i=0; i<numPMTS; i++){
            numHits[i] = 0;
            hitTimes[i] = 0;
            hitErrors[i] = 0;

        }
        cout << "Performing Time cuts" << endl;
        timeCuts(fibreNum);
        cout << "Completed Time cuts" << endl;
        TMinuit min(1);
        min.SetFCN(funcn);
        min.SetErrorDef(4);
        min.SetPrintLevel(1);
        string fibreNum;
        ss << fibreNumbers[i];
        ss >> fibreNum;
        ss.clear();
        string fibreRad = "fibre "+fibreNum;
        string fibreRadOffset = fibreRad+" Offset";
        cout << fibreRad <<endl;
        min.DefineParameter(0,fibreRadOffset.c_str(),0,50,200,-200);
        min.Migrad();
        free(numHits);
        free(hitTimes);
        free(hitErrors);
    }
    return 0;
}

//Method to perform the time and distance cuts and fill up the hitTimes and hitError arrays with the data to be fitted to
void timeCuts(int fibreNumber){
    //Calculating time cut limits Using fibre FT003A and PMT LCN 2755
    //Upper time is pmt just below dist cut seperation 2086mm
    //double upperTime = trialFunction(14,6459,5500);
    //Lower time is pmt close to fibre seperation 152mm
    //double lowerTime = trialFunction(3,2329,6500);
    //double rPSUP = 8900;
    //double rAV = 6000;
    //Big number is 14.5 deg in rad
    //double phi =asin((rPSUP-rAV)*tan(0.253072742)/rAV);
    //Phi value for upper time limit assuming rav shifited 50cm away
    //double phiUpper =asin((rPSUP-rAV+500)*tan(0.253072742)/rAV);
    //double distCut = 2*sin(phi)*rPSUP;
    //Angle of reflection off AV see logbook
    //double xAngle = 0.253072742+phiUpper;
    //Minus 500 term for AV being 50cm closer as minimum path
    //double lowerTime = 2*(rPSUP-rAV-500)/vg;
    //Have to divide by 2 to get hypotenuse as dist cut double, but also have a factor of 2 for light to av and light reflected off AV factors cancel
    //double upperTime = distCut/(vg*sin(xAngle));
    double distCut = 2500;
    double lowerTime = 10;
    double upperTime = 30;
    cout << "Distance Cut is :"<< distCut << endl;
    cout << "Upper Time "<<upperTime<<" Lower Time "<<lowerTime<<endl;
    int events = ntuple->GetEntries();
    for(int i=0; i<events; i++){
        ntuple->GetEntry(i);
        int  fibre  = (int)ntuple->GetArgs()[0];
        int    sub  = (int)ntuple->GetArgs()[1];
        int    lcn  = (int)ntuple->GetArgs()[2];
        double time = (double)ntuple->GetArgs()[3];
        double dis  = (double)ntuple->GetArgs()[4]; 
        //Fibre number cut
        //cout << "Upper Time "<<upperTime<<" Lower Time "<<lowerTime<<" Actual Time"<<time<<endl;
        if(fibre == fibreNumber){
            //Distance cuts
            if(dis < distCut){
                //Time Cuts
                if(lowerTime<time && time<upperTime){
                    //cout << "Fibre: "<<fibre<<"  Fibre Number:"<<fibreNumber<<"  PMT LCN: "<<lcn<<"  Distance Cuts:  "<<dis<<endl;
                    //cout << "Got past cuts LCN: "<<lcn<<endl;
                    numHits[lcn]++;
                    hitTimes[lcn]+=time;
                }
            }
        }
    }
    //Iterating over the hitTimes array and dividing by numhits to get average
    for(int i=0; i<numPMTS; i++){
        hitTimes[i]/=numHits[i];
        //If less than 30 hits cant really do statistics
        if(numHits[i]<=30){
            hitTimes[i] = 0;
        }
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
            //Distance cuts
            if(dis < distCut){
                //Time Cuts
                if(lowerTime<time && time< upperTime){
                    hitErrors[lcn]+= (time-hitTimes[lcn])*(time-hitTimes[lcn]);
                }
            }
        }
    }
    //Dividing the hitErrors array by N-1 and sqrting to get errors
    for(int i=0; i<numPMTS; i++){
        //If less than 30 hits cannot really do statistics
        if(numHits[i]<=30){
            hitErrors[i] = 0;
            numHits[i] = 0;
            continue;
        }
        //cout << "Hit errors before  nomalization and sqrt: "<<hitErrors[i]<<endl;
        hitErrors[i]/=numHits[i]-1;
        hitErrors[i] = sqrt(hitErrors[i]);
        //cout << "Num Hits "<<numHits[i]<<endl;
     //  cout << "Hit errors after  nomalization and sqrt: "<<hitErrors[i]<<endl;
        if(hitErrors[i]==0){
            numHits[i] = 0;
        }
    };

};

double trialFunction(int fibreNumber, int LCN,double AVOffset ){
    led = GetLEDInfoFromFibreNr(fibreNumber,0);
    PMTPos =  TVector3(pmts.x_pos[LCN],pmts.y_pos[LCN],pmts.z_pos[LCN]);
    
    /*TVector3 norm = (PMTPos+LEDPos);
    norm.SetMag(1.0);
    double calcTime = ((RAV*norm)-LEDPos).Mag()+(PMTPos-(RAV*norm)).Mag();
    
    //cout <<" Distance Travelled "<<calcTime<<endl;
    calcTime = calcTime/vg;
    */
    lp.SetAVOffset(AVOffset);
    double energy = 0.00000243658;
    double localityVal = 20;
    lp.CalcByPosition(led.position, PMTPos, energy, localityVal);
    double distInWater = lp.GetDistInWater();
    double distInScint = lp.GetDistInScint();
    double distInAV = lp.GetDistInAV();
    double timeOfFlight = gv.CalcByDistance(distInScint,distInAV,distInWater,energy);
    //cout << "calcTime: "<<calcTime<< " Fibre number: "<<fibreNumber<<" PMT LCN: "<<LCN<<endl;
    return timeOfFlight;
};







