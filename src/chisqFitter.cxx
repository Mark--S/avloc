#include <include/AVLocTools.h>
#include <RAT/DB.hh>
#include <TMath.h>
#include <TFile.h>
//#include <TNTuple.h>
#include <TVector3.h>
#include <TMinuit.h>
#include <RAT/DU/Utility.hh>
#include <RAT/DU/GroupVelocity.hh>
#include <RAT/DB.hh>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <TF1.h>
#include <sstream>
#include <RAT/DU/GroupVelocity.hh>
#include <RAT/DU/LightPathCalculator.hh>
using namespace std;
int fibre_nr;
int sub_nr;
//Distance cut for the PMTS
double distCut = 1500.;
//Systematic errors from LightPathCalculator
//First term 20mm/c due to locality value 
//Second Term due to 
//PMT and LED info being loaded in from db
LEDInfo led;     
PMTInfo pmts;  
TNtuple * ntuple;
TVector3  LEDPos;
TVector3  PMTPos;
RAT::DU::LightPathCalculator lp;
RAT::DU::GroupVelocity gv;
double vg;
//Number of times pmt is hit
double * numHits;
//Average hit time for each PMT 
double * hitTimes;
//Errors on each hit hime
double * hitErrors;
TH1D ** hitHistos;
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

int main(int argc, char ** argv){
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
    pmtfile += "/data/pmt/airfill2.ratdb";
    string geofile = rat;
    geofile += "/data/geo/snoplus.geo";
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
    cout << "Obtained PMT positions "<< pmts.x_pos[1000]<<" " <<pmts.y_pos[1000]<<" "<<pmts.z_pos[1000]<< endl;
    cout << "Direction "<<pmts.x_dir[1000]<<" "<<pmts.y_dir[1000]<<"  "<<pmts.z_dir[1000]<<endl;
    //Loading up root file
    string ntuple_filename = argv[1];
    string plot_filename = argv[2];
    TFile * ntuple_file = new TFile(ntuple_filename.data(),"READ");
    if ( !ntuple_file->IsOpen() ) {
        cerr << "Could not open file " << ntuple_filename << endl;
        return 0;
    }

    TFile * plot_file = new TFile(plot_filename.data(),"RECREATE");
    if ( !plot_file->IsOpen() ) {
        cerr << "Could not open file " << plot_filename << endl;
        return 0;
    }
   // Histogram to store fit values for offset and errors
  TH1D * offsetAndErrors = new TH1D("offsetAndErrors","offsetAndErrors",100,0,100);
  ntuple = (TNtuple*)ntuple_file->Get("avloctuple");
    //Obtaining the number of fibres
    cout << "Getting fibre Numbers" << endl;
    /*TH1I * fibreHisto = new TH1I("fibreNr","fibreNr",100,0,100);
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
    }*/
    fibreNumbers.push_back(79);
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
        hitHistos = (TH1D**)malloc(sizeof(TH1D*)*numPMTS);
        for(int i=0; i<numPMTS; i++){
            numHits[i] = 0;
            hitTimes[i] = 0;
            hitErrors[i] = 0;
            //delete hitHistos[i];
            //hitHistos[i]=0;
            char name[20];
            sprintf(name,"pmt%d",i);
            hitHistos[i] = new TH1D(name,name,51,0,50);
        }
        cout << "Performing Time cuts" << endl;
        timeCuts(fibreNum);
        cout << "Completed Time cuts" << endl;
        TMinuit min(1);
        min.SetFCN(funcn);
        min.SetErrorDef(1.0);
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
        double value;
        double error;
        min.GetParameter(0,value,error);
        offsetAndErrors->SetBinContent(fibreNumbers[i],value);
        offsetAndErrors->SetBinError(fibreNumbers[i],error);
        free(numHits);
        free(hitTimes);
        free(hitErrors);
    }
    offsetAndErrors->Write();
    plot_file->Close();
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
                    //hitTimes[lcn]+=time;
                    hitHistos[lcn]->Fill(time);
                }
            }
        }
    }
    //Iterating over the hitTimes array and dividing by numhits to get average
   /* for(int i=0; i<numPMTS; i++){
        hitTimes[i]/=numHits[i];
        //If less than 30 hits cant really do statistics
        if(numHits[i]<=30){
            hitTimes[i] = 0;
            numHits[i] = 0;
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
    */
    for(int i=0; i<numPMTS; i++){
        if(numHits[i]<=30){
            numHits[i]=0;
        }
        else{
            hitHistos[i]->Fit("gaus");
            TF1 * f = hitHistos[i]->GetFunction("gaus");
            hitTimes[i]=f->GetParameter(1);
            hitErrors[i]=f->GetParError(1);
        }
    }

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
    double localityVal = 10;
    double energy = lp.WavelengthToEnergy(506.787e-6);
    //double energy = 0.00000243658;
    lp.CalcByPosition(led.position, PMTPos, energy, localityVal);
    double distInWater = lp.GetDistInWater();
    double distInScint = lp.GetDistInInnerAV();
    double distInAV = lp.GetDistInAV();
    double timeOfFlight = gv.CalcByDistance(distInScint,distInAV,distInWater,energy);
    //Adding time spent in pmt
    TVector3 pmtDir(0,0,-1.0);
    TVector3 entryDir = lp.GetIncidentVecOnPMT();
    double angleOfEntry = entryDir.Angle(pmtDir)*TMath::RadToDeg();
    timeOfFlight += gv.PMTBucketTime(angleOfEntry);
    //cout << "calcTime: "<<calcTime<< " Fibre number: "<<fibreNumber<<" PMT LCN: "<<LCN<<endl;
    return timeOfFlight;
};







