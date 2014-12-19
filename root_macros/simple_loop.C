//
// simple_loop.C : simple example of looping over SNO+ data
// S.J.M.Peeters@sussex.ac.uk - June 2014
//

// function to load the root file
void LoadRootFile(string filename, TTree **tree, RAT::DS::Root **rDS, RAT::DS::Run **rRun)
{
  TFile *file = new TFile(filename.data());
  (*tree) = (TTree*)file->Get( "T" );
  TTree *runTree = (TTree*)file->Get("runT");
  assert(runTree);
  *rDS = new RAT::DS::Root();
  (*tree)->SetBranchAddress( "ds", &(*rDS) );
  assert(rDS);
  *rRun = new RAT::DS::Run();
  assert(rRun);
  runTree->SetBranchAddress( "run", &(*rRun) );
  runTree->GetEntry();
}

// function to load database
void LoadDataBase(string logname, RAT::DB **db)
{
  RAT::Log::Init(logname.data());
  (*db) = RAT::DB::Get();
  assert(db);	
  char* glg4data = getenv("GLG4DATA");
  if (glg4data == static_cast<char*>(NULL)) {
    cerr << "ratzdab::ratdb: Environment variable $GLG4DATA must be set" << endl;
    assert(glg4data);
  }
  string data = string(glg4data);
  cout << "LoadDataBase: loading defaults..." << endl;
  (*db)->LoadDefaults();
}	

// process event
bool ProcessEvent(RAT::DS::Root * rDS)
{
  for( int iEV = 0; iEV < rDS->GetEVCount(); ++iEV) {
    RAT::DS::EV* rEV = rDS->GetEV(iEV);
    for( int ipmt = 0; ipmt < rEV->GetPMTCalCount(); ++ipmt) {
      int pmtID = rEV->GetPMTCal(ipmt)->GetID();
      Double_t PMTTime = rEV->GetPMTCal(ipmt)->GetTime();
      cerr << pmtID << ": " << PMTTime << endl;
    }
  }
  return true;
}

// main function
void simple_loop(string filename = "50Z_FT079B.root") {
  
  RAT::DS::Root * rDS  = NULL;
  RAT::DS::Run  * rRun = NULL;
  TTree         * tree = NULL;
  LoadRootFile(filename,&tree,&rDS,&rRun);
  RAT::DB       * db   = NULL;
  LoadDataBase("simple_loop.log",&db);
  
  for( int iEvent = 0; iEvent < 10 /*tree->GetEntries()*/ ; ++iEvent) {
    tree->GetEntry(iEvent);
    ProcessEvent(rDS);
  }
  return;
} 
