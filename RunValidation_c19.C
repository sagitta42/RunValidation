/* Any questions, complaints and suggestions: sssagitta@gmail.com */

/********************************************
 * 
 * ~~~~~ TIPS FOR SEARCHING THROUGH THE CODE ~~~~~
 * "Structure": HELPER FUNCTIONS, BxValidation functions, MAIN FUNCTION
 *
 *  Canvases and histograms:
 *  	Precalibration: CANVAS NO 0
 *  	Calibration: CANVAS NO 1, CANVAS NO 2, CANVAS NO 3
 *  	Run validation: canvas 00, canvas 01, ..., canvas 16
 * 			search forward --> initialization of the histograms for this canvas in the class
 * 			search forward x 2 --> histogram assignment (so, names, binning etc.)
 * 			search backward --> drawing the canvas
 * 
 *  Type of events:
 *  	LASER, PULSER, RANDOM, NEUTRINO, INTERNAL MUON, EXTERNAL MUON, NEUTRON
 * 		also TT8, TT32, TT64, TT1, TT1&BTB0, TT1&BTB1, TT2, TT128
 * 		also case 8, case 32, case 64, case 1, case 2, case 128 (switch on trigger type inside the for loop on events)
 *
********************************************/

#include <iostream>
#include <vector>
#include "TPaveLabel.h"
#include "TNetFile.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH2F.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TMath.h"
#include "TROOT.h"
#include <cstdlib>
#include <cstdio>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TTreeResult.h>
#include <TSQLRow.h>
#include <string>
#include <sstream>
#include <time.h>
#include "TPaveText.h"

#include <fstream>
#include <math.h>
#include "TImage.h"
#include "TPaveStats.h"

//#include "event/BxEvent.hh"
#include "./event/BxEvent.hh"

#include "/home/production/RunValidation/RunStability_c19.C"
//#include "RunStability_c19_mylocal.C"
//STAB

using namespace std;

////// HELPER FUNCTIONS //////
long RunEventsBxdb(int RunNumber);
string StartDate(int RunNumber);
string TimeNow();
long int StartDate_int(int RunNumber);
string GetPath (int runnum, int cycle);
string GetWeek (int runnum);
string GetYear(int runnum);
void SetGeneralStyle();
float VerticalMean(TH1F* h, float rms_limit, bool IgnoreZeros);
double SkewGaus(double *x, double *p);
bool FileCheck(string filepath);

// the class with all the parameters, histos and canvases that will be accessed by all the other functions 
class BxValidation{
	/*
	 * This class contains all the histograms that are filled and canvases that are drawn by the Run Validation function and later saved by Save Validation
	 */
	
	public:

		BxValidation(){run_int=0; clean_cin=false;};
		BxValidation(string filepath, int Run, int Cycle);
		~BxValidation(){};

		// weird quasi destructors
		void CleanShelf();
		void Clean100();

		// weird quasi constructor (can be used only after StartEndInfo)
		void InitializeHistos();

		// getters
		int GetRun();
		int GetNumEvents();

		// plotters and other
		int ElectronicsCanvases();
		void StartEndInfo(int endpoint);
		void FillHistos();
		void DrawCanvases();
		void OuterDetectorCheck();
		void ElectronicsRescale();
		bool SaveValidation();

	private:
		// Echidna file
		string filepath;
		TFile* frun;
		// info about the run number and cycle
		string runnum, cycle, week, year;
		Int_t run_int;
		Int_t cycle_int;
		// tree
		TTree* t;
		BxEvent* ev;
		
		// info extracted from the run
		int num_events, start_evnum, end_evnum,  RunEvents, t_diff, btb_thresh;
		double start_time, end_time;
		double previous_neu_time;
		
		// statistics from the histos
		Float_t VerticalMean1, VerticalMean2, VerticalMean3, VerticalMean4, VerticalMean5;  //laben
		Float_t MuonRateID, MuonRateOD, MuonRateOD25; //muon
		Int_t   LastMcrEvent, LastAlignedEvent;  //muon
		Int_t good_laben_channels, good_muon_channels, good_fadc_channels;
	
		// other stuff to be saved in the info file
		string loc_validation_time, loc_valid, loc_group, loc_root_files, loc_start_time;


		//// stuff used when filling histos ////
		
		// counting TT2 events to check OD behaviour later
		int TT2count;
		
		// counting pulser and TT1 events to rescale electronics calibration histograms later
		int TT32count;
		int TT1count;
		
		// counting TT1 muons (TT1&BTB4) to later draw sigma lines
		int muonTT1count;
		int TT128count;
		
		// big bin size for the histos with 200 bins (canvas 02 and canvas 10)
		int binstep;

		bool clean_cin;
		
		//// precalibration and calibration ////
	
		// canvas NO 0: precalibration
		TCanvas* c_precalib;
		TH2F* charge_vs_channel;
		TH2F* time_vs_channel; 
		TH2F* time_vs_channel_good;
		
		// canvases NO 1, 2 and 3: calibration
//		TCanvas *raw_lg;
		TCanvas *trigger, *laser;
		TCanvas* raw_lg_renorm;
		// histos
		TH1F* trigref;
		TH1F* lasref;
		TH1F* raw_lg_pulser; // canvas no 1
		TH1F* raw_lg_neutrino;
		TH1D* trgr[16]; // canvas no 2
		TH1D* las[16]; // canvas no 3
	
		//// run validation ////
		
		// canvases 00 to 17
		TCanvas* c[19];

		// histos
	
		// note: event size = decoded hits, pulse shape = decoded hits start times

		// ***** INNER DETECTOR *****//

		// event size VS evnum: canvas 01 
		TH1F* IDevsize[6];
		int pulserIDylim[2]; // limit on the pulser y-axis in canvases 01 and 02
		// event size VS evnum in 200 bins: canvas 02 
		TH1F* IDevsize200[6];
		// neutrino trigger rate VS time: canvas 03
		// max time: 10 hours, 1-minute bins
		TH1F* IDtrgRate[2];	
		// hits in neutrino trigger: canvas 04, canvas 07
		TH1F* IDneutrg[9];
		TH1F* IDc7sample[4];
		// pulse shape: canvas 05
		TH1F* IDtimes[6];
		// event size VS channels: canvas 06
		TH1F* IDchannels[4];
		// hits in muon triggers: canvas 12
		TH1F* IDhitsMuon[5];
		// neutron trigger: canvas 16, canvas 17
		TH1F* IDtrg128[11];

	
		// ***** OUTER DETECTOR *****//
	
		// OD flag stability: canvas 08
		TH1F* ODmuonFlag[6];
		// event size VS evnum: canvas 09
		TH1F* ODevsize[6];
		int pulserODylim[2]; // limit on the pulser y-axis in canvases 01 and 02
		// event size VS evnum in 200 bins: canvas 10 
		TH1F* ODevsize200[6];
		// muon trigger rate VS time: canvas 11
		// max time: 10 hours, 10-minute bins
		TH1F* ODtrgRate[4];
		// hits in muon triggers: canvas 13
		TH1F* ODhitsTT1[4];
		TH1F* ODhitsTT2[4];
		// pulse shape: canvas 14
		TH1F* ODtimes[6];
		// event size VS channels: canvas 15
		TH1F* ODchannels[6];


};



RunStability *run_stability; 
// STAB

 
// Has to be a pointer declared outside to have access from both RunValidation() and SaveValidation()
// when root prompt opened the first time: this line creates a new object and initializes counters to 0. if somebody runs SV without RV it will prevent it.
BxValidation* bx = new BxValidation();


//////////////////////////////////////////////////////////
/// **************** MAIN FUNCTION **************** ///
//////////////////////////////////////////////////////////

void RunValidation(string filepath){
	/*
	 * INPUT: path to the file e.g. "Run029114_c18.root" or "/bxstorage/rootfiles/cycle_18/2017/Jul_16/Run029114_c18.root"
	 * OUTPUT: mmm... A LOT of canvases :)
	 */


	// ----------------------- stability.lock ---------------------- //

	ifstream temp("stability.lock");

	if(temp){
		cout << endl << "!! WARNING: stability.lock exists !!" << endl << endl;
		bool ask = true;
		string answer;
		
		while(ask){
			cout << "Are you sure nobody else is running RunValidation now? Do you want to continue RunValidation (yes/no)? ";
			cin >> answer;
			if(answer == "yes"){ask = false;}
			else if(answer == "no"){temp.close(); return;}
		}
	}
	temp.close();


	// ----------------------- do checks and initialize ---------------------- //

	// ********************* check if the file exists
	if(!FileCheck(filepath)){return;}

	// if somebody does RV+RV, we need to clean the shelf first --> kind of a lame destructor that cleans the memory but doesn't delete the pointer that was declared outside of this function
	if(bx->GetRun()){bx->CleanShelf();}
	// now that we have freed the memory in case somebody ran RV+RV in a row, we can declare a new one
	cout << endl << "~~~ All set up for a new Run Validation ~~~" << endl;
	
	// common canvas style
	SetGeneralStyle();
	
	
	/// find the run number, cycle number and week from the root file name ///
	int length = filepath.size();
	int run_int, cycle_int;

	// read and check the format of the filename
	if( sscanf( filepath.substr(length-18,18).c_str(), "Run%6d_c%2d.root", &run_int, &cycle_int) != 2){ 
		cerr << "The file format should be Run0XXXXX_cXX.root, otherwise cannot extract run number and cycle. Sorry." << endl;
		return;
	}
	
	// the output of the StartDate_int is meaningless to a human, but it is used to determine the Week (e.g. "Jul_16" later)
	if( StartDate_int(run_int) < 0 ){
		cerr << "Start date is not valid. Later when I study the code written in 2008 better, I can also explain what to do in this case and why it's not valid." << endl;
		return;
	}
	
	// If everything is correct, we can proceed
	bx = new BxValidation(filepath, run_int, cycle_int);
	
	// ----------------------- ELECTRONICS CALIBRATION & PRECALIBRATION ---------------------- //

	cout << "......proceeding to the calibration and precalibration......" << endl << endl;
	
	// if it returns 1 means something was wrong
	if( bx->ElectronicsCanvases() ){ return; }
	

	// ----------------------- RUN VALIDATION ---------------------- //

	cout << "......proceeding to the validation of this run......" << endl << endl;
	cout << "#######################################" << endl;
	cout << "Total number of events: " << bx->GetNumEvents() << endl;
	cout << "Enter the last event number to be validated (0 = use all the events, -100 = quit): ";
	
	/// ********************* user input for event number to validate
	int endpoint;
	bool flag = true;
	
	while(flag){
		cin >> endpoint;
		
		// if it wasn't int
		while(cin.fail()){
		cout << "Please enter an integer: ";
		cin.clear(); // clean the fail flag
		cin.ignore(256,'\n'); // cleans the input
		cin >> endpoint;
		}

		// -100 means exit
		if(endpoint == -100){
			// if somebody does RV + -100 + RV need to clean
			bx->Clean100();
			bx = new BxValidation();
			return;
		}

		if(endpoint < 0){	cout << "A positive integer...: "; }
		
		else if( (endpoint < 20000) && (endpoint != 0)){
			cout << "The event number has to be above 20000. Enter another one: ";
		}

		else if (endpoint > bx->GetNumEvents()){
			cout << "The number you entered is larger than the total number of events (" << bx->GetNumEvents() << "). Enter another one: ";
		}
		else{flag = false;}
	}
	
	// ********************* main steps

	// find out info about start/end evnum/time before initializing histograms
	
#ifdef debug
	gROOT->SetBatch(kTRUE);
# endif

	bx->StartEndInfo(endpoint);
	bx->InitializeHistos();
	bx->FillHistos();
	bx->DrawCanvases();
	bx->OuterDetectorCheck(); // check for amazing OD behaviour 
	bx->ElectronicsRescale();

	cout << "#######################################" << endl << endl;
	cout << "If the run is acceptable, save it by running SaveValidation()" << endl << endl;
	cout << "#######################################" << endl << endl;

  run_stability = new RunStability (bx->run_int);
	// STAB
//	bx->clean_cin = true; // not needed?...
  return; 
} // end of RunValidation macro




void SaveValidation(){

	bool success = bx->SaveValidation();

	// if SaveValidation was successful
	if(success){
		// in order to do RV+SV+RV we need everything cleaned and initialized. or in case RV+SV+SV happens to prevent it.
		bx->CleanShelf();
		bx = new BxValidation();
		cout << endl << "~~~ All set up for a new Run Validation ~~~" << endl << endl;
	}
	// if it wasn't, don't remove the old info yet, what if they change their mind

}



////////////////////////////////////////////////////////////////
/// **************** BxValidation functions **************** ///
////////////////////////////////////////////////////////////////



// --------------- constructor --------------- //
BxValidation::BxValidation(string path, int Run, int Cycle){

	filepath = path;
	run_int = Run;
	cycle_int = Cycle;
	
	int length = filepath.size();
	
	// string with all the zeros e.g. Run009114_c18.root -> "009114"
	runnum = filepath.substr(length - 15, 6); 
	// string e.g. "18"
	cycle = filepath.substr(length - 7, 2);
	

	// If the date is valid, we can get the week out of it
	week = GetWeek(run_int); // string
	year = GetYear(run_int); // string
	
	cout << endl << "#######################################" << endl;
	cout << "File: " << filepath  << endl;
	cout << "Run number: " << runnum << endl;
	cout << "Cycle: " << cycle << endl;
	cout << "Week: " << week << endl;
	cout << "Year: " << year << endl;
	cout << "#######################################" << endl << endl;

	/// read the tree ///	
	frun = TFile::Open(filepath.c_str());
	frun->GetObject("bxtree",t);
	num_events = t->GetEntries();
	ev = new BxEvent();
	t->SetBranchAddress("events",&ev);
	
	// important to initialize counters to 0 not to have garbage value
	TT2count = 0;
	TT32count = 0;
	TT1count = 0;
	muonTT1count = 0;
	TT128count = 0;
	// this always stays zero, not implemented, but need to avoid garbage value
	good_fadc_channels = 0;

	clean_cin = true;

}

// ------------------------------------------------------------

// --------------- weird quasi destructors --------------- //


void BxValidation::CleanShelf(){
	// delete all the histos and pointers
	    delete frun;
			delete c_precalib;
//			delete raw_lg; // we're not plotting it now
			delete raw_lg_renorm;
			delete trigger;
			delete laser;
			
			for(int i = 1; i < 18; i++){ delete c[i];}

			for(int i = 0; i < 6; i++){
				delete IDevsize[i];
				delete IDevsize200[i];
				delete IDtimes[i];
				delete ODmuonFlag[i];
				delete ODevsize[i];
				delete ODevsize200[i];
				delete ODtimes[i];
				delete ODchannels[i];
			}
		
			for(int i = 0; i < 2; i++){delete IDtrgRate[i];}
			
			for(int i = 0; i < 9; i++){delete IDneutrg[i];}
			
			for(int i = 0; i < 4; i++){
				delete IDchannels[i];
				delete ODtrgRate[i];
				delete ODhitsTT1[i];
				delete ODhitsTT2[i];
			}
			
			for(int i = 0; i < 5; i++){delete IDhitsMuon[i];}
			
			for(int i = 0; i < 11; i++){delete IDtrg128[i];}
}

// ------------------------------------------------------------

// destructor only for the case if -100 was the input
void BxValidation::Clean100(){
	delete c_precalib;
	delete trigger;
	delete laser;
}

// ------------------------------------------------------------

// --------------- getters --------------- //
int BxValidation::GetRun(){ return run_int; }
int BxValidation::GetNumEvents(){ return num_events; }

// ------------------------------------------------------------


// --------------- plotters --------------- //
int BxValidation::ElectronicsCanvases(){
	/*
	 * Read and draw the electronics calibrations and precalibrations related canvases
	 */

//	SetGeneralStyle();

	// the common path to calibrations and precalibrations
	string ancillary = "http://bxmaster-data.lngs.infn.it//bxstorage/rootfiles/cycle_" + cycle;
	
	// if the run from January this year belongs to a week from December last year

	if( (week.substr(0,3) == "Dec") && (StartDate(run_int).substr(5,2) == "01") ){
		int yr;
		sscanf(year.c_str(),"%4d",&yr);
		// go back to last year
		yr--;
		char temp[4];
		sprintf(temp,"%d",yr);
		year = temp;
	}

	ancillary += "/" + year + "/" + week + "/ancillary/Run" + runnum;

	/// calibration ///
	string el_calib_file = ancillary + "_electronics_calibrations_c" + cycle + ".root";
	
	TFile* fcalib = TFile::Open(el_calib_file.c_str());
	
	// check whether file exists
	if(!fcalib ){
		cerr << endl << "File " << el_calib_file << " does not exist." << endl;
		return 1;
	}
	// check whether file is Zombie
	else if( fcalib->IsZombie() ){
		cerr << endl << "File " << el_calib_file << " is a zombie." << endl;
		return 1;
	}

	/// precalibration ///
	string el_precalib_file = ancillary + "_precalibrations_c" + cycle + ".root";

	TFile* fprecalib = TFile::Open(el_precalib_file.c_str());
	
	// check whether file exists
	if(!fprecalib ){
		cerr << endl << "File " << el_precalib_file << " does not exist." << endl;
		return 1;
	}
	// check whether file is Zombie
	else if( fprecalib->IsZombie() ){
		cerr << endl << "File " << el_precalib_file << " is a zombie." << endl;
		return 1;
	}
		
	/// if everything is fine with the file, we can initialize histos and proceed ///
	cout << "#######################################" << endl;
	cout << "Electronics calibration file:" << endl;
	cout << el_calib_file << endl;
	cout << endl;
	cout << "Electronics precalibration file: " << endl;
	cout << el_precalib_file << endl;
	cout << "#######################################" << endl;
	cout << endl;
	
	/// ************************* DRAW ************************* ///

	
	// el. precalib.: canvas no 0
	charge_vs_channel = (TH2F*)fprecalib->Get("barn/bx_calib_laben_decoding/charge_vs_channel");
	time_vs_channel = (TH2F*)fprecalib->Get("barn/bx_calib_laben_decoding/time_vs_channel");
	time_vs_channel_good = (TH2F*)fprecalib->Get("barn/bx_calib_laben_decoding/time_vs_channel_good");

	// el. calib.
	// these are not plotted but are used to get lg channel numbers
	trigref = (TH1F*)fcalib->Get("barn/bx_calib_laben_electronics/trigref");
	lasref = (TH1F*)fcalib->Get("barn/bx_calib_laben_electronics/lasref");

	// these are plotted on canvases no 1 and 2
	TH2F *pulser_charge_vs_lg = (TH2F*)fcalib->Get("barn/bx_calib_laben_electronics/pulser_charge_vs_lg");
	TH2F *laser_charge_vs_lg = (TH2F*)fcalib->Get("barn/bx_calib_laben_electronics/laser_charge_vs_lg");

	// these are not plotted now, but plotted later on canvas no 3 after normalization
	raw_lg_pulser = (TH1F*)fcalib->Get("barn/bx_calib_laben_electronics/raw_Nhits_pulser_vs_lg");
	raw_lg_neutrino = (TH1F*)fcalib->Get("barn/bx_calib_laben_electronics/raw_Nhits_neutrino_vs_lg");

	// variables to recycle
	int lgnum;
	string title; // canvas titles
	TPaveLabel *labels[16]; // labels
	char name[13]; // name of the histo projection 

	int ww = 900;
	int hh = 600;

	// ************ CANVAS NO 0: precalibration
	title = "[0] Run" + runnum + " Precalibration";
	c_precalib = new TCanvas("c_precalib", title.c_str(), ww, hh);
	c_precalib->Divide(2,2);

	// pad 1
	c_precalib->cd(1);
	TH1* frame = c_precalib->GetPad(1)->DrawFrame(0, 0, charge_vs_channel->GetXaxis()->GetXmax(), charge_vs_channel->GetYaxis()->GetXmax());
	frame->SetTitle("Raw charge [ADC channel]");
	frame->SetXTitle(charge_vs_channel->GetXaxis()->GetTitle());
	frame->SetYTitle(charge_vs_channel->GetYaxis()->GetTitle());
	frame->SetTitleSize(0.055);
	frame->SetTitleSize(0.045,"X");
	charge_vs_channel->Draw("sames");
	
	// pad 2
	c_precalib->cd(2);
	c_precalib->GetPad(2)->SetLogy();
	
	frame = c_precalib->GetPad(2)->DrawFrame(-10,0.1, charge_vs_channel->GetYaxis()->GetXmax(), (charge_vs_channel->ProjectionY())->GetMaximum()*2);
	frame->SetTitle("Raw charge (y-projection)");
	frame->SetXTitle(charge_vs_channel->GetYaxis()->GetTitle());
	frame->SetTitleSize(0.055);
	frame->SetTitleSize(0.045,"X");
	frame->Draw();

	(charge_vs_channel->ProjectionY())->Draw("sames");

	c_precalib->cd(3);
	time_vs_channel->GetYaxis()->SetTitle("Hit time after precalibrations [ns]");
	time_vs_channel->SetTitle("Time alignment of all channels");
	time_vs_channel->SetTitleSize(0.055);
	time_vs_channel->SetTitleSize(0.045,"X");
	time_vs_channel->SetTitleSize(0.045,"Y");
	time_vs_channel->Draw();
	
	
	c_precalib->cd(4);
	time_vs_channel_good->GetYaxis()->SetTitle("Hit time after precalibrations [ns]");
	time_vs_channel_good->SetTitleSize(0.045,"Y");
	(time_vs_channel_good->ProjectionY())->SetTitle("Correctly precalibrated channels (y-projection)");
	(time_vs_channel_good->ProjectionY())->SetTitleSize(0.055);
	(time_vs_channel_good->ProjectionY())->Draw();
	c_precalib->GetPad(4)->SetLogy();
	
	c_precalib->Modified(); c_precalib->Update();


	// ************ CANVAS NO 0.5 :D not plotted :D
	//raw hits non normalizzati vanno visualizzati qua solo nel caso production = 0
	//altrimenti histo si normalizza e plotta later (dopo contare Nevents)
//	title = "[1] Run" + bx->runnum + " lg distribution of raw hits";
//	bx->raw_lg = new TCanvas( "raw_lg",title.c_str() );
//	bx->raw_lg->Divide(1,2);
//	bx->raw_lg->cd(1); 
//	bx->raw_lg_pulser->SetTitleSize(0.055);
//	bx->raw_lg_pulser->SetTitleSize(0.045,"X");
//	bx->raw_lg_pulser->Draw();
//	
//	bx->raw_lg->cd(2);
//	bx->raw_lg_neutrino->SetTitleSize(0.055);
//	bx->raw_lg_neutrino->SetTitleSize(0.045,"X");
//  bx->raw_lg_neutrino->Draw(); 
//
////	bx->raw_lg->Modified();
//	bx->raw_lg->Update(); // otherwise the canvas is blank

	// ************ CANVAS NO 1
	// reycle the variable title
	title = "[1] Run" + runnum + " Trigger reference channels, pulser trigger";
	trigger = new TCanvas("trigger", title.c_str(), ww, hh);
	
	trigger->Divide(4,4);    
	
	for(int i = 1; i < 17; i ++){
		lgnum = (int) trigref->GetBinContent(i);
		
		if(lgnum) {
			sprintf(name, "trgref%d", lgnum);
			trgr[i] = pulser_charge_vs_lg->ProjectionY(name, lgnum, lgnum);
			trigger->cd(i);
			trgr[i]->SetAxisRange(0,100);
			trgr[i]->GetXaxis()->SetTitle("Raw charge [ADC channels]");
			trgr[i]->SetTitleSize(0.055);
			trgr[i]->SetTitleSize(0.050,"X");
			trgr[i]->Draw();
			labels[i] = new TPaveLabel(0.15, 0.6, 0.45, 0.8, name, "NDC");
			labels[i]->SetBorderSize(0);
			labels[i]->SetFillColor(kWhite);
			labels[i]->Draw();
			trigger->Modified();
		}		
	}
	
	trigger->Update();

	// ************ CANVAS NO 2
	// reycle the variable title
	title = "[2] Run" + runnum + " Laser reference channels, laser trigger";
	laser = new TCanvas("laser", title.c_str(), ww,hh);
	laser->Divide(3,3);    
	
	for(int i = 1; i < 17; i ++){
		lgnum = (int) lasref->GetBinContent(i);
			
		if(lgnum){
			sprintf(name, "lasref%d", lgnum);
			las[i] = laser_charge_vs_lg->ProjectionY(name,lgnum, lgnum);
			laser->cd(i);
			las[i]->SetAxisRange(0,100);
			las[i]->GetXaxis()->SetTitle("Raw charge [ADC channels]");
			las[i]->SetTitleSize(0.055);
			las[i]->SetTitleSize(0.050,"X");
			las[i]->Draw();
			labels[i] = new TPaveLabel(0.15, 0.6, 0.45, 0.8, name, "NDC");
			labels[i]->SetBorderSize(0);
			labels[i]->SetFillColor(kWhite);
			labels[i]->Draw();
			laser->Modified();
		}
	}

	laser->Update();

	return 0;
}

// ------------------------------------------------------------

void BxValidation::StartEndInfo(int endpoint){
	
	// if zero given by user, take the last event in the run
	if(endpoint == 0){ end_evnum = num_events; } 	
	// otherwise truncate at given endpoint
	else{ end_evnum = endpoint; }

	// end event time
// oink
	//	BxEvent* ev = new BxEvent();
//	t->SetBranchAddres("events",&ev);
	t->GetEntry(end_evnum - 1);
	end_time = ev->GetTrigger().GetGpsTimeSec();
	
	// some info taken from the last event and saved for later
	btb_thresh = ev->GetTrigger().GetBtbThreshold(); // int
	good_laben_channels = ev->GetLaben().GetNLivePmts(); // int
	good_muon_channels = ev->GetMuon().GetNPmts(); // int

	// error: 'class BxEvent' has no member named 'GetFadc'
	//	bx->good_fadc_channels = ev->GetFadc().GetNRawChannels(); // int
	
	// start time of the first "ID neutrino" (i.e. TT1)
	// CURRENTLY first neytrino after 1000 events
	bool neuflag = true; int evno = 1000;
//	double previous_neu_time; // will need for canvas 4

	while(neuflag){
		t->GetEntry(evno);
		// check if neutrino
		if(ev->GetTrigger().GetTrgType() == 1){
			neuflag = false;
			start_time= ev->GetTrigger().GetGpsTimeSec();
			previous_neu_time = ev->GetTrigger().GetGpsTimeNs();
			start_evnum = ev->GetEvNum();
		}
		// if not, keep going
		evno++;
	}

	// save for later
	t_diff = (int) (end_time - start_time);

	// end time of the last "ID neutrino" (i.e. TT1)
	// NOT USED??
	float last_neu_time = 0;
	neuflag = true; evno = end_evnum;

	while(neuflag){
		t->GetEntry(evno);
		// check if neutrino
		if(ev->GetTrigger().GetTrgType() == 1){
			neuflag = false;
			last_neu_time = ev->GetTrigger().GetGpsTimeSec();
			// initial value is of the first neutrino
		}
		// if not, keep going
		evno--;
	}

}

// ------------------------------------------------------------

void BxValidation::InitializeHistos(){
	/*
	 * Is called from Fill Histos to initialize the histos for the BX_VALIDATION object (bx) which is global.
	 * To avoid cluttering the main function
	 */
	

	// note: event size = decoded hits, pulse shape = decoded hits start times

	// ***** INNER DETECTOR *****//


	// event size VS evnum: canvas 01 
	IDevsize[0] = new TH1F("IDlaserEvSize", "Laser (TT8)", end_evnum - start_evnum, start_evnum, end_evnum);
	IDevsize[1] = new TH1F("IDpulserEvSize", "Pulser (TT32)", end_evnum - start_evnum, start_evnum, end_evnum);
	IDevsize[2] = new TH1F("IDrandomEvSize", "Random (TT64)", end_evnum - start_evnum, start_evnum, end_evnum);
	IDevsize[3] = new TH1F("IDneutrinoEvSize", "Neutrino (TT1&BTB0)", end_evnum - start_evnum, start_evnum, end_evnum);
	IDevsize[4] = new TH1F("IDmuonTT1EvSize", "Internal Muon (TT1&BTB4)", end_evnum - start_evnum, start_evnum, end_evnum);
	IDevsize[5] = new TH1F("IDmuonTT2EvSize", "External Muon (TT2)", end_evnum - start_evnum, start_evnum, end_evnum);
	pulserIDylim[0] = 100000; // y-min
	pulserIDylim[1] = -10; // y-max

	// event size VS evnum in 200 bins: canvas 02 
	IDevsize200[0] = new TH1F("IDlaserEvSize200", "Laser (TT8)", 200, 0, 200);
	IDevsize200[1] = new TH1F("IDpulserEvSize200", "Pulser (TT32)", 200, 0, 200);
	IDevsize200[2] = new TH1F("IDrandomEvSize200", "Random (TT64)", 200, 0, 200);
	IDevsize200[3] = new TH1F("IDneutrinoEvSize200", "Neutrino (TT1&BTB0)", 200, 0, 200);
	IDevsize200[4] = new TH1F("IDmuonTT1EvSize200", "Internal Muon (TT1&BTB4)", 200, 0, 200);
	IDevsize200[5] = new TH1F("IDmuonTT2EvSize200", "External Muon (TT2)", 200, 0, 200);

	// neutrino trigger rate VS time: canvas 03
	// max time: 10 hours, 1-minute bins
	IDtrgRate[0] = new TH1F("IDtrgRateAll", "Internal events (TT1) all", 3600*10/60, 0, 3600*10);
	IDtrgRate[0]->Sumw2();
	IDtrgRate[1] = new TH1F("IDtrgRate100", "Internal events (TT1) ev. size > 100", 3600*10/60, 0, 3600*10);
	IDtrgRate[1]->Sumw2();
	
	// hits in neutrino trigger: canvas 04, canvas 07
	// canvas 07
	IDneutrg[0] = new TH1F("IDneutrinoRaw", "Raw hits in neutrino events (TT1&BTB0)",1000,0,1000);
	IDneutrg[1] = new TH1F("IDneutrinoDecoded", "Decoded hits in neutrino events (TT1&BTB0)",400,0,400);
	IDneutrg[2] = new TH1F("IDneutrinoClustered1", "Clustered hits in 1-cluster events",400,0,400);
	IDneutrg[3] = new TH1F("IDneutrinoCharge","Charge in 1-cluster events",400,0,400);
	// extra in canvas 04
	IDneutrg[4] = new TH1F("IDneutrinoClusters", "Clusters in neutrino events (TT1&BTB0)",4,0,4);
	IDneutrg[5] = new TH1F("IDtimeDifference", "Time difference between neutrino events",400,0,400); // in ms
	IDneutrg[6] = new TH1F("IDclusterTime","Cluster start time distribution", 4500, -18500, -18500+4500);
	IDneutrg[7] = new TH1F("IDneutrinoDecoded1", "Decoded hits in 1-cluster events",85,0,85);
	IDneutrg[8] = new TH1F("IDtimeDifferenceZoom", "Time difference on a smaller scale",(int)(50/0.2),0,50); // in mus

	// sample run for comparison with canvas 07
	TFile* can7 = TFile::Open("/home/production/run_plots/29114.root");
	// the historical file has different histo names as following

	IDc7sample[0] = (TH1F*)can7->Get("n_raw_hits");
	IDc7sample[0]->SetTitle(IDneutrg[0]->GetTitle());
	
	IDc7sample[1] = (TH1F*)can7->Get("n_dec_hits");
	IDc7sample[1]->SetTitle(IDneutrg[1]->GetTitle());
	
	IDc7sample[2] = (TH1F*)can7->Get("n_clus_hits");
	IDc7sample[2]->SetTitle(IDneutrg[2]->GetTitle());
	
	IDc7sample[3] = (TH1F*)can7->Get("charge");
	IDc7sample[3]->SetTitle(IDneutrg[3]->GetTitle());
	for(int i = 0; i < 4; i++){IDc7sample[i]->SetDirectory(0);}
	can7->Close();

	
	// pulse shape: canvas 05
	IDtimes[0] = new TH1F("IDlaserStartTime", "Laser (TT8)", 9000/5, -1000, 8000);
	IDtimes[1] = new TH1F("IDneutrinoStartTime", "Neutrino (TT1&BTB0)", 18000/5, -17000, 1000);
	IDtimes[2] = new TH1F("IDpulserStartTime", "Pulser (TT32)", 18000/5, -17000, 1000);
	IDtimes[3] = new TH1F("IDmuonTT1StartTime", "Internal muon (TT1&BTB4)", 18000/5, -17000, 1000);
	IDtimes[4] = new TH1F("IDrandomStartTime", "Random (TT64)", 18000/5, -17000, 1000);
	IDtimes[5] = new TH1F("IDmuonTT2StartTime", "External muon (TT2)", 18000/20, -17000, 1000); // different binning so that it's visible that it's flat when truncated

	// event size VS channels: canvas 06
	IDchannels[0] = new TH1F("IDlaserChannels", "Laser (TT8)", 2240, 1, 2241);
	IDchannels[1] = new TH1F("IDpulserChannels", "Pulser (TT32)", 2240, 1, 2241);
	IDchannels[2] = new TH1F("IDrandomChannels", "Random (TT64)", 2240, 1, 2241);
	IDchannels[3] = new TH1F("IDneutrinoChannels", "Internal events (TT1)", 2240, 1, 2241);
	
	// hits in muon triggers: canvas 12
	IDhitsMuon[0] = new TH1F("IDmuonTT1raw", "Internal muons (TT1&BTB4) raw hits", 200, 0, 50000);
	IDhitsMuon[1] = new TH1F("IDmuonTT1decoded", "Internal muons (TT1&BTB4) decoded hits", 200, 0, 50000);
	IDhitsMuon[2] = new TH1F("IDmuonTT1clustered", "Internal muons (TT1&BTB4) clustered hits", 200, 0, 50000);
	IDhitsMuon[3] = new TH1F("IDmuonTT2clusters", "Clusters muon TT2", 4, 0, 4);
	IDhitsMuon[4] = new TH1F("IDmuonTT1clusters", "Clusters muonTT1", 4, 0, 4);

	// neutron trigger: canvas 16
	IDtrg128[0] = new TH1F("IDneutronRaw", "Raw hits",100,0,2000);
	IDtrg128[1] = new TH1F("IDneutronDecoded", "Decoded hits",100,0,2000);
	IDtrg128[2] = new TH1F("IDneutronClustered1", "Clustered hits in 1-cluster events",200,0,400);
	IDtrg128[3] = new TH1F("IDneutronClustered", "Hits in clusters of multi-cluster events",200,0,400);
	IDtrg128[4] = new TH1F("IDneutronCharge","Charge in each cluster",200,0,400);
	IDtrg128[5] = new TH1F("IDclusterTimeNn","Cluster start time distribution", 1880/20, -1800, 80); // us
	IDtrg128[6] = new TH1F("IDneutronClusters", "Clusters",12,0,12);
	IDtrg128[7] = new TH1F("IDemptyBoards","Empty boards", 200,0,200);
	
	// neutron trigger, canvas 17 
	IDtrg128[8] = new TH1F("IDneutronStartTime","Neutron (TT128) pulse shape", 1800/5, -1700, 10); //us
	IDtrg128[9] = new TH1F("IDneutronStartTimeSec","Neutron-muon pulse shape", 30000/2, -17000, -17000+30000); //ns
	IDtrg128[10] = new TH1F("IDneutronMuonTime","Neutron-muon pulse shape", 30000/2, -17000, -17000+30000); //ns


	// ***** OUTER DETECTOR *****//
	
	// OD flag stability: canvas 08
	ODmuonFlag[0] = new TH1F("ODmuonFlagTT1MTB","MTB",10, 0,10);
	ODmuonFlag[1] = new TH1F("ODmuonFlagTT1MCR","MCR",10, 0,10);
	ODmuonFlag[2] = new TH1F("ODmuonFlagTT1ID","IDF",10, 0,10);
	ODmuonFlag[3] = new TH1F("ODmuonFlagTT2MTB","MTB",10, 0,10);
	ODmuonFlag[4] = new TH1F("ODmuonFlagTT2MCR","MCR",10, 0,10);
	ODmuonFlag[5] = new TH1F("ODmuonFlagTT2MTB25","MTB ev.size.>25",10, 0,10);

	// event size VS evnum: canvas 09
	ODevsize[0] = new TH1F("ODlaserEvSize", "Laser (TT8)", end_evnum - start_evnum, start_evnum, end_evnum);
	ODevsize[1] = new TH1F("ODpulserEvSize", "Pulser (TT32)", end_evnum -  start_evnum, start_evnum, end_evnum);
	ODevsize[2] = new TH1F("ODrandomEvSize", "Random (TT64)", end_evnum -  start_evnum, start_evnum, end_evnum);
	ODevsize[3] = new TH1F("ODneutrinoEvSize", "Neutrino (TT1&BTB0)", end_evnum - start_evnum, start_evnum, end_evnum);
	ODevsize[4] = new TH1F("ODmuonTT1EvSize", "Internal Muon (TT1&BTB4)", end_evnum - start_evnum, start_evnum, end_evnum);
	ODevsize[5] = new TH1F("ODmuonTT2EvSize", "External Muon (TT2)", end_evnum - start_evnum, start_evnum, end_evnum);
	pulserODylim[0] = 100000; // y-min
	pulserODylim[1] = -10; // y-max
	
	// event size VS evnum in 200 bins: canvas 10 
	ODevsize200[0] = new TH1F("ODlaserEvSize200", "Laser (TT8)", 200, 0, 200);
	ODevsize200[1] = new TH1F("ODpulserEvSize200", "Pulser (TT32)", 200, 0, 200);
	ODevsize200[2] = new TH1F("ODrandomEvSize200", "Random (TT64)", 200, 0, 200);
	ODevsize200[3] = new TH1F("ODneutrinoEvSize200", "Neutrino (TT1&BTB0)", 200, 0, 200);
	ODevsize200[4] = new TH1F("ODmuonTT1EvSize200", "Internal muon (TT1&BTB4)", 200, 0, 200);
	ODevsize200[5] = new TH1F("ODmuonTT2EvSize200", "External muon(TT2)", 200, 0, 200);

	// muon trigger rate VS time: canvas 11
	// max time: 10 hours, 10-minute bins
	ODtrgRate[0] = new TH1F("ODtrgRateAlltt1", "Internal muon (TT1&BTB4) all", 3600*10/600, 0, 3600*10);
	ODtrgRate[1] = new TH1F("ODtrgRateAlltt2", "External muon (TT2) all", 3600*10/600, 0, 3600*10 );
	ODtrgRate[2] = new TH1F("ODtrgRate100tt1", "Internal muon (TT1&BTB4) ev.size > 25", 3600*10/600, 0, 3600*10 );
	ODtrgRate[3] = new TH1F("ODtrgRate100tt2", "External muon (TT2) ev.size > 25", 3600*10/600, 0, 3600*10 );
	for(int i = 0; i < 4; i++){ ODtrgRate[i]->Sumw2(); }

	// hits in muon triggers: canvas 13
	ODhitsTT1[0] = new TH1F("ODrawTT1","Raw hits", 200, 0, 200);
	ODhitsTT1[1] = new TH1F("ODdecodedTT1","Decoded hits", 200, 0, 200);
	ODhitsTT1[2] = new TH1F("ODclusteredTT1","Clustered hits", 200, 0, 200);
	ODhitsTT1[3] = new TH1F("ODclustersTT1","Number of clusters", 14, 0, 14);
	
	ODhitsTT2[0] = new TH1F("ODrawTT2","Raw hits", 200, 0, 200);
	ODhitsTT2[1] = new TH1F("ODdecodedTT2","Decoded hits", 200, 0, 200);
	ODhitsTT2[2] = new TH1F("ODclusteredTT2","Clustered hits", 200, 0, 200);
	ODhitsTT2[3] = new TH1F("ODclustersTT2","Clusters in muons", 14, 0, 14);

	// pulse shape: canvas 14
	ODtimes[0] = new TH1F("ODlaserStartTime", "Laser (TT8)", 9000/5, 0, 9000);
	ODtimes[1] = new TH1F("ODneutrinoStartTime", "Neutrino (TT1&BTB0)", 18000/5, -18000, 0);
	ODtimes[2] = new TH1F("ODpulserStartTime", "Pulser (TT32)", 18000/5, -18000, 0);
	ODtimes[3] = new TH1F("ODmuonTT1StartTime", "Internal muon (TT1&BTB4)", 18000/5, -18000, 0);
	ODtimes[4] = new TH1F("ODrandomStartTime", "Random (TT64)", 18000/5, -18000, 0);
	ODtimes[5] = new TH1F("ODmuonTT2StartTime", "External muon (TT2)", 18000/5, -18000, 0);

	// event size VS channels: canvas 15
	ODchannels[0] = new TH1F("ODlaserChannels","Laser (TT8)",255,1,256);
	ODchannels[1] = new TH1F("ODneutrinoChannels","Neutrino (TT1&BTB0)",255,1,256);
	ODchannels[2] = new TH1F("ODpulserChannels","Pulser (TT32)",255,1,256);
	ODchannels[3] = new TH1F("ODmuonTT1Channels","Internal muon (TT1&BTB4)",255,1,256);
	ODchannels[4] = new TH1F("ODrandomChannels","Random (TT64)",255,1,256);
	ODchannels[5] = new TH1F("ODmuonTT2Channels","External muon (TT2)",255,1,256);

}


// ------------------------------------------------------------

void BxValidation::FillHistos(){
	/*
	 * Fill the histograms related to the validation of the run and save some useful values to the BX VALIDATION object for later
	 */

//	BxEvent* ev = new BxEvent();
//	t->SetBranchAddress("events",&ev);
//	oink


	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	cout << "Start event (first neutrino): " << start_evnum << endl;
	cout << "End event: "<< end_evnum  << endl;
	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	cout << "Filling the histograms..." << endl;
	

	//********************** variables that go throughout the loop
	// muon time to remember for TT128 events
	double muon_time = 0;
	// for filling muon rates
	double muonFactor = 10./(end_time - start_time)*60.*1440.;
	// normalizing the 200-bin histograms number of small bins in one big bin
	binstep = (int) ceil((end_evnum - start_evnum)*1./200);

//	cout << "Canvases 2 and 10: big bin every " << binstep << " events" << endl;
	
	double bincontent = 0; // bin content for the 200-bin histos to be normalized later

	// 200-bin histo normalization factors:
	
	// counting how many non-zero small bins fall in a big bin
	int countOD[6];
	for(int i = 0; i < 6; i++){countOD[i] = 0;}
	// counting how many non-zero small bins fall in a big bin
	int countID[6]; 
	for(int i = 0; i < 6; i++){countID[i] = 0;}
	
	// counting TT2 events to check OD behaviour later
	TT2count = 0;
	// counting pulser and TT1 events to rescale electronics calibration histograms later
	TT32count = 0; bx->TT1count = 0;
	// counting TT1 muons (TT1&BTB4) to later draw sigma lines
	muonTT1count = 0;
	TT128count = 0;

	static const long int gray_window = (1 << 16) * 50;
	
	//********************** LOOP OVER EVENTS

	for(int evnum = start_evnum; evnum < end_evnum; evnum++){
//		cout << evnum << endl;
		t->GetEntry(evnum); // index is equivalent to the event number here
	
		// ************ COMMON VARIABLES
	
		const BxTrigger& trg = ev->GetTrigger();
		const BxLaben& laben = ev->GetLaben();
		const BxMuon& muon = ev->GetMuon();
		const BxLabenCluster& cluster = laben.GetCluster(0);
		
		int trgtype = trg.GetTrgType(); 
		double trgtime = laben.GetTriggerTime(); // ns
		double lasertime = laben.GetLaserTime(); // time of laser reference
    double ev_time = trg.GetGpsTimeSec();
		double ev_time_ns = trg.GetGpsTimeNs();
		double cluster_time = cluster.GetStartTime();
		
		int n_dec_hits_ID = laben.GetNDecodedHits();
		int n_raw_hits_ID = laben.GetNRawHits();
		int n_clustered_hits_ID = laben.GetNClusteredHits();
		int n_clusters_ID = laben.GetNClusters();
		float charge_ID = cluster.GetCharge();

		int n_dec_hits_OD = muon.GetNDecodedHits();
		int n_raw_hits_OD = muon.GetNRawHits();
		int n_clustered_hits_OD = muon.GetNClusteredHits();
		int n_clusters_OD = muon.GetNClusters();

		float muonBin = (ev_time - start_time)/(end_time - start_time)*9.999; // if x10, the last even will go into overflow

		// ************ NORMALIZATION OF 200 BIN HISTOS

		// current bin number for the histo with 200 bins
		// note: bindex starts from 0 (what we fill), bin number is bindex+1 (0 gets filled into bin #1)
		int bindex = (evnum - start_evnum) / binstep;
		if(evnum == end_evnum-1){ bindex = 200; }
		
		// normalize the big binned histos 
		// from canvas 02 and canvas 10
		// if it's not the first bin and it's the first time we touch this bin, normalize the previous bin unless it's empty
		// or, when we hit the end, bindex is set to 200 to normalize the last possibly unfull bin
		if( (bindex && !( (evnum - start_evnum) % binstep)) || bindex == 200 ){
			// inner detector
			for(int i = 0; i < 6; i++){
				if(countID[i]){
					// normalized bin content
					bincontent = IDevsize200[i]->GetBinContent(bindex)/countID[i];
					IDevsize200[i]->SetBinContent(bindex, bincontent);
					// poisson error
					IDevsize200[i]->SetBinError(bindex, sqrt(bincontent/countID[i]));
					// pulser is not poisson
					if(i == 1){
						IDevsize200[i]->SetBinError(bindex,sqrt(bincontent)/5./sqrt(countID[i]));
					}
				}
				// reset the counter for the next bin
				countID[i] = 0;
			}

			// outer detector
			for(int i = 0; i < 6; i++){
				if(countOD[i]){
					// normalized bin content
					bincontent = ODevsize200[i]->GetBinContent(bindex)/countOD[i];
					ODevsize200[i]->SetBinContent(bindex, bincontent);
					// poisson error
					ODevsize200[i]->SetBinError(bindex, sqrt(bincontent/countOD[i]));
					// pulser is not poisson
					if(i == 1){
						ODevsize200[i]->SetBinError(bindex,sqrt(bincontent)/5./sqrt(countOD[i]));
					}
				}
				// reset the counter for the next bin
				countOD[i] = 0;
			}
		} // if bindex

		// ************ TRIGGER TYPE SWITCH
		
		switch(trgtype){
			
			// LASER // TT8 //
			case 8: {
								// event size VS evnum
								countID[0]++;
								// ID: canvas 01
								IDevsize[0]->Fill(evnum, n_dec_hits_ID);
							
								// ID: canvas 02
								IDevsize200[0]->Fill(bindex,n_dec_hits_ID);
								
								for(int i = 0; i < n_dec_hits_ID;i++){
									// ID pulse shape: canvas 05
									IDtimes[0]->Fill(laben.GetDecodedHits()[i].GetRawTime() - lasertime);
									// ID event size VS channels: canvas 06
									IDchannels[0]->Fill(laben.GetDecodedHits()[i].GetLg());
								}
							
								// OD: canvas 09
								ODevsize[0]->Fill(evnum, n_dec_hits_OD);
								
								// OD: canvas 10
								ODevsize200[0]->Fill(bindex,n_dec_hits_OD);
								countOD[0]++;


								for(int i = 0; i < n_dec_hits_OD; i++){
									// OD pulse shape: canvas 14
									ODtimes[0]->Fill(muon.GetDecodedHits()[i].GetTime());
									// OD event size VS channels: canvas 15
									ODchannels[0]->Fill(muon.GetDecodedHits()[i].GetMch());
								}
								
								break;
							} // laser TT8

			// PULSER // TT32 //
			case 32: {
								 // counting pulser events to rescale electronics calibration histograms later
								 TT32count++;
								 countID[1]++;
							
								 // event size VS evnum
								 // ID: canvas 01
								 IDevsize[1]->Fill(evnum, n_dec_hits_ID);
								 // y-axis limit for later
								 if(n_dec_hits_ID > pulserIDylim[1]){
									 pulserIDylim[1] = n_dec_hits_ID;
								 }
								 else if(n_dec_hits_ID < pulserIDylim[0]){
									 pulserIDylim[0] = n_dec_hits_ID;
								 }

								 // ID: canvas 02
								 
								 IDevsize200[1]->Fill(bindex,n_dec_hits_ID);
								
								 for(int i = 0; i < n_dec_hits_ID;i++){
									 // ID pulse shape: canvas 05
									 IDtimes[2]->Fill(laben.GetDecodedHits()[i].GetRawTime() - trgtime);
									 // ID event size VS channels: canvas 06
									 IDchannels[1]->Fill(laben.GetDecodedHits()[i].GetLg());
								 }

								 // OD: canvas 09
								 ODevsize[1]->Fill(evnum, n_dec_hits_OD);
								 
								 // y-axis limit for later
								 if(n_dec_hits_OD > pulserODylim[1]){
									 pulserODylim[1] = n_dec_hits_OD;
								 }
								 else if(n_dec_hits_OD < pulserODylim[0]){
									 pulserODylim[0] = n_dec_hits_OD;
								 }
								 
								 // OD: canvas 10
								 ODevsize200[1]->Fill(bindex,n_dec_hits_OD);
								 countOD[1]++;
							

								for(int i = 0; i < n_dec_hits_OD; i++){
									// OD pulse shape: canvas 14
									ODtimes[2]->Fill(muon.GetDecodedHits()[i].GetTime());
									// OD event size VS channels: canvas 15
									ODchannels[2]->Fill(muon.GetDecodedHits()[i].GetMch());
								}
								
								break;
							 } // pulser TT32

			// RANDOM // TT64 //
			case 64: {
								 // event size VS evnum
								 countID[2]++;
								 // ID: canvas 01
								 IDevsize[2]->Fill(evnum, n_dec_hits_ID);
								 
								 // ID: canvas 02
								 IDevsize200[2]->Fill(bindex,n_dec_hits_ID);								
								 
								 for(int i = 0; i < n_dec_hits_ID;i++){
									 // ID pulse shape: canvas 05
									 IDtimes[4]->Fill(laben.GetDecodedHits()[i].GetRawTime() - trgtime);
									 
									 // ID event size VS channels: canvas 06
									 IDchannels[2]->Fill(laben.GetDecodedHits()[i].GetLg());
								 }
								 
								 // OD: canvas 09
								 ODevsize[2]->Fill(evnum, n_dec_hits_OD);
								 
								 // OD: canvas 10
								 ODevsize200[2]->Fill(bindex,n_dec_hits_OD);
								 countOD[2]++;
	

								 for(int i = 0; i < n_dec_hits_OD; i++){
									 // OD pulse shape: canvas 14
									 ODtimes[4]->Fill(muon.GetDecodedHits()[i].GetTime());
									 // OD event size VS channels: canvas 15
									 ODchannels[4]->Fill(muon.GetDecodedHits()[i].GetMch());
								}
								
								break;
							 } // random TT64


			// TT1 // ID yes OD whatever //
			case 1: {
								// counting TT1 events for electronics rescaling later
								TT1count++;

								// ID trigger rate: canvas 03
								// not to fill for the events from the last bin (if I do not have the whole minute of data)
								if( (ev_time - start_time) < (int)((end_time - start_time)/60.)*60 ){
									IDtrgRate[0]->Fill(ev_time - start_time,1/60.);
									if(n_dec_hits_ID > 100){
										IDtrgRate[1]->Fill(ev_time - start_time,1/60.);
									}
								}
					
								for(int i = 0; i < n_dec_hits_ID;i++){
									// ID event size VS channels: canvas 06
									IDchannels[3]->Fill(laben.GetDecodedHits()[i].GetLg());
								}
								
								// OD muon flag stability: canvas 08
								bool IDflag = false;
//									cout << n_clusters_ID << " " << n_clustered_hits_ID << " " << cluster.GetPeakTime(0) << endl;
								// if n_clusters_ID = 0, n_clustered_hits_ID is automatically 0 --> IDflag stays false
								if(n_clustered_hits_ID > 2100){
//									cout << n_clusters_ID << endl;
									IDflag = (cluster.GetMeanTime() > 100) && (laben.GetRecClusters()[0].GetGatti() < 0.55);
								}

								else if(n_clustered_hits_ID > 900){
//									cout << n_clusters_ID << " " << n_clustered_hits_ID << endl;
//									IDflag = cluster.GetPeakTimes()[0] > 30.;
									IDflag = cluster.GetPeakTime(0) > 30.;
//									return;
								}

								else if( n_clustered_hits_ID > 100){
//									cout << n_clusters_ID << " " << n_clustered_hits_ID << endl;
//									IDflag = cluster.GetPeakTimes()[0] > 40.;
									IDflag = cluster.GetPeakTime(0) > 40.;
//									return;
								}
								
								// TT1 ID
								if(IDflag){
									ODmuonFlag[2]->Fill(muonBin, muonFactor); 
								}
								
								// MCR
								if(n_clusters_OD > 0){
									// OD muon flag stability: canvas 8
									// TT1 MCR
									ODmuonFlag[1]->Fill(muonBin, muonFactor); 
								}
							

								// NEUTRINO event // TT1 & BTB0 // ID yes OD no //
								if(trg.GetBtbInputs() == 0){						
									// event size vs evnum
									countID[3]++;
									// ID event size: canvas 01
									IDevsize[3]->Fill(evnum, n_dec_hits_ID);
								
									// ID: canvas 02								 
									IDevsize200[3]->Fill(bindex,n_dec_hits_ID);
									// ID neutrino trigger: canvas 04, canvas 07
									IDneutrg[0]->Fill(n_raw_hits_ID);
									IDneutrg[1]->Fill(n_dec_hits_ID);
									IDneutrg[4]->Fill(n_clusters_ID);
									// time difference between neutrinos
									IDneutrg[5]->Fill((ev_time_ns - previous_neu_time)/1000000.); // in ms
									IDneutrg[8]->Fill((ev_time_ns - previous_neu_time)/1000.); // in mus

									previous_neu_time = ev_time_ns;

									if(n_clusters_ID == 1){
										IDneutrg[2]->Fill(n_clustered_hits_ID);
										IDneutrg[7]->Fill(n_dec_hits_ID);
										IDneutrg[3]->Fill(charge_ID);
										IDneutrg[6]->Fill(cluster_time - trgtime);
									}

									// ID pulse shape: canvas 05
									for(int i = 0; i < n_dec_hits_ID;i++){
										IDtimes[1]->Fill(laben.GetDecodedHits()[i].GetRawTime() - trgtime);
									}

									// OD event size vs evnum: canvas 09
									ODevsize[3]->Fill(evnum, n_dec_hits_OD);
								
									// OD: canvas 10
									ODevsize200[3]->Fill(bindex,n_dec_hits_OD);
									countOD[3]++;
									
									for(int i = 0; i < n_dec_hits_OD; i++){
									  // OD pulse shape: canvas 14
										ODtimes[1]->Fill(muon.GetDecodedHits()[i].GetTime());
									  // OD event size VS channels: canvas 15
										ODchannels[1]->Fill(muon.GetDecodedHits()[i].GetMch());
									}
								} // if for NEUTRINO i.e. TT1&BTB0


								// INTERNAL MUON // TT1 & BTB4 // ID yes OD yes //
								else if(trg.GetBtbInputs() & 4){

								 // event size VS evnum
								 countID[4]++;
								 // ID: canvas 01
								 IDevsize[4]->Fill(evnum, n_dec_hits_ID);
								 
								 // ID: canvas 02
								 IDevsize200[4]->Fill(bindex,n_dec_hits_ID);
									
								 for(int i = 0; i < n_dec_hits_ID;i++){
								   // ID pulse shape: canvas 05
									 IDtimes[3]->Fill(laben.GetDecodedHits()[i].GetRawTime() - trgtime);
									 // ID neutron-muon pulse shape: canvas 17
									 IDtrg128[10]->Fill(laben.GetDecodedHits()[i].GetRawTime() - trgtime);
								 }

									// OD muon flag stability: canvas 08
									// count muons
									muonTT1count++;
									// TT1 MTB
									ODmuonFlag[0]->Fill(muonBin, muonFactor); 
									
									// OD event size vs evnum: canvas 09
									ODevsize[4]->Fill(evnum, n_dec_hits_OD);
								
									// OD: canvas 10
									ODevsize200[4]->Fill(bindex,n_dec_hits_OD);
									countOD[4]++;

									// OD trigger rate: canvas 11
									// rate is per day -> *24; 10-minute bins -> *6)
									//if( (ev_time - start_time) < (int)((end_time - start_time)/60.)*60){
									if( (ev_time - start_time) < (int)( (end_time - start_time)/600. ) * 600  ){
										ODtrgRate[0]->Fill(ev_time - start_time, 24*6.);
										if (n_dec_hits_OD > 25){
											ODtrgRate[2]->Fill(ev_time - start_time, 24*6.);
										}
									}
									
									// ID hits in muon triggers: canvas 12
									IDhitsMuon[0]->Fill(n_raw_hits_ID);
									IDhitsMuon[1]->Fill(n_dec_hits_ID);
									IDhitsMuon[2]->Fill(n_clustered_hits_ID);
									IDhitsMuon[4]->Fill(n_clusters_ID);
									
									// OD hits in muon triggers: canvas 13
									ODhitsTT1[0]->Fill(n_raw_hits_OD);
									ODhitsTT1[1]->Fill(n_dec_hits_OD);
									// the zero-clustered hits are crowding the space
									if(n_clustered_hits_OD){ ODhitsTT1[2]->Fill(n_clustered_hits_OD); }
									ODhitsTT1[3]->Fill(n_clusters_OD);
									
									for(int i = 0; i < n_dec_hits_OD; i++){
									  // OD pulse shape: canvas 14
										ODtimes[3]->Fill(muon.GetDecodedHits()[i].GetTime());
									  // OD event size VS channels: canvas 15
										ODchannels[3]->Fill(muon.GetDecodedHits()[i].GetMch());
									}

									 // canvas 17: remember the muon trigger time for the neutron afterwards 
									 muon_time = trgtime;

								} // if for INTERNAL MUON i.e. TT1&BTB4

								break;
							} // TT1
			
							// EXTERNAL MUON // TT2 // ID no OD yes //
			case 2: {
								// counting TT2 events for checking OD behaviour later
								TT2count ++;
								
								// ID event size VS evnum
								countID[5]++;
								// canvas 01
								IDevsize[5]->Fill(evnum, n_dec_hits_ID);								 
								// canvas 02
								IDevsize200[5]->Fill(bindex,n_dec_hits_ID);
								// ID pulse shape: canvas 05
								for(int i = 0; i < n_dec_hits_ID;i++){
									IDtimes[5]->Fill(laben.GetDecodedHits()[i].GetRawTime() - trgtime);
								}

								// OD muon flag stability: canvas 08
								// TT2 MTB
								ODmuonFlag[3]->Fill(muonBin, muonFactor); 
								// TT2 MCR
								if(n_clusters_OD > 0){
									ODmuonFlag[4]->Fill(muonBin, muonFactor);
								}
								// TT2 MTB > 25 hits
								if(n_dec_hits_OD > 25){
									ODmuonFlag[5]->Fill(muonBin, muonFactor);
								}

								// OD event size VS evnum: canvas 09
								ODevsize[5]->Fill(evnum, n_dec_hits_OD);
								
								// OD: canvas 10
								ODevsize200[5]->Fill(bindex,n_dec_hits_OD);
								countOD[5]++;
								
								// OD trigger rate: canvas 11
								// rate is per day -> *24; 10-minute bins -> *6)
									//if( (ev_time - start_time) < (int)((end_time - start_time)/60.)*60){
									if( (ev_time - start_time) < (int)( (end_time - start_time)/600. ) * 600  ){
										ODtrgRate[1]->Fill(ev_time - start_time, 24*6.);
										if (n_dec_hits_OD > 25){
											ODtrgRate[3]->Fill(ev_time - start_time, 24*6.);
										}
									}
									
								// ID hits in muon triggers: canvas 12
								IDhitsMuon[3]->Fill(n_clusters_ID);
									
							  // OD hits in muon triggers: canvas 13
								ODhitsTT2[0]->Fill(n_raw_hits_OD);
								ODhitsTT2[1]->Fill(n_dec_hits_OD);
								if(n_clustered_hits_OD){ ODhitsTT2[2]->Fill(n_clustered_hits_OD); }
								ODhitsTT2[3]->Fill(n_clusters_OD);

								for(int i = 0; i < n_dec_hits_OD; i++){
									// OD pulse shape: canvas 14
									ODtimes[5]->Fill(muon.GetDecodedHits()[i].GetTime());
									// OD event size VS channels: canvas 15
									ODchannels[5]->Fill(muon.GetDecodedHits()[i].GetMch());
								}
									
								break;
							} // case 2

			// NEUTRON // TT128
			case 128: {

								TT128count++;
									
								if(muon_time == 0){ cout << "?? Neutron (TT128) without a muon ??" << endl;}

								// neutron trigger: canvas 16
								
								IDtrg128[0]->Fill(n_raw_hits_ID);
								IDtrg128[1]->Fill(n_dec_hits_ID);
								IDtrg128[6]->Fill(n_clusters_ID);
								IDtrg128[7]->Fill(laben.GetEmptyBoards());
								
								// if 0 clusters won't even go to this loop
								for(int i = 0; i < n_clusters_ID; i++){
									// for each cluster
									IDtrg128[4]->Fill(charge_ID);
									IDtrg128[5]->Fill( (laben.GetClusters()[i].GetStartTime() - trgtime)/1000.); // us
							
									// one cluster
									if(n_clusters_ID == 1){
										IDtrg128[2]->Fill(n_clustered_hits_ID);
									}
								
									// multicluster
									else{
										IDtrg128[3]->Fill(laben.GetClusters()[i].GetNHits());
									}
								}
								
								// canvas 17
								// correction if the clock failed
								int correction = trgtime - muon_time < 0 ? gray_window : 0; // ns
								
								for(int i = 0; i < n_dec_hits_ID; i++){
									double hittime = laben.GetDecodedHits()[i].GetRawTime(); //ns
									IDtrg128[8]->Fill( (hittime - trgtime)/1000.); //us
									IDtrg128[9]->Fill(hittime - muon_time + correction); //ns
								}

								muon_time = 0;

								break;
							}

		} // switch trigger type
		
		// occasional print statement
		if(evnum % 10000 == 0){
			cout << "Event " << evnum << "/" << end_evnum << endl;
		}
	
	} // for loop over events

	// last even done
	cout << "Event " << end_evnum << "/" << end_evnum  << endl;
	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;


	// TT1 muons = TT128 neutrons check
	if(muonTT1count != TT128count){
		cerr << "WARNING: # TT1 muons (" << bx->muonTT1count << ") != # TT128 (" << bx->TT128count << ")!" << endl;
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	}

}


// ------------------------------------------------------------

void BxValidation::DrawCanvases(){
	/*
	 * Draw the canvases from the histos that were filled with Fill Histos
	 */

	cout << "......prepare for canvas overload......" << endl;
	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	

	// recyclable variables
	// fitting
	TF1 *fit; 
	// axis range
	TH1* frame;
	// sigma lines
	TLine* sigma;
	// label
	TPaveLabel *label, *label1;
	// label text
	string lbl;
	// text
	TPaveText* text;
	// any kind of char that's needed, can be reused
	char txt[100];

	// canvas size for PNG
//	int ww = 1000; int hh = 800;
	int ww = 800; int hh = 600;
	//int ww = 900; int hh = 600;
	//	int ww = 1200; int hh = 800;




	//******************** canvas 01 ********************//
	//ID: event size VS evnum
	
	cout << "canvas 01" << endl;
		
	c[1] = new TCanvas("c1","[01] ID: event size VS evnum",ww,hh);
	c[1]->Divide(3,2); 

	// pad 2: pulser y-axis limit
	IDevsize[1]->SetAxisRange(pulserIDylim[0]*0.95, pulserIDylim[1]*1.05, "Y");

	for(int i = 0; i < 6; i++){
		c[1]->cd(i+1);
		
		IDevsize[i]->SetXTitle("Event number");
		IDevsize[i]->SetYTitle("Event size in ID (N decoded hits)");
	
		if(IDevsize[i]->GetMaximum() > 999){
//			bx->IDevsize[i]->SetLabelOffset(-0.02,"Y");
			IDevsize[i]->SetTitleOffset(1.05,"Y");
		}
		
		IDevsize[i]->Draw("simple");
		
		// fit
		fit = new TF1("line", "pol1", 0., IDevsize[i]->GetSize());
		IDevsize[i]->Fit("line","q");	
	}
	
	//**************************************************//


	//******************** canvas 02 ********************//
	//  ID: event size VS evnum (200 bins)

	cout << "canvas 02" << endl;
		
	c[2] = new TCanvas("c2","[02] ID: event size VS evnum (200 bins)",ww,hh);
	c[2]->Divide(3,2);	
	
	sprintf(txt,"Ev. number in 200 bins (1 bin = %d events)", binstep);
	

	for(int i = 0; i < 6; i++){
		c[2]->cd(i+1);

		IDevsize200[i]->SetXTitle(txt);
		IDevsize200[i]->SetYTitle("Avg event size in ID (N decoded hits)");
		
		if(IDevsize200[i]->GetMaximum() > 999){
//			bx->IDevsize200[i]->SetLabelOffset(-0.02,"Y");
			IDevsize200[i]->SetTitleOffset(1.05,"Y");
		}

		IDevsize200[i]->Draw("simple"); 

		// fit
		fit = new TF1("line", "pol1", 0., IDevsize200[i]->GetSize());
		IDevsize200[i]->Fit("line","q");
	}

	//**************************************************//


	//******************** canvas 03 ********************//
	// ID: neutrino trigger rate VS time

	cout << "canvas 03" << endl;
		
	c[3] = new TCanvas("c3","[03] ID: trigger rate",ww,hh);
	c[3]->Divide(1,2);
	
	for(int i = 0; i < 2; i++){
		c[3]->cd(i+1);
	
		IDtrgRate[i]->SetAxisRange(0,end_time - start_time);
		IDtrgRate[i]->SetXTitle("Event time (1 bin = 1 min) [s]");
		IDtrgRate[i]->SetYTitle("Trigger rate (per second)");
		IDtrgRate[i]->Draw("simple");
		// fit
		fit = new TF1("line", "pol1", 0., IDtrgRate[i]->GetSize());
		IDtrgRate[i]->Fit("line","q");		
	}
	
	//**************************************************//

	
	//******************** canvas 04 ********************//
	// ID: hits in neutrino trigger 
	cout << "canvas 04" << endl;
		
	c[4] = new TCanvas("c4","[04] ID: neutrino (TT1&BTB0)",ww,hh);
	c[4]->Divide(3,2);

	//// ------- pad 1
	c[4]->cd(1);
	IDneutrg[4]->SetXTitle("Number of clusters in ID"); 
	IDneutrg[4]->Draw(); 
	c[4]->GetPad(1)->SetLogy();

	// percentage of 0-clusters
	sprintf(txt,"0-clusters: %.1f%%",IDneutrg[4]->GetBinContent(1)*1.0/IDneutrg[4]->GetEntries()*100);
	text = new TPaveText(0.7,0.6,0.75,0.8, "NDC");
	text->AddText(txt);
	text->SetFillColor(kWhite);
	text->SetTextSize(0.04);
	text->SetBorderSize(0);
	text->Draw();
	/// -------------

	// common label for pads 2 and 3
	lbl = "Decoded Hits";
	label1 = new TPaveLabel(0.7,0.55,0.8,0.65,lbl.c_str(),"NDC");
	label1->SetBorderSize(0);
	label1->SetTextSize(0.4);
	label1->SetFillColor(kWhite);
	label1->SetTextColor(kRed);

	////------- pad 2: Raw and decoded
	c[4]->cd(2);

	// axes limits
	frame = c[4]->GetPad(2)->DrawFrame(0,0,120, IDneutrg[1]->GetMaximum()*1.05);
	frame->SetXTitle("Number of hits in ID");
	frame->SetTitle("Raw and decoded hits in neutrino (TT1&BTB0) events");	

	TH1F* htemp = (TH1F*)IDneutrg[1]->Clone();
	htemp->SetDirectory(0);
	htemp->SetLineColor(kRed);
	htemp->Draw("sames"); // not to overwrite the frame
//	bx->IDneutrg[1]->SetLineColor(kRed);
//	bx->IDneutrg[1]->Draw("sames"); // not to overwrite the frame
	IDneutrg[0]->Draw("sames");
		
	lbl = "Raw hits";
	label = new TPaveLabel(0.7,0.45,0.8,0.55,lbl.c_str(),"NDC");
	label->SetBorderSize(0);
	label->SetTextColor(kBlue+3);
	label->SetTextSize(0.4);
	label->SetFillColor(kWhite);
	label->Draw();

	label1->Draw();
	/// -------------

	////------- pad 3: Decoded and clustered
	c[4]->cd(3);

	// axes limits: use DrawFrame rather than SetAxisRange not to change canvas 7 which uses the same histos
	frame = c[4]->GetPad(3)->DrawFrame(0,0,100, IDneutrg[2]->GetMaximum()*1.05);
	frame->SetXTitle("Number of hits in ID");
	frame->SetTitle("Hits in 1-cluster events");	
	
	IDneutrg[2]->Draw("sames"); // not to overwrite the frame
	IDneutrg[7]->SetLineColor(kRed);
	IDneutrg[7]->Draw("sames");
		
	lbl = "Clustered hits";
	label = new TPaveLabel(0.7,0.45,0.8,0.55,lbl.c_str(),"NDC");
	label->SetBorderSize(0);
	label->SetTextSize(0.4);
	label->SetFillColor(kWhite);
	label->SetTextColor(kBlue+3);
	label->Draw();

	label1->Draw();
	/// -------------
	

	////------- pad 4: cluster time

	c[4]->cd(4);
	IDneutrg[6]->SetXTitle("Cluster start time (trigger time = 0) [ns]");
	IDneutrg[6]->Draw();
	c[4]->GetPad(4)->SetLogy();

	// fit
	fit = new TF1("fit", SkewGaus, -15500, -14000,4);
	float c0mean = IDneutrg[6]->GetMean();
	float c0rms = bx->IDneutrg[6]->GetRMS();

	fit->SetParameter(0, IDneutrg[6]->GetEntries());
	fit->SetParameter(1, c0mean); // mu
	fit->SetParameter(2, c0rms); // sigma left
	fit->SetParameter(3, c0rms); // sigma right
	fit->SetParLimits(1, c0mean - 3*c0rms, c0mean + 3*c0rms);
	fit->SetParLimits(2, 0., 1000.);
	fit->SetParLimits(3, 0., 1000.);
//	line->SetParNames("normalization","mean value","left sigma","right sigma");
	fit->SetLineColor(kRed); 
//	fit->SetLineWidth(2);      
	IDneutrg[6]->Fit(fit, "QEMR+");   

	// text on the plot
	text = new TPaveText(0.15,0.25,0.7,0.6, "NDC");
	// clusters before the gate
	int clus = (int) IDneutrg[6]->Integral(0, int(-16000 - IDneutrg[6]->GetXaxis()->GetXmin()));
	// fraction
	float fraction = clus*100./IDneutrg[6]->GetEntries();

	sprintf(txt, "Before gate start: %d", clus);
	text->AddText(txt);
	sprintf(txt, "Fraction: %.2f%%",fraction);
	text->AddText(txt);
	
	if(fraction > 0.1){
		sprintf(txt,"=> ABOVE THRESHOLD!");
		text->SetTextColor(kRed);
	}
	else{
		sprintf(txt,"=> OK");
		text->SetTextColor(kGreen+3); // dark green
	}

	text->AddText(txt);
	text->SetFillColor(kWhite);
	text->SetTextSize(0.05);
	text->SetBorderSize(0);
	text->Draw();
	/// -------------


	////------- pad 5
	c[4]->cd(5);
	IDneutrg[5]->SetXTitle("Time difference [ms]");
	IDneutrg[5]->Draw(); 
	c[4]->GetPad(5)->SetLogx();
	
	// fit
	fit = new TF1("fit", "[0]*exp(-x*[1])", IDneutrg[5]->GetXaxis()->GetXmin(), IDneutrg[5]->GetXaxis()->GetXmax());
	fit->SetLineColor(kRed); 
	fit->SetLineWidth(2);      
	IDneutrg[5]->Fit(fit, "QEMR+");   
	
	// text on the plot
	float r = fit->GetParameter(1)*1000; // from ms^-1 to Hz
	sprintf(txt,"%.2f+-%.2f Bq", r ,fit->GetParError(1)*1000);
	text = new TPaveText(0.7,0.75,0.8,0.8, "NDC");
	text->AddText(txt);
	text->SetFillColor(kWhite);
	text->SetBorderSize(0);
	text->SetTextSize(0.045);
	text->Draw();
	
	/// -------------

	
	////------- pad 6: 
	c[4]->cd(6);
	
	IDneutrg[8]->SetXTitle("Time difference [us]"); 
	IDneutrg[8]->Draw(); 
	
	//**************************************************//


	//******************** canvas 05 ********************//
	// ID: pulse shape
	cout << "canvas 05" << endl;
		
	c[5] = new TCanvas("c5","[05] ID: pulse shape",ww,hh);
	c[5]->Divide(2,3);
	
	for(int i = 0; i < 6; i++){	
		c[5]->cd(i+1);
		IDtimes[i]->SetXTitle("ID decoded hit raw time (trigger time = 0) [ns]");
		IDtimes[i]->Draw();
		c[5]->GetPad(i+1)->SetLogy();
	}

	c[5]->Update();
	
	//**************************************************//


	//******************** canvas 06 ********************//
	//ID: event size VS channels
	cout << "canvas 06" << endl;
		
	c[6] = new TCanvas("c6","[06] ID: event size VS channels",ww,hh);
	c[6]->Divide(2,2);
	
	for(int i = 0; i < 4; i++){	
		c[6]->cd(i+1);
		IDchannels[i]->SetXTitle("Logical channel");
		IDchannels[i]->SetYTitle("Event size in ID (N decoded hits)");
		IDchannels[i]->SetTitleOffset(1.1,"Y");
		IDchannels[i]->Draw();
	}
	
	//**************************************************//

	
	//******************** canvas 07 ********************//
  // ID: hits in neutrino trigger, diff
	cout << "canvas 07" << endl;
	// note: uses the same histos as canvas 4
	
	c[7] = new TCanvas("c7","[07] ID: hits in neutrinos (TT1&BTB0)",ww,hh);	
	c[7]->Divide(2,2);

	IDneutrg[1]->SetLineColor(kBlue+2);

//	TPaveStats* st;
	
	for(int i = 0; i < 4; i++){
		c[7]->cd(i+1);
		c[7]->GetPad(i+1)->SetLogy();
		
		// scale if the run being validated was truncated
		IDc7sample[i]->Scale(bx->IDneutrg[i]->GetEntries()*1./bx->IDc7sample[i]->GetEntries()); 
		IDc7sample[i]->SetLineColor(kGreen+1); // darker green
		IDc7sample[i]->SetXTitle("Number of hits in ID");
		IDc7sample[i]->SetAxisColor(kBlack); 
		IDc7sample[i]->GetXaxis()->SetTitleColor(kBlack); 
		IDc7sample[i]->SetAxisColor(kBlack, "Y");
		IDc7sample[i]->SetLabelColor(kBlack);  
		IDc7sample[i]->SetLabelColor(kBlack, "Y");  
		IDc7sample[i]->SetStats(0);
		
		IDc7sample[i]->Draw();
		IDneutrg[i]->Draw("sames");

		// label
		lbl = "Run029114 2017 Jul 16";
		label = new TPaveLabel(0.6,0.6,0.7,0.7,lbl.c_str(),"NDC");
		label->SetTextColor(kGreen+1);
		label->SetBorderSize(0);
		label->SetFillColor(kWhite);
		label->SetTextSize(0.5);
		label->Draw("same");
	}
		
	IDc7sample[3]->SetXTitle("Charge in ID");
	
	//bx->c[7]->Modified();
	c[7]->Update();
	
	//**************************************************//


	//******************** canvas 08 ********************//
	// OD muon flag stability
	cout << "canvas 08" << endl;
		
	c[8] = new TCanvas("c8","[08] OD: muon flag stability",ww,hh);
	c[8]->Divide(2,1);
	
	////////// calculate all kinds of stuff //////////

	// determine maxima for the axes
	float yMax[2]; float binmax;
	yMax[0] = yMax[1] = 0;
	
	for(int i = 0; i < 6; i++){
		binmax = ODmuonFlag[i]->GetMaximum();
		if(binmax > yMax[i/3]){ yMax[i/3] = binmax; }
	}

	// NOTE: assuming the histos are the same size
	float xMax = ODmuonFlag[0]->GetNbinsX();

	////////// draw additional stuff //////////
	float c8mean[2]; c8mean[0] = 4250; c8mean[1] = 4000;
	float c8sigma;

	sprintf(txt,"Time in 10 bins (total %.2f hours)",(end_time-start_time)/60./60.);

	for(int i = 0; i < 2; i++){
		// calculate sigma
		c8sigma = c8mean[i] / sqrt(muonTT1count/10.);
	
		// draw frames to limit axes
		c[8]->cd(i+1);
		frame = c[8]->GetPad(i+1)->DrawFrame(0,0,xMax,yMax[i]*1.1);
		// titles
		if(i == 0){frame->SetTitle("Internal muons");}
		else{frame->SetTitle("External muons (TT2)");}

		//frame->SetXTitle("Time in 10 bins");
		frame->SetXTitle(txt);
		frame->SetYTitle("Events/day");	
		frame->SetTitleOffset(1.1,"Y");
		
		// draw mean lines
		sigma = new TLine(0,c8mean[i],xMax,c8mean[i]);
		sigma->SetLineColor(33); // light grey
		sigma->SetLineWidth(2);
		sigma->Draw("L, same");
		
		// draw sigma lines
		for(int j = 1; j < 3; j++){
			// + j*s
			sigma = new TLine(0, c8mean[i] + j*c8sigma, xMax, c8mean[i] + j*c8sigma);
			sigma->SetLineWidth(2);
			sigma->SetLineStyle(kDashed);
			sigma->SetLineColor(33); // light grey
			sigma->Draw("L, same");
			// - j*s
			sigma = new TLine(0, c8mean[i] - j*c8sigma, xMax, c8mean[i] - j*c8sigma);
			sigma->SetLineWidth(2);
			sigma->SetLineStyle(kDashed);
			sigma->SetLineColor(33); // light blue?
			sigma->Draw("L, same");
		}
	}

	// legends for later
	TLegend *legend[2];
	legend[0]  = new TLegend(.75,.2,.975,.4);
	legend[1]  = new TLegend(.6,.4,.975,.6);

	////////// draw the actual histos//////////
	for(int i = 0; i < 6; i++){
		c[8]->cd(i/3+1);
		ODmuonFlag[i]->SetLineColor(3-i%3);
		ODmuonFlag[i]->SetMarkerColor(3-i%3);
		ODmuonFlag[i]->Draw("P0,L,same");
		legend[i/3]->AddEntry( ODmuonFlag[i], ODmuonFlag[i]->GetTitle());
	}
	
	ODmuonFlag[1]->SetLineStyle(7); // long dashed line, cause MTB and MCR overlap
	ODmuonFlag[1]->SetLineWidth(2);
	ODmuonFlag[4]->SetLineStyle(7);
	ODmuonFlag[4]->SetLineWidth(2);

	// draw legends
	for(int i = 0; i < 2; i++){
		c[8]->cd(i+1);
		legend[i]->SetTextSize(0.04);
		legend[i]->Draw();
	}

	c[8]->Update();
	c[8]->Modified();
	
	//**************************************************//

	
	//******************** canvas 09 ********************//
	// OD: event size VS evnum 
	cout << "canvas 09" << endl;
		
	c[9] = new TCanvas("c9","[09] OD: event size VS evnum",ww,hh);
	c[9]->Divide(3,2); 
	
	// pad 2: pulser y-axis limit
	ODevsize[1]->SetAxisRange(pulserODylim[0]*0.8, pulserODylim[1]*1.2, "Y");
	
	for(int i = 0; i < 6; i++){
		c[9]->cd(i+1);
		ODevsize[i]->SetXTitle("Event number");
		ODevsize[i]->SetYTitle("Event size in OD (N decoded hits)");

		ODevsize[i]->Draw("simple");
		// fit
		fit = new TF1("line", "pol1", 0., ODevsize[i]->GetSize());
		ODevsize[i]->Fit("line","q");
	}
	
	//**************************************************//

	
	//******************** canvas 10 ********************//
	// OD: event size VS evnum (200 bins) 
	cout << "canvas 10" << endl;
		
	c[10] = new TCanvas("c10","[10] OD: event size VS evnum (200 bins)",ww,hh);
	c[10]->Divide(3,2); 
	
	// x-axis title
	sprintf(txt,"Ev. number in 200 bins (1 bin = %d events)", binstep);
	
	for(int i = 0; i < 6; i++){
		c[10]->cd(i+1);
		ODevsize200[i]->SetXTitle(txt);
		ODevsize200[i]->SetYTitle("Avg event size in OD (N decoded hits)");
		ODevsize200[i]->Draw("simple");
		// fit
		fit = new TF1("line", "pol1", 0., ODevsize200[i]->GetSize());
		ODevsize200[i]->Fit("line","q");
	}
	
	//**************************************************//
	
	
	//******************** canvas 11 ********************//
	// OD: muon trigger rate VS time
	cout << "canvas 11" << endl;
		
	c[11] = new TCanvas("c11","[11] OD: trigger rate",ww,hh);
  c[11]->Divide(2,2); 
	
	for(int i = 0; i < 4; i++){
		c[11]->cd(i+1);
		ODtrgRate[i]->SetAxisRange(0, end_time - start_time);
		ODtrgRate[i]->SetXTitle("Event time (1 bin = 10 min) [s]");
		ODtrgRate[i]->SetYTitle("Muon trigger rate (per day)");
		ODtrgRate[i]->SetTitleOffset(1.1,"Y");
		ODtrgRate[i]->Draw("simple");
		fit = new TF1("line", "pol1", 0., ODtrgRate[i]->GetSize());
		ODtrgRate[i]->Fit("line","q");
	}
	
	//**************************************************//
	
	
	//******************** canvas 12 ********************//
	// ID: hits in muon triggers
	cout << "canvas 12" << endl;
		
	c[12] = new TCanvas("c12","[12] ID: hits in muon triggers",ww,hh);
	c[12]->Divide(2,2);
	
	for(int i = 0; i < 3; i++){
		c[12]->cd(i+1);
		IDhitsMuon[i]->SetXTitle("Number of hits in ID");
		IDhitsMuon[i]->Draw();
		c[12]->GetPad(i+1)->SetLogy();
	}
	
	// pad 4
	c[12]->cd(4);
	c[12]->GetPad(4)->SetLogy();
	// draw TT2 first to have correct axis ranges
	IDhitsMuon[3]->SetXTitle("Number of clusters in ID");
	IDhitsMuon[3]->SetTitle("Clusters in muons");
	IDhitsMuon[3]->SetLineColor(kRed);
	IDhitsMuon[3]->Draw();
	IDhitsMuon[4]->Draw("same");
	
	// labels
	lbl = "External muons (TT2)";
	label = new TPaveLabel(0.65,0.75,0.75,0.85,lbl.c_str(),"NDC");
	label->SetBorderSize(0);
	label->SetTextColor(kRed);
	label->SetFillColor(kWhite);
	label->SetTextSize(0.4);
	label->Draw();
	
	lbl = "Internal muons (TT1&BTB4)";
	label1 = new TPaveLabel(0.65,0.65,0.76,0.75,lbl.c_str(),"NDC");
	label1->SetBorderSize(0);
	label1->SetFillColor(kWhite);
	label1->SetTextColor(kBlue);
	label1->SetTextSize(0.4);
	label1->Draw();
	
	c[12]->Update();
	
	//**************************************************//
	
	
	//******************** canvas 13 ********************//
	// OD: hits in muon triggers
	cout << "canvas 13" << endl;
		
	c[13] = new TCanvas("c13","[13] OD: hits in muon triggers",ww,hh);
	c[13]->Divide(2,2);
	
	for(int i = 0; i < 4; i++){
		c[13]->cd(i+1);
		ODhitsTT2[i]->SetLineColor(kRed);
		ODhitsTT2[i]->SetXTitle("Number of hits in OD");
		ODhitsTT2[i]->Draw();
		ODhitsTT1[i]->Draw("same");
		c[13]->GetPad(i+1)->SetLogy();

		//labels: reuse from c12
		label->Draw("same");
		label1->Draw("same");
	}

	ODhitsTT2[3]->SetXTitle("Number of clusters in OD");
	
	c[13]->Update();
	
	//**************************************************//
	
	
	//******************** canvas 14 ********************//
	// OD: pulse shape
	cout << "canvas 14" << endl;
		
	c[14] = new TCanvas("c14","[14] OD: pulse shape",ww,hh);
	c[14]->Divide(2,3);
	
	for(int i = 0; i < 6; i++){
		c[14]->cd(i+1);
		ODtimes[i]->Draw();
		ODtimes[i]->SetXTitle("OD decoded hit raw time (trigger time = 0) [ns]");
		c[14]->GetPad(i+1)->SetLogy();
	}

	c[14]->Update();
	
	//**************************************************//
	

	//******************** canvas 15 ********************//
	// OD: event size VS channels
	cout << "canvas 15" << endl;
		
	c[15] = new TCanvas("c15", "[15] OD: event size VS channels",ww,hh);
	c[15]->Divide(2,3);
	
	for(int i = 0; i < 6; i++){
		c[15]->cd(i+1);
		ODchannels[i]->SetXTitle("Logical channel");
		ODchannels[i]->SetYTitle("Event size in OD (N decoded hits)");
		ODchannels[i]->Draw();
		c[15]->GetPad(i+1)->SetLogy();
	}

	// pad 3: pulser
	// otherwise the numbers on the y-axis are not displayed for some reason just for this one histo:
	ODchannels[2]->SetMinimum(1);
	ODchannels[2]->SetMaximum( pow(10, log10(ODchannels[2]->GetMaximum())+1) );

	c[15]->Update();
	
	//**************************************************//


	//******************** canvas 16 ********************//
	// ID: hits in neutron trigger 

	cout << "canvas 16" << endl;
		
	c[16] = new TCanvas("c16","[16] ID: neutron (TT128)",ww,hh);
	c[16]->Divide(3,2);

	// pad 1: raw and decoded
	c[16]->cd(1);
	IDtrg128[0]->SetXTitle("Number of hits in ID");
	IDtrg128[0]->SetTitle("Raw and decoded hits");
	IDtrg128[0]->Draw();

	IDtrg128[1]->SetLineColor(kRed);
	IDtrg128[1]->Draw("same");
		
	// labels
	lbl = "Raw hits";
	label = new TPaveLabel(0.6,0.75,0.7,0.85,lbl.c_str(),"NDC");
	label->SetBorderSize(0);
	label->SetTextColor(kBlue+3);
	label->SetFillColor(kWhite);
	label->SetTextSize(0.4);
	label->Draw();
	
	lbl = "Decoded hits";
	label = new TPaveLabel(0.6,0.65,0.7,0.75,lbl.c_str(),"NDC");
	label->SetBorderSize(0);
	label->SetTextColor(kRed);
	label->SetFillColor(kWhite); 
	label->SetTextSize(0.4);
	label->Draw();

	// pad 2: hits in clusters
	c[16]->cd(2);
	c[16]->GetPad(2)->SetLogy();
	IDtrg128[2]->SetTitle("Hits in clusters");
	IDtrg128[2]->SetXTitle("Number of hits in ID");
	IDtrg128[2]->SetAxisRange(0,100,"X");
	IDtrg128[2]->Draw();

	IDtrg128[3]->SetLineColor(kRed);
	IDtrg128[3]->Draw("same");
	
	// labels
	lbl = "1-cluster events";
	label = new TPaveLabel(0.6,0.65,0.7,0.75,lbl.c_str(),"NDC");
	label->SetBorderSize(0);
	label->SetTextColor(kBlue+3);
	label->SetFillColor(kWhite);
	label->SetTextSize(0.5);
	label->Draw();

	lbl = "Multi-cluster events";
	label = new TPaveLabel(0.6,0.55,0.7,0.65,lbl.c_str(),"NDC");
	label->SetBorderSize(0);
	label->SetTextColor(kRed);
	label->SetFillColor(kWhite); 
	label->SetTextSize(0.5);
	label->Draw();

	// pad 3: charge in each cluster
	c[16]->cd(3);
	IDtrg128[4]->SetXTitle("Charge in ID");
	IDtrg128[4]->SetAxisRange(0,100,"X");
	IDtrg128[4]->Draw();
	c[16]->GetPad(3)->SetLogy();

	// pad 4: cluster start time
	c[16]->cd(4);
	IDtrg128[5]->SetXTitle("Cluster start time (trigger time = 0) [us]");
	IDtrg128[5]->Draw();
	c[16]->GetPad(4)->SetLogy();

	// pad 5: clusters
	c[16]->cd(5);
	IDtrg128[6]->SetXTitle("Number of clusters in ID");
	IDtrg128[6]->Draw();
	c[16]->GetPad(5)->SetLogy();

	// pad 6: empty boards
	c[16]->cd(6);
	IDtrg128[7]->SetXTitle("Number of boards");
	IDtrg128[7]->Draw();
	c[16]->GetPad(6)->SetLogy();

	c[16]->Update();

	//**************************************************//


	//******************** canvas 17 ********************//
	
	cout << "canvas 17" << endl;
	
	c[17] = new TCanvas("c17","[17] ID: neutron (TT128) pulse shape",1000,600);	
	c[17]->Divide(1,2);

	IDtrg128[9]->SetStats(0);
	IDtrg128[10]->SetStats(0);
	
	c[17]->cd(1);
	c[17]->GetPad(1)->SetLogy();
	IDtrg128[8]->SetXTitle("ID decoded hit raw time (trigger time = 0) [us]");
	IDtrg128[8]->Draw();

	c[17]->cd(2);
	c[17]->GetPad(2)->SetLogy();
	
	// muon
	IDtrg128[10]->SetXTitle("ID decoded hit raw time (muon trigger time = 0) [ns]");
	IDtrg128[10]->Draw();
	
	// neutron
	IDtrg128[9]->SetLineColor(kRed);
	IDtrg128[9]->Draw("same");

	// labels
	
	lbl = "Internal muon (TT1&BTB4)";
	label = new TPaveLabel(0.2,0.75,0.4,0.85,lbl.c_str(),"NDC");
	label->SetBorderSize(0);
	label->SetTextColor(kBlue+3);
	label->SetFillColor(kWhite); 
	label->SetTextSize(0.4);
	label->Draw();
	
	lbl = "Neutron (TT128)";
	label = new TPaveLabel(0.3,0.65,0.5,0.75,lbl.c_str(),"NDC");
	label->SetBorderSize(0);
	label->SetTextColor(kRed);
	label->SetFillColor(kWhite); 
	label->SetTextSize(0.4);
	label->Draw();

	// zoom
	TPad* p = new TPad("p","",.5,.6,0.95,0.95);
	p->Draw();
	p->SetFillStyle(4000);
	p->cd();
	p->SetLogy();
	frame = p->DrawFrame(0,0.0001,500,IDtrg128[10]->GetMaximum()*1.1);
	frame->SetLabelSize(0.1,"X");
	frame->SetLabelSize(0.1,"Y");
	IDtrg128[9]->Draw("same");
	IDtrg128[10]->Draw("same");
	
	//**************************************************//
		
	
}

// ------------------------------------------------------------

void BxValidation::OuterDetectorCheck(){
	LastMcrEvent = end_evnum;
	
	TTreeResult* res1 = (TTreeResult*)t->Query("evnum", "!((enabled_crates>>14)&1)");
	
	// check if mcr died at some point
	if( res1->GetRowCount() ){
		TSQLRow* row = res1->Next();
		LastMcrEvent = atoi( row->GetField(0) ) - 1;
	}
	// check if mcr was disabled from profile and still looks enabled
	else{
		TTreeResult* res2 = (TTreeResult*)t->Query("evnum", "muon.n_decoded_hits>0 && trigger.trgtype==32");
		if( !res2->GetRowCount() ){ LastMcrEvent = -2;}
		res2->Delete();
	}
	res1->Delete();

	// case where electronics was working fine but PMTs were off or had problems
	if(LastMcrEvent > 0){
		TTreeResult* res3 = (TTreeResult*)t->Query("evnum", "trigger.trgtype!=32", "", LastMcrEvent);
		TTreeResult* res4 = (TTreeResult*)t->Query("evnum", "muon.n_decoded_hits==0 && trigger.trgtype!=32", "", LastMcrEvent);
		if( res4->GetRowCount() > .4*res3->GetRowCount() ){LastMcrEvent = -3;}
		res3->Delete();
		res4->Delete();
	}

	LastAlignedEvent = LastMcrEvent;
	
	if(LastAlignedEvent > 0){
		TTreeResult* res3 = (TTreeResult*)t->Query("evnum", "trigger.trgtype==32 && muon.n_decoded_hits<180", "", LastMcrEvent);
		
		if(res3->GetRowCount() >= 5){
			TSQLRow* row = (TSQLRow*)res3->GetRows()->First();
			int first_disaligned_calib = atoi(row->GetField(0));
			TTreeResult* res4 = (TTreeResult*)t->Query("evnum", "trigger.trgtype==32 && muon.n_decoded_hits>180", "", first_disaligned_calib);
			
			if( res4->GetRowCount() ){
				TSQLRow* row = (TSQLRow*)res4->GetRows()->Last();
				LastAlignedEvent = atoi( row->GetField(0) );
			}
			else{ LastAlignedEvent = 0; }
			
			res4->Delete();
		}
		res3->Delete();
	}

	cout << endl << "CHECK OF THE OUTER DETECTOR:" << endl << endl;  

	// check MCR
	if(LastMcrEvent == end_evnum){ cout << "MCR was working fine for the whole run." << endl;}

	else if(LastMcrEvent > 0){
		cerr << "!! MCR was down after event " << LastMcrEvent << " (of " << end_evnum << " events, " << int(LastMcrEvent*100./end_evnum) << "%) !!" << endl;
	}
	else if(LastMcrEvent == -2){ cerr << "!! MCR was disabled from profile !!" << endl; }
	else if(LastMcrEvent == -3){ cerr << "!! Muon PMTs had problems !!" << endl; }

	// check ID-OD
	if(LastAlignedEvent > 0.99 * LastMcrEvent){ cout << "No problems with ID/OD-disalignments." << endl;}
	else if(LastMcrEvent >= 0){ cerr << "!! Disalignment of ID/OD events, beginning from event " << LastAlignedEvent + 1 << " !!" << endl; }

	cout << endl;

	// check trigger behaviour
	if(TT2count == 0){ cerr << "!! No trigger type 2. Check if MTB is enabled !!" << endl; }

}

// ------------------------------------------------------------

void BxValidation::ElectronicsRescale(){
/* rescale electronics calibrations histograms */
		
		// ************ CANVAS NO 3
		string title = "[3] Run" + runnum + " Raw hits VS channels norm. to # of triggers";
		
		raw_lg_renorm = new TCanvas("raw_lg_renorm", title.c_str());	
		raw_lg_renorm->Divide(1,2);
	
		raw_lg_renorm->cd(1);
		raw_lg_pulser->Scale( 1.0 / (float)bx->TT32count);
		raw_lg_pulser->SetXTitle("Logical channel");
		raw_lg_pulser->Draw ();
		
		raw_lg_renorm->cd(2);
		raw_lg_neutrino->Scale( 1.0 / (float)bx->TT1count);
		raw_lg_neutrino->SetXTitle("Logical channel");
		raw_lg_neutrino->Draw(); 
		
		raw_lg_renorm->Modified();
		raw_lg_renorm->Update();
}

// ------------------------------------------------------------

bool BxValidation::SaveValidation(){

	if (run_int < 5000 || run_int > 5e6) {
		cerr << "Invalid run number: " << bx->run_int << endl;
		cerr << "Did you run RunValidation() before this?" << endl;
		return false;
	}

	cout << endl;
	cout << "#######################################" << endl << endl;
	cout << "Run Validation saving procedure starting for run: " << bx->runnum << endl << endl;
	cout << "#######################################" << endl << endl;

	/// **************************** ///
	/// COLLECT THE DATA TO BE SAVED ///
	/// **************************** ///

	// part of it is saved in the Fill Histos function
		
	// path to where the root file should be moved
	// e.g. "http://bxmaster-data.lngs.infn.it//bxstorage/rootfiles/cycle_18/2017/Jul_16/Run029114_c18.root"
	loc_root_files = GetPath ( run_int, atoi(cycle.c_str()) ) + "Run" + runnum + "_c" + cycle + ".root";

	RunEvents = (int) RunEventsBxdb(run_int); // int // ?? why not num entries????
		
	// from canvas # 03: ID: neutrino trigger rate VS time
	VerticalMean1 = VerticalMean(IDtrgRate[0], 3., true); // all neutrino

	// from canvas # 01: ID: event size VS evnum
	VerticalMean2 = VerticalMean(IDevsize[0], 3., true); // LASER (TT8)
	VerticalMean3 = VerticalMean(IDevsize[2], 3., true); // RANDOM (TT64)
	VerticalMean4 = VerticalMean(IDevsize[1], 3., true); // PULSER (TT32)
	VerticalMean5 = VerticalMean(IDevsize[3], 3., true); // NEUTRINO (TT1) --> now TT1&BTB0

	// from canvas # 11: OD: muon trigger rate VS time
	MuonRateID = VerticalMean(ODtrgRate[0], 3., true); // muon TT1 all
	MuonRateOD25 = VerticalMean(ODtrgRate[3], 3., true); // muon TT2 ev.size > 25
		
	// from the part above called "CHECK FOR AMAZING OD BEHAVIOUR"

	// format "2017-07-16 17:44:22"
	loc_start_time = StartDate(run_int); 
	loc_validation_time = TimeNow(); // format
	
	loc_valid = "true";
	loc_group = week.c_str();

	///// comment /////
	
	string comment;

	
	cout << "Please enter your run comment (maximum 256 characters, everything else will be cut off). If everything is fine, just enter \"ok\". If you want to quit, enter -100: ";
	
	// clean the input to record new comment
	cin.clear(); // clean the fail flag
	if(clean_cin){cin.ignore(256, '\n');} // cleans the input only if run after RV; SV after SV will not glitch this way
	
	// this way it truncates on the first space and saves only the first word
//	cin >> comment;
	// this way we save the whole sentenece
	getline(cin, comment);
	comment = comment.substr(0,256);

	cout << "Your comment: " << comment << endl;
	if( atoi(comment.c_str()) == -100){ clean_cin = false; return false; }

	//****************** WRITE HISTOS  ****************//
	
	// the root file in which all the histograms which are drawn will be saved: run_plots/xxxxxx.root (run number)
	char output[100];
	// the run number in the name as to be without zero (why?). the easiet way to do it is convert to integer which automatically removes all the zeros

	sprintf(output, "/home/production/run_plots/%d.root", bx->run_int);
	// LOCAL
//	sprintf(output, "/home/production/Echidna_c18_mariia/run_plots/%d.root", bx->run_int);

	cout << endl << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	cout << "...... writing the histograms ......"<< endl << endl;
	
	TFile runplots(output, "recreate");

	// calibration and precalibration

	// canvas no 0 (precalib)
	runplots.mkdir("precalibration/");
	runplots.cd("precalibration/");

	cout << "~~~ precalibration ~~~" << endl;
	charge_vs_channel->Write();
	cout << charge_vs_channel->GetName() << endl;
	time_vs_channel->Write(); 
	cout << time_vs_channel->GetName() << endl;
	time_vs_channel_good->Write();
	cout << time_vs_channel_good->GetName() << endl;

	// calibration
	runplots.mkdir("calibration/");
	runplots.cd("calibration/");
	
	cout << "~~~ calibration ~~~" << endl;

	// canvas no 1
	raw_lg_pulser->Write();
	cout << raw_lg_pulser->GetName() << endl;
	raw_lg_neutrino->Write();
	cout << raw_lg_neutrino->GetName() << endl;

	// canvas no 2
	int lgnum;
	for(int i = 0; i < 14; i++){
		lgnum = (int) trigref->GetBinContent(i);
		if(lgnum){
			trgr[i]->Write();
			cout << trgr[i]->GetName() << endl;
		}
	}

	// canvas no 3
	for(int i = 0; i < 6; i++){
		lgnum = (int) lasref->GetBinContent(i);
		if(lgnum){
			las[i]->Write();
			cout << las[i]->GetName() << endl;
		}
	}

	// validation canvases
	runplots.cd();
	cout << "~~~ validation ~~~" << endl;

	// canvas # 01
	for(int i = 0; i < 6; i++){
		IDevsize[i]->Write();
		cout << IDevsize[i]->GetName() << endl;
	} 
	// canvas # 02
	for(int i = 0; i < 6; i++){
		IDevsize200[i]->Write();
		cout << IDevsize200[i]->GetName() << endl;
	}
	// canvas # 03
	for(int i = 0; i < 2; i++){
		IDtrgRate[i]->Write();
		cout << IDtrgRate[i]->GetName() << endl;
	}
	// canvas # 04 & 07
	for(int i = 0; i < 9; i++){
		IDneutrg[i]->Write();
		cout << IDneutrg[i]->GetName() << endl;
	}
	// canvas # 05
	for(int i = 0; i < 6; i++){
		IDtimes[i]->Write();
		cout << IDtimes[i]->GetName() << endl;
	}
	// canvas # 06
	for(int i = 0; i < 4; i++){
		IDchannels[i]->Write();
		cout << IDchannels[i]->GetName() << endl;
	}
	// canvas # 08
	for(int i = 0; i < 6; i++){
		ODmuonFlag[i]->Write();
		cout << ODmuonFlag[i]->GetName() << endl;
	}
	// canvas # 09
	for(int i = 0; i < 6; i++){
		ODevsize[i]->Write();
		cout << ODevsize[i]->GetName() << endl;
	}
	// canvas # 10
	for(int i = 0; i < 6; i++){
		ODevsize200[i]->Write();
		cout << ODevsize200[i]->GetName() << endl;
	}
	// canvas # 11
	for(int i = 0; i < 4; i++){
		ODtrgRate[i]->Write();
		cout << ODtrgRate[i]->GetName() << endl;
	}
	// canvas # 12
	for(int i = 0; i < 5; i++){ 
		IDhitsMuon[i]->Write(); 
		cout << IDhitsMuon[i]->GetName() << endl;
	}
	// canvas # 13
	for(int i = 0; i < 4; i++){
		ODhitsTT2[i]->Write();
		cout << ODhitsTT2[i]->GetName() << endl;
		ODhitsTT1[i]->Write();
		cout << ODhitsTT1[i]->GetName() << endl;
	}
	// canvas # 14
	for(int i = 0; i < 6; i++){ 
		ODtimes[i]->Write();
		cout << ODtimes[i]->GetName() << endl;
	}
	// canvas # 15
	for(int i = 0; i < 6; i++){
		ODchannels[i]->Write();
		cout << ODchannels[i]->GetName() << endl;
	}
	// canvas # 16, canvas # 17
	for(int i = 0; i < 11; i++){
		IDtrg128[i]->Write();
		cout << IDtrg128[i]->GetName() << endl;
	}
	
	runplots.Close();

	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	cout << "Histograms saved to: " << endl;
	cout << output << endl;
	
	//***************************************************//
	//****************** WRITE CANVASES  ****************//
	//***************************************************//
	
	TImage *img = TImage::Create();

	char imgname[100];
	// format C01 rather than C02
	string zero = "0";

	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	cout << "...... writing the canvases ......" << endl << endl;
	
	for(int i = 1; i < 18; i++){
		// after 9 we don't need that extra 0
		if(i == 10){zero = "";}

		sprintf(imgname, "/home/production/run_plots/C%s%d_%d.png", zero.c_str(), i, run_int);
		// LOCAL
//		sprintf(imgname, "/home/production/Echidna_c18_mariia/run_plots/C%s%d_%d.png", zero.c_str(), i, run_int);

		cout << imgname << endl;
		// write as PNG
		img->FromPad(c[i]);
		img->WriteImage(imgname);
	}

	
	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	cout << "Saving info to: " << endl;
	cout << "/home/production/run_validation_out.txt" << endl;

	// ofstream using space or "\t" in between gives different format anyway
	FILE* fp = fopen("/home/production/run_validation_out.txt","a+");

	// formatting wrong this way:
	//ofstream fp("/home/production/run_validation_out.txt",ofstream::app);

	// LOCAL
//	FILE* fp = fopen("/home/production/Echidna_c18_mariia/run_plots/run_validation_out.txt","a+");


	fprintf(fp, "%5d %3d %7d %7d %7d %6d %6d %8.1f %6.0f %6.0f %7d %6d %6d %6.0f %6.0f %7.0f %6.0f %7d %7d %20s %20s %5s %20s %s     %s\n",
			run_int,
			cycle_int,
			RunEvents,
			start_evnum,
			end_evnum,
			t_diff,
			btb_thresh,
			VerticalMean1,
			MuonRateID,
			MuonRateOD25,
			good_laben_channels,
			good_muon_channels,
			good_fadc_channels,
			VerticalMean2,
			VerticalMean3,
			VerticalMean4,
			VerticalMean5,
			LastMcrEvent,
			LastAlignedEvent,
			loc_start_time.c_str(),
			loc_validation_time.c_str(),
			loc_valid.c_str(),
			loc_group.c_str(),
			loc_root_files.c_str(),
			comment.c_str());   

	fclose(fp);

	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	cout << "Run Stability:" << endl;

  run_stability->Valid();
//  STAB
	
	cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	cout << "Run Validation saving procedure ended." << endl;
	cout << "#######################################" << endl << endl;

	return true;
}


//////////////////////////////////////////////////////////
/// **************** HELPER FUNCTIONS **************** ///
//////////////////////////////////////////////////////////

long RunEventsBxdb(int RunNumber) {
	/*
	 * INPUT: run number (without zeros)
	 * OUTPUT: number of events in this run
	 * EXAMPLE:
	 * root [0] RunEvents(29114)
	 * (long)439465
	 */
	
	TSQLServer* db = TSQLServer::Connect("pgsql://bxdb.lngs.infn.it/daq_config","borex_guest","xyz");
	if (db) {
		TSQLResult* result = db->Query(Form("SELECT \"Events\" FROM \"Run\" WHERE \"RunNumber\" = %d", RunNumber) );
		if (result){
			if (result->GetRowCount() == 1) {
				TSQLRow* row = result->Next(); 
				return atol(row->GetField(0));
			}
		}
		db->Close();
	}
	std::cerr << "Error connecting to DB in RunEvents()\n";
	exit(11);
}

// ------------------------------------------------------------


string StartDate(int RunNumber) {
	/*
	 * INPUT: run number
	 * OUTPUT: string containing date and time
	 * EXAMPLE:
	 * root [1] StartDate(29114)
	 * (class TString)"2017-07-16 17:44:22"
	 */
	
	TSQLServer* db = TSQLServer::Connect("pgsql://bxdb.lngs.infn.it/daq_config","borex_guest","xyz");
	if (!db){	
		std::cerr << "Error connecting to DB in StartDate()\n";
		exit(20);
	}
	TSQLResult* result = db->Query(Form("SELECT \"StartTime\" FROM \"Run\" WHERE \"RunNumber\" = %d", RunNumber) );
	if (!result){
		std::cerr << "Error querying DB in StartDate()\n";
		db->Close();
		exit(21);
	}
	if( result->GetRowCount() == 1) {
		TSQLRow* row = result->Next();
		std::ostringstream oss;
		oss << row->GetField(0);
		return (oss.str());
	}  else {
		std::cerr << "Error in quering DB in StartDate() \n";
		const std::string bla;;
		return bla;
	}
	db->Close();
}

// ------------------------------------------------------------

string GetYear(int runnum){
	/*
	 * INPUT: run number 
	 * OUTPUT: year
	 * EXAMPLE:
	 * root [1] GetYear(29114)
	 * (class TString)"2017"
	 */

	// string of the type "2017-07-16 17:44:22
	string start_date = StartDate(runnum);

	return start_date.substr(0,4);
}

// ------------------------------------------------------------

string TimeNow() {
	/*
	 * OUTPUT: gives current date and time
	 * EXAMPLE:
	 * root [2] TimeNow()
	 * (class TString)"2017-09-12 14:45:59"
	 */
	
	time_t t = time(0);
	static char s[200];
	struct tm *ts = localtime(&t);
	strftime(s,200,"%Y-%m-%d %H:%M:%S",ts);
	//        printf("Tempo: %s\n",s);
	const std::string s_time = s;
	return s_time;
}

// ------------------------------------------------------------

long int StartDate_int(int RunNumber) {
	/*
	 * INPUT: run number (integer)
	 * OUTPUT: number that contains info about the start date
	 
	 * EXAMPLE:
	 * root [0] StartDate_int(29114)
	 * (long)1500219862
	 
	 * This number is not human-friendly, but we use it in GetWeek(), and GetPath(), and some checks in the Validation code
	 */
	
	TSQLServer* db = TSQLServer::Connect("pgsql://bxdb.lngs.infn.it/daq_config","borex_guest","xyz");
	if (!db){
		std::cerr << "Error connecting to DB in StartDate_int()\n";
		exit(31);
	}
	TSQLResult* result = db->Query(Form("SELECT \"StartTime\" FROM \"Run\" WHERE \"RunNumber\" = %d", RunNumber) );
	if (!result){
		std::cerr << "Error getting Query in StartDate_int()\n";
		db->Close();
		exit(32);
	}
	if( result->GetRowCount() == 1) {
		TSQLRow* row = result->Next(); 
		result = db->Query(Form("SELECT EXTRACT(EPOCH FROM TIMESTAMP WITH TIME ZONE\'%s\')", row->GetField(0)));

		if (result->GetRowCount() == 1) {
			row = result->Next();
			return ::strtol (row->GetField(0), 0, 10);
		}
	}
	db->Close();
	return -1;
}

// ------------------------------------------------------------

string GetPath (int runnum, int cycle) {
	/*
	 * INPUT: run number (integer), cycle number
	 * OUTPUT: string that represents path on bxmaster where the file should be stored
	 * EXAMPLE:
	 * root [2] GetPath(29114, 19)
	 * (class TString)"http://bxmaster-data.lngs.infn.it//bxstorage/rootfiles/cycle_18/2017/Jul_16/"
	 */

	long int start_time = StartDate_int(runnum);

	struct tm *start_date = ::localtime (&start_time);
	//time_t week_second = start_time - ((start_date->tm_wday * 24 + start_date->tm_hour) * 60 + start_date->tm_min) * 60;
	time_t week_second = start_time - ((start_date->tm_wday * 24 + start_date->tm_hour - start_date->tm_isdst) * 60 + start_date->tm_min) * 60;
	struct tm *week_date = ::localtime (&week_second);
	char tmp_str[100];
	::strftime (tmp_str, 99, "/%Y/%b_%d/", week_date);
	std::ostringstream str;
	str << "http://bxmaster-data.lngs.infn.it//bxstorage/rootfiles/cycle_" << cycle << tmp_str;
	return (str.str ());
}

// ------------------------------------------------------------

string GetWeek (int runnum) {
	/*
	 * INPUT: run number (integer)
	 * OUTPUT: a string month_day representing the week of the shift 
	 * EXAMPLE:
	 * root [1] GetWeek(29114)
	 * (class TString)"Jul_16"
	 */


	long int start_time = StartDate_int(runnum);
	
	struct tm *start_date = ::localtime (&start_time);
	//time_t week_second = start_time - ((start_date->tm_wday * 24 + start_date->tm_hour) * 60 + start_date->tm_min) * 60;
	time_t week_second = start_time - ((start_date->tm_wday * 24 + start_date->tm_hour - start_date->tm_isdst) * 60 + start_date->tm_min) * 60;

	struct tm *week_date = ::localtime (&week_second);
	char tmp_str[100];
	::strftime (tmp_str, 99, "%b_%d", week_date);
	std::ostringstream str;
	str << tmp_str;
	return str.str ();
}

// ------------------------------------------------------------

void SetGeneralStyle() {

	gStyle->SetOptStat("ne"); 
//	gStyle->SetOptFit(0);
//	gStyle->SetOptStat(111111);
//	gStyle->SetOptStat(000000011);
//	gStyle->SetStatFontSize(0.1);

// Set position (fraction of pad size)
//	gStyle->SetStatY(0);                
//	gStyle->SetStatX(1);
	// Set width and height of stat-box (fraction of pad size)
	gStyle->SetStatW(0.23);                
	gStyle->SetStatH(0.17);                
	
	gStyle->SetTitleSize(0.055);
	gStyle->SetTitleSize(0.044,"X");
	gStyle->SetTitleSize(0.044,"Y");

	// bullet
	gStyle->SetMarkerStyle(7);    // small squared bullet
	gStyle->SetMarkerColor(50);   // red

}

/// Mathematical helper functions

// ------------------------------------------------------------

float VerticalMean(TH1F* h, float rms_limit, bool IgnoreZeros) {
	/* 
	 * compute mean of vertical axis in one-dim histos 
	 * through an option, you can include zero bins or not
	 * points that are out mean more than 5 sigmas are skipped in second iteration
	 * 
	 */

	int size = h->GetSize() - 2; //number of bins, without under/over flow
	float mean = 0.,rms=0;
	float S = 0.;
	int cnt=0;
	for(int i=1; i<=size; i++) {
		double value= h->GetBinContent (i); 
		if (value >0.0000001 || !IgnoreZeros) {
			cnt ++;
			double delta = value - mean;   
			if (cnt) mean += delta / cnt;
			S += delta * (value - mean);
		}
	}
	if(cnt) rms = TMath::Sqrt(double(S / cnt));

	//std::cout << "mean " << mean << " rms " << rms << std::endl;


	double mean_old = mean;
	mean = 0.;
	S = 0.;
	cnt=0;

	for(int i=1; i<=size; i++) {
		double value= h->GetBinContent (i); 
		//cout << "value " << value << " mean_old" << mean_old << " rms_limit " << rms_limit << " rms " << rms << endl; 
		if ( TMath::Abs(value - mean_old) < rms_limit* rms && (value >0.0000001 || !IgnoreZeros)) {
			cnt ++;
			double delta = value - mean;   
			if(cnt) mean += delta / cnt;
			S += delta * (value - mean);
		}
	}	
	if(cnt) rms = ::sqrt(S / cnt);
//	std::cout << "mean " << mean << " rms " << rms << " (truncated) " << std::endl;
	return mean;
}

// ------------------------------------------------------------

double SkewGaus(double *x, double *p){
	/* 
	 * fit function for cluster start time distribution asymmetric Gaussian
	 */
	
	double N = p[0];
	double mu = p[1];
	double sL = p[2];
	double sR = p[3];
	double sAvg = (sL + sR) / 2.; 


	if(x[0] <= mu) return 2. * N  / (TMath::Sqrt(2. * TMath::Pi()) * sAvg) * TMath::Exp(-0.5 * TMath::Power((x[0] - mu) / sL, 2));
	else return 2. * N  / (TMath::Sqrt(2. * TMath::Pi()) * sAvg) * TMath::Exp(-0.5 * TMath::Power((x[0] - mu) / sR, 2));
}

// ------------------------------------------------------------

bool FileCheck(string filepath){

	// check whether file exists
	TFile* frun = TFile::Open(filepath.c_str());
	if(!frun ){
//		cerr << "File " << filepath << " does not exist." << endl;
		return false;
	}

	// check whether file is Zombie
	else if( frun->IsZombie() ){
		cerr << "File " << filepath << " is a zombie." << endl;
		return false;
	}
	delete frun;
	return true;
}
