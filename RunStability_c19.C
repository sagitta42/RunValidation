#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <string>
#include <iterator>
#include <TROOT.h>
#include <TObjectTable.h>
#include <TNamed.h>
#include <TSystem.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TKey.h>
#include <TMemStat.h>

#define N_HISTOS_CANVAS 4
//#define SIZEOF_VSRUN 25200
#define SIZEOF_VSRUN 35200


//simple option list for stability_check class. Put here all histos of interest before processing stability class
std::vector<std::string> v_histo_draw;
std::vector<std::string> v_histo_style;
std::vector<std::string> v_histo_style_single;
std::vector<std::string> v_ytitles;

void init_style() {
  //clear global vectors with histogram styles
  v_histo_draw.clear();
  v_histo_style.clear();
  v_histo_style_single.clear();
  v_ytitles.clear();

	// simply add new histograms here if you want, no need to modify anything else at all in the code
  v_histo_draw.push_back(std::string("stability_internal_mu"));
  v_histo_style.push_back(std::string("E"));
  v_histo_style_single.push_back(std::string("E"));
	v_ytitles.push_back("Rate (counts/day)");
  
	v_histo_draw.push_back(std::string("stability_external_mu"));
  v_histo_style.push_back(std::string("E"));
  v_histo_style_single.push_back(std::string("E"));
	v_ytitles.push_back("Rate (counts/day)");
  
	v_histo_draw.push_back(std::string("nlive_vs_run"));
  v_histo_style.push_back(std::string("HIST P"));
  v_histo_style_single.push_back(std::string("E"));
	v_ytitles.push_back("Number of PMTs");
 
	v_histo_draw.push_back(std::string("stability_2"));
  v_histo_style.push_back(std::string("E"));
  v_histo_style_single.push_back(std::string("E"));
	v_ytitles.push_back("Rate (counts/day)");


  v_histo_draw.push_back(std::string("rate1_vs_run_IV"));
  v_histo_style.push_back(std::string("E"));
  v_histo_style_single.push_back(std::string("E"));
	v_ytitles.push_back("Rate (counts/day)");
  
	v_histo_draw.push_back(std::string("rate1_vs_run_FV"));
  v_histo_style.push_back(std::string("E"));
  v_histo_style_single.push_back(std::string("E"));
	v_ytitles.push_back("Rate (counts/day)");

  v_histo_draw.push_back(std::string("rate3_vs_run_IV"));
  v_histo_style.push_back(std::string("E"));
  v_histo_style_single.push_back(std::string("E"));
	v_ytitles.push_back("Rate (counts/day)");
  
	v_histo_draw.push_back(std::string("rate3_vs_run_FV"));
  v_histo_style.push_back(std::string("E"));
  v_histo_style_single.push_back(std::string("E"));
	v_ytitles.push_back("Rate (counts/day)");
  

	v_histo_draw.push_back(std::string("rate4_vs_run_IV"));
  v_histo_style.push_back(std::string("E"));
  v_histo_style_single.push_back(std::string("E"));
	v_ytitles.push_back("Rate (counts/day)");
  
	v_histo_draw.push_back(std::string("rate4_vs_run_FV"));
  v_histo_style.push_back(std::string("E"));
  v_histo_style_single.push_back(std::string("E"));
	v_ytitles.push_back("Rate (counts/day)");
  
  v_histo_draw.push_back(std::string("rate6_vs_run_IV"));
  v_histo_style.push_back(std::string("E"));
  v_histo_style_single.push_back(std::string("E"));
	v_ytitles.push_back("Rate (counts/day)");
  
	v_histo_draw.push_back(std::string("rate6_vs_run_FV"));
  v_histo_style.push_back(std::string("E"));
  v_histo_style_single.push_back(std::string("E"));
	v_ytitles.push_back("Rate (counts/day)");
  

	v_histo_draw.push_back(std::string("stability_radon_events"));
  v_histo_style.push_back(std::string("E"));
  v_histo_style_single.push_back(std::string("E"));
	v_ytitles.push_back("Rate (counts/day)");
  
	v_histo_draw.push_back(std::string("stability_radon_events_FV"));
  v_histo_style.push_back(std::string("E"));
  v_histo_style_single.push_back(std::string("E"));
	v_ytitles.push_back("Rate (counts/day)");

  v_histo_draw.push_back(std::string("startT_sRight"));
  v_histo_style.push_back(std::string("E"));
  v_histo_style_single.push_back(std::string("E"));
	v_ytitles.push_back("Sigma (ns)");
  
	v_histo_draw.push_back(std::string("startT_sLeft"));
  v_histo_style.push_back(std::string("E"));
  v_histo_style_single.push_back(std::string("E"));
	v_ytitles.push_back("Sigma (ns)");


}

//class definition, implementation at the EOF
class stability_check {
 private:
  stability_check() {};
  ~stability_check();
 public:
  stability_check(char* history_file, char* single_file, int);
  void draw_histos();
  void add_histos();
  void delete_single_file();
  void save_histos();
 private:
  //pointers to history and single stability root file
  TFile *history_file;
  TFile *single_file;
  //vector of all histograms in stability root file
  std::vector<std::string> v_histo_names;
  //vectors with histos to check & merge
  std::vector<TH1F*> v_history_histos;
  std::vector<TH1F*> v_single_histos;
  //canvases also handled
  std::vector<TCanvas*> v_canvases;
  //run_number under check
  int run_number;

};

class RunStability : public TNamed {
 public:
  RunStability(int run_number);
  void static Valid();
	void static UpdateBkgMonitor();
  void static Clean(char* remove_single_file = "no");
  int GetRun();
 private:
  char* history_filename;
  char* single_filename; 
  static char* lock_filename; 
  static stability_check* framework;
  static RunStability *me;
  int run_number;
  ~RunStability();
  ClassDef(RunStability,0)
};

RunStability::RunStability(int _run_number) : TNamed("Run_Stability", "Run_Stability object"), run_number(_run_number) {
  //check if RunStability exists already 
  if(me != 0) {
    throw std::runtime_error("Run_Stability object already in memory ! Please clean (RunStability::Clean()) or validate (RunStability::Valid()) stability!" ) ;
  }

  //check if stability.lock file exists ?
  std::ifstream in;
  in.open(Form("%s", lock_filename ), std::ifstream::in );
  if (in.good() ) {
    throw std::runtime_error("stability.lock file exists! Please remove and check if no one is using RunStability macro right now." ) ;
  }
  
	//create lock file
  gSystem->Exec(Form("touch %s", lock_filename ) );

  //init style of histogram 
	init_style();

  //OPTIONS, the only place to configure the script
  history_filename = "/home/production/stability_c19/bxfilter/stability_c19.root";
	// LOCAL
//  history_filename = "/home/production/Echidna_c18_mariia/stability_c18.root";
  single_filename  = Form("/home/production/Echidna/stability_%i_c19.root", run_number);
	// LOCAL
//  single_filename  = Form("/home/production/Echidna_c18_mariia/stability_%i_c18.root", run_number);

  //MAIN create framework -> load all histos from histo and single stability files
  framework = new stability_check(history_filename, single_filename, run_number);

  //present historical histos + new single run in separate canvases
  framework->draw_histos();

  //remember static pointer to Run_Stability singleton
  me = this;

  //add to object table in memory
  gObjectTable->AddObj(me);
}

///////////////////////////////////////////////////////////////// SOURCE CODE ///////////////////////////////////////////////////////////////////////////////

void RunStability::Valid() { 
  if(me == 0) {
    throw std::runtime_error("Run_Stability must be created, e.g. RunStability(22458)" ) ;
  }
	me->framework->add_histos();
  me->framework->save_histos();
  //clean memory
  RunStability::Clean("remove");
	//RunStability::UpdateBkgMonitor();
  //RunStability::Clean();

	// LOCAL
	gSystem->Exec("/home/production/Echidna/monitor_dir/bkg_mon/background");
	gSystem->Exec("/home/production/Echidna/monitor_dir/bkg_mon/homepage");
}

void RunStability::Clean(char* remove) {
  if(!gObjectTable->PtrIsValid(me)) return;
  if(me == 0) {
    throw std::runtime_error("Run_Stability must be created, e.g. RunStability(22458)" ) ;
  }
  if(std::string(remove) == "remove") framework->delete_single_file();
  delete me;
}

RunStability::~RunStability() {
  delete framework;
  gSystem->Exec(Form("rm -f %s", lock_filename ) );
  gObjectTable->Remove(me);
  me = 0;
}

//static pointers
RunStability* RunStability::me = 0;

stability_check* RunStability::framework = 0;

char* RunStability::lock_filename = "/home/production/Echidna/stability.lock";  
// LOCAL
//char* RunStability::lock_filename = "/home/production/Echidna_c18_mariia/stability.lock";  

int RunStability::GetRun() {
  return run_number;
}

void RunStability::UpdateBkgMonitor() {

}

//source code for stability_check class
stability_check::stability_check(char *history_filename, char *single_filename, int _run_number) : history_file(0), single_file(0), run_number(_run_number) {
  history_file = new TFile(history_filename, "UPDATE");
  single_file  = new TFile(single_filename,  "READ");

  if(! history_file->IsOpen() ) {
    throw std::runtime_error("empty pointer to history_file in stability macro !" );
  }
  if(! single_file->IsOpen() ) {
    throw std::runtime_error("empty pointer to single_file in stability macro ! probably stability root file was not produced." );
  }

  //get names of all histograms in stability root file and put to v_histo_names vector
  TList* list = history_file->GetListOfKeys();
  TIter next(list) ;
  TKey* key ;
  TObject* obj ;
  while ( (key = (TKey*) next()) ) {
    obj = key->ReadObj() ;
    if (obj->InheritsFrom("TH1F")) {
      std::string str(obj->GetName());
      str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
      v_histo_names.push_back(str);
    }
  }
  
	//load histos from history and single file
  for(std::vector<std::string>::const_iterator it = v_histo_names.begin(); it != v_histo_names.end(); it++) {
    //load to memory history histograms
    TH1F* hist = (TH1F*)history_file->Get(it->c_str());
    if(hist == 0) throw std::runtime_error(Form("null pointer to histogram %s in history file !", it->c_str()));
    v_history_histos.push_back(new TH1F(*hist));
		v_history_histos.back()->SetDirectory(0);
    delete hist;
    
		//single histograms
    hist = (TH1F*) single_file->Get(it->c_str());
    if(hist == 0) throw std::runtime_error(Form("null pointer to histogram %s in single file ! probably stability root file is corrupted. ", it->c_str()));
    //std::cout << "Entries in single histogram " << hist->GetName() << " : " << hist->GetEntries() << std::endl;
    v_single_histos.push_back(new TH1F(*hist));
		v_single_histos.back()->SetDirectory(0);
    delete hist;
  }

}

stability_check::~stability_check() {
  history_file->cd();
  //for(std::vector<TH1F*>::iterator    it = v_history_histos.begin(); it != v_history_histos.end(); it++) delete *it;
  //for(std::vector<TH1F*>::iterator    it = v_single_histos.begin() ; it != v_single_histos.end() ; it++) delete *it;
  for(std::vector<TCanvas*>::iterator it = v_canvases.begin()      ; it != v_canvases.end()      ; it++) delete *it;
  delete single_file;
  delete history_file;
}

void stability_check::draw_histos() {

	// canvas size
	int ww = 900;
	int hh = 600;

  history_file->cd();

	// canvas pad counter
	int pos = 5;

	// go through the list of histo names to be drawn
	int n_histos = v_histo_draw.size();
	for(int n = 0; n < n_histos; n++){
//		std::cout << v_histo_draw[n] << std::endl;

		// create new canvas if we're out of 4 pads
		if( pos > N_HISTOS_CANVAS ) {
      TCanvas* canvas = new TCanvas(Form("canvas_stab_%d_%d", run_number, v_canvases.size() + 1), Form("Stability check for run %d canvas n. %d", run_number,  v_canvases.size() + 1), ww, hh);
      canvas->Divide(2,2,0.001,0.001);
      v_canvases.push_back(canvas);
      pos = 1;
    }

//		std::cout << "pos " << pos << std::endl;
    
		// choose the last available canvas
		v_canvases.back()->cd( pos );

		// history histos
    TH1F* hist = (TH1F*)history_file->Get(v_histo_draw[n].c_str());
    if(hist == 0) throw std::runtime_error(Form("null pointer to histogram %s in history file !", v_histo_draw[n].c_str()));
//		std::cout << v_history_draw[n]->GetName() << std::endl;//<< " " << v_histo_style[n] << " " << v_history_draw[n]->GetMarkerSize() << std::endl;
		hist->SetDirectory(0);
		hist->Draw(v_histo_style[n].c_str());
    hist->GetXaxis()->UnZoom();
    hist->GetYaxis()->UnZoom();
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(0.6);
		hist->SetXTitle("Run number");
		hist->SetYTitle(v_ytitles[n].c_str());
		hist->SetTitleSize(0.055);
		hist->SetTitleSize(0.048,"X");
		hist->SetTitleSize(0.048,"Y");
		hist->GetXaxis()->SetRangeUser(run_number - 200, run_number + 80);

		// single file histos
    hist = (TH1F*)single_file->Get(v_histo_draw[n].c_str());
    if(hist == 0) throw std::runtime_error(Form("null pointer to histogram %s in history file !", v_histo_draw[n].c_str()));

		if(hist->GetEntries() > 0){
			hist->Draw((v_histo_style_single[n] + std::string(" same")).c_str());
			hist->SetMarkerStyle(20);			 
			hist->SetMarkerSize(1.1);
			hist->SetMarkerColor(kRed);
			hist->SetLineColor(kRed);
		}
    
		pos++;
	}// for loop on n histos

	//close empty canvas
  if(pos == 1) v_canvases.back()->Close();

}

void stability_check::add_histos() {
  if(v_canvases.size() < 1) throw std::runtime_error("Please draw all histograms before adding !" );
  
	int sz = v_histo_names.size();
	for(int n = 0; n < sz; n++){
    //adding histos general
    //v_history_histos[n]->Add(v_single_histos[n]);
    //stupid adding, works only if run_number - 5002 is the first bin and X axis is alwasy run_number class
		//std::cout << v_history_histos[n]->GetSize() << " SIZEOF_VSRUN " << SIZEOF_VSRUN << std::endl;
    
		if(v_history_histos[n]->GetSize() == SIZEOF_VSRUN) {
//			float cnt = v_single_histos[n]->GetBinContent(run_number - 5002);
//			float err = v_single_histos[n]->GetBinError  (run_number - 5002);
//			v_history_histos[n]->SetBinContent(run_number - 5002, cnt);
//			v_history_histos[n]->SetBinError(run_number - 5002, err);
      v_history_histos[n]->SetBinContent(run_number - 5002, v_single_histos[n]->GetBinContent(run_number - 5002));
      v_history_histos[n]->SetBinError  (run_number - 5002, v_single_histos[n]->GetBinError  (run_number - 5002));
//			if(v_histo_names[n].substr(0,6) == "events") {
			 // std::cout << v_histo_names[n] << "events: " << v_single_histos[n]->GetBinContent(run_number - 5002) << std::endl;
//		  }
    }
  }
}


void stability_check::save_histos() {
  //backup stability histo file
  gSystem->Exec(Form("cp -f %s %s", history_file->GetName(), Form("%s.bak", history_file->GetName() ) ) );
  
  history_file->cd();

  int n = 0;
  for(std::vector<std::string>::const_iterator it = v_histo_names.begin(); it != v_histo_names.end(); it++) {
    gDirectory->Delete(Form("%s;1",it->c_str()));
    gDirectory->Delete(Form("%s;2",it->c_str()));
    gDirectory->Delete(Form("%s;3",it->c_str()));
    v_history_histos[n]->Write();
    n++;
  }
}

void stability_check::delete_single_file() {
  //std::cout << Form("rm -f %s", single_file->GetName() ) << std::endl;
  gSystem->Exec(Form("rm -f %s", single_file->GetName() ) );
}
