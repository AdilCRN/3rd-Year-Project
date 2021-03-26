#define L1CaloForward_cxx
#include "L1CaloForward.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdlib.h> 
#include <iostream>
#include <TLorentzVector.h>
#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <TGraphAsymmErrors.h>
#include <valarray>

using namespace std;

int main(int argvc, char *argv[]) {

  // Read arguments

  cout<<"Start "<<endl;
  string fname = "datafile.root";
  if (argvc>1){
    fname = argv[1];
  }
  cout<<"Filename: "<<fname<<endl;

  string oname = "histos";
  if (argvc>2){
    oname = argv[2];
  }
  cout<<"Output name: "<<oname<<endl;

  int option = 0;
  if (argvc>3){
    option  = atoi(argv[3]);
  }
  cout<<"Option: "<<option<<endl;
  
  // Setup the process
  TROOT myRoot("L1CaloForward","L1CaloForward");

  // Tree name
  string tname("ntuple");

  // Define TChain and add input file
  TChain *feChain  = new TChain(tname.data());
  feChain->Add(fname.data());
  
  // Check number of entries
  Long64_t nEntries = Long64_t(feChain->GetEntries());
  cout<< "Number of Events: "<< nEntries<<endl;

  // Build analysis instance and run the Loop method 
  L1CaloForward* Ana = new L1CaloForward((TChain *) feChain, oname, option);
  cout<<"Start Loop "<<endl;
    Ana->Loop();
}




void L1CaloForward::Loop()
{

  double PI = 3.1415927;
  TCanvas *can1 = new TCanvas("Efficiency1","Histograms",600,600); 

  //strhisfile seems to create a path for the output file
  string strhisfile = string("./") + outf + string(".root");
  TFile hisfile(strhisfile.data(), "recreate");

  int debug = opt;

  cout<<"DEbuglevel: "<<debug<<", "<<"Outfile: "<<strhisfile<<endl<<endl;


  // Define histograms
  // L1 Calo matched to Truth
  //TH1F Class - 1-d Histogram with a float per channel
  TH1F *hL1EleEt   = new TH1F("hL1EleEt",  "hL1EleEt",  40,5.,100.);
  TH1F *hL1EleEta  = new TH1F("hL1EleEta", "hL1EleEta", 40,-5.,5.);
  TH1F *hL1ElePhi  = new TH1F("hL1ElePhi", "hL1ElePhi", 40,-5.,5.);
  //eFEX electrons
  TH1F *hL1EfexEleEt   = new TH1F("hL1EfexEleEt",  "hL1EfexEleEt",  40,5.,100.);
  TH1F *hL1EfexEleEta  = new TH1F("hL1EfexEleEta", "hL1EfexEleEta", 40,-5.,5.);
  TH1F *hL1EfexElePhi  = new TH1F("hL1EfexElePhi", "hL1EfexElePhi", 40,-5.,5.);
  // All Truth
  TH1F *hTrEleEt   = new TH1F("hTrEleEt",  "hTrEleEt",  40,5.,100.);
  TH1F *hTrEleEta  = new TH1F("hTrEleEta", "hTrEleEta", 40,-5.,5.);
  TH1F *hTrElePhi  = new TH1F("hTrElePhi", "hTrElePhi", 40,-5.,5.);
  // Truth matched to L1Calo
  TH1F *hTrmEleEt  = new TH1F("hTrmEleEt", "hTrmEleEt", 40,5.,100.);
  TH1F *hTrmEleEta  = new TH1F("hTrmEleEta", "hTrmEleEta", 40,-5.,5);
  TH1F *hTrmElePhi  = new TH1F("hTrmElePhi", "hTrmElePhi", 40,-5.,5.);

  // Truth matched to L1Calo eFEX
  TH1F *hTrmEfexEleEt  = new TH1F("hTrmEfexEleEt", "hTrmEfexEleEt", 40,5.,100.);
  TH1F *hTrmEfexEleEta  = new TH1F("hTrmEfexEleEta", "hTrmEfexEleEta", 40,-5.,5.);
  TH1F *hTrmEfexElePhi  = new TH1F("hTrmEfexElePhi", "hTrmEfexElePhi", 40,-5.,5.);

  // Hadronic Jet
//   TH1F *jEleEt   = new TH1F("hL1EleEtHad",  "hL1EleEtHad",  20,10.,100.);
//   TH1F *jEleEta  = new TH1F("hL1EleEtaHad", "hL1EleEtaHad", 20,-5.,5.);
//   TH1F *jElePhi  = new TH1F("hL1ElePhiHad", "hL1ElePhiHad", 20,-5.,5.);

  //
//   TH1F *hL1EleEt   = new TH1F("hjEleEt",  "hjEleEt",  20,10.,100.);
//   TH1F *hjEleEta  = new TH1F("hjEleEta", "hjEleEta", 20,-5.,5.);
//   TH1F *hjElePhi  = new TH1F("hjElePhi", "hjElePhi", 20,-5.,5.);

  //Jet rejection variables
  TH1F *hL1EleIso   = new TH1F("hL1EleIso",  "hL1EleIso",  40,-10000.,40000.);
  TH1F *hL1EleHad   = new TH1F("hL1EleHad",  "hL1EleHad",  40,-1.,1.);

  TH1F *hL1EleRange   = new TH1F("hL1EleRange",  "hL1EleRange",  40,0.,1.);
  TH1F *hL1HadRange   = new TH1F("hL1HadRange",  "hL1HadRange",  40,-2.,2.);

  //cuts histograms - unusable
  // TH1F *hL1EleEtaCut   = new TH1F("hL1EleEtaCut",  "hL1EleEtaCut",  20,-2.,2.);
  // TH1F *hL1ElePtCut   = new TH1F("hL1ElePtCut",  "hL1ElePtCut",  20,-2.,2.);

  //Histo wrt truth eta
  TH1F *hTrEleEtFab   = new TH1F("hTrEleEtFab",  "hTrEleEtFab", 40,5.,100.);
  TH1F *hTrEleEtaFab  = new TH1F("hTrEleEtaFab", "hTrEleEtaFab", 40,-5,5.);
  TH1F *hTrElePhiFab  = new TH1F("hTrElePhiFab", "hTrElePhiFab", 40,-5,5.);
  TH1F *hTrEleEtaCheck  = new TH1F("hTrEleEtaCheck", "hTrEleEtaCheck", 40,-5.,5.);
  
  TH1F *hTrmEleEtFab   = new TH1F("hTrmEleEtFab",  "hTrmEleEtFab", 40,5.,100.);
  TH1F *hTrmEleEtaFab = new TH1F("hTrmEleEtaFab", "hTrmEleEtaFab", 40,-5,5.);
  TH1F *hTrmElePhiFab = new TH1F("hTrmElePhiFab", "hTrmElePhiFab", 40,-5,5.); 
  TH1F *hTrmEleEtaCheck = new TH1F("hTrmEleEtaCheck", "hTrmEleEtaCheck", 40,-5.,5.);

  TH1F *hTrmEleEtFabIso = new TH1F("hTrmEleEtFabIso",  "hTrmEleEtFabIso", 40,5.,100.);
  TH1F *hTrmEleEtaFabIso = new TH1F("hTrmEleEtaFabIso", "hTrmEleEtaFabIso", 40,-5,5.);
  TH1F *hTrmElePhiFabIso = new TH1F("hTrmEleEtaIso", "hTrmEleEtaIso", 40,-5,5.);
  TH1F *hTrmEleEtaCheckIso = new TH1F("hTrmEleEtaCheckIso", "hTrmEleEtaCheckIso", 40,-5.,5.);

  TH1F *hL1EtaCut1   = new TH1F("hL1EtaCut1",  "hL1EtaCut1",  10, -6., 6.);
  TH1F *hL1EtaCut2   = new TH1F("hL1EtaCut2",  "hL1EtaCut2",  10, -6., 6.);
  TH1F *hL1EtaCut3   = new TH1F("hL1EtaCut3",  "hL1EtaCut3",  10, -6., 6.);
  TH1F *hL1EtaCut4   = new TH1F("hL1EtaCut4",  "hL1EtaCut4",  10, -6., 6.);
  TH1F *hL1EtaCut5   = new TH1F("hL1EtaCut5",  "hL1EtaCut5",  10, -6., 6.);

  TH1F *hL1PtCut1   = new TH1F("hL1PtCut1",  "hL1PtCut1",  20, 0., 110.);
  TH1F *hL1PtCut2   = new TH1F("hL1PtCut2",  "hL1PtCut2",  20, 0., 110.);
  TH1F *hL1PtCut3   = new TH1F("hL1PtCut3",  "hL1PtCut3",  20, 0., 110.);
  TH1F *hL1PtCut4   = new TH1F("hL1PtCut4",  "hL1PtCut4",  20, 0., 110.);
  TH1F *hL1PtCut5   = new TH1F("hL1PtCut5",  "hL1PtCut5",  20, 0., 110.);
  TH1F* h[10] = {hL1EtaCut1, hL1EtaCut2, hL1EtaCut3, hL1EtaCut4, hL1EtaCut5,
		hL1PtCut1, hL1PtCut2, hL1PtCut3, hL1PtCut4, hL1PtCut5};

  //2D Histograms, named after x-axis
  TH2F *hL1Treta = new TH2F("hL1Treta",  "hL1Treta",  40, 2.3, 5., 40, 0., 2.5);
  hL1Treta->SetTitle("L1Calo/Truth et against abs #eta ; Electron abs truth #eta ; Electron L1Calo/Truth #et");
  TH2F *hL1Tret = new TH2F("hL1Tret",  "hL1Tret",  40, 5., 100., 40, 5, 100.);
  hL1Tret->SetTitle("Reconstructed against Truth #et ; Electron truth #et ; Electron recon #et");
  //primitive var declarations

  //nentries = no of events
  //divide pjets/nentries
  int pJets = 0;

  //rate counters
  int elecPtCounter, emFlagCounter, elecPtRate, elecPt2Rate = 0;

  //flag for matched elec
  bool hasMatched;
  bool hasMatchedIso;
  
  //declaring list of isocuts
  //int  effiso[5] ={0}; 
  std::valarray<int> effiso(0,5);

  //Min et/eta cut values
  float etMin = 15000.;
  float etMax = 20000.;
  float etaMin = 2.4;
  float etaMax = 5.0;
  
  //container for all L1Calo objects that satisfy et cut
  std::vector<float> DRvals;

   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     //loop over all events
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(jentry%10000 == 0){
	cout<<"Processing "<<jentry<<endl;

	// additional debug printout
	if(debug ==1){
	  for (uint i = 0;i<Run3_jEle_Et->size();i++){
	    if(Run3_jEle_Et->at(i)>5000.){
	      cout<<"L1 electron: "<<Run3_jEle_Et->at(i)<<", "<<Run3_jEle_eta->at(i)<<", "<<Run3_jEle_phi->at(i)
		  <<", "<<Run3_jEle_iso->at(i)<<", "<<Run3_jEle_Had->at(i)<<endl;
	    }
	  }
	  for (uint i = 0;i<TruthEles_pt->size();i++){
	    if(TruthEles_motherID->at(i)==23 && TruthEles_pt->at(i)>10000. && TruthEles_barcode->at(i)<200000){
	      cout<<"Truth electron: "<<TruthEles_pt->at(i)<<", "<<TruthEles_eta->at(i)<<", "<<TruthEles_phi->at(i)
		  <<", "<<TruthEles_barcode->at(i)<<endl;
	    }
	  }
	}
      }// info for selected events

      int nl1el = Run3_jEle_Et->size();
      int ntel  = TruthEles_pt->size();
      int nl1gel = Run3_L1Ele_Et->size();
      int nl1iel = Run3_jEle_iso->size();
      int nl1hel = Run3_jEle_Had->size();
      
      vector<float> Run3_L1Ele_Et_New;
      for(int i = 0; i < nl1gel; i++){
	Run3_L1Ele_Et_New.push_back(Run3_L1Ele_Et->at(i)*1000.);
	//not a pointer, so access is through '.',
      }

      vector<float> Run3_jEle_Had_New;
      for(int i = 0; i < nl1hel; i++){
	Run3_jEle_Had_New.push_back(Run3_jEle_Had->at(i)*1000.);
      }

      for(int i = 0; i < Run3_jEle_Et->size(); i++){
	elecPtCounter++;
      }
      if(elecPtCounter > 1) elecPtRate++;
      if(elecPtCounter > 2) elecPt2Rate++;
      if(trig_L1_EM22VHI) emFlagCounter++;

      // Loop over truth electrons
      for (int i = 0;i<ntel;i++){
      	DRvals.clear();
	hasMatched = 0;
        hasMatchedIso = 0;
	//float TDEt = 0;
	float TDEthold = 0;
        float TDEtratio = 0;
	float TDEta = 0;
	std::vector<float> L1CaloEt_vals;

	if(TruthEles_motherID->at(i)==23 && TruthEles_pt->at(i)>10000.  && TruthEles_barcode->at(i)<200000  && fabs(TruthEles_eta->at(i))>etaMin){ 
	  
	  // Loop over L1Calo electrons
	  for (int j = 0;j<nl1el;j++){
            if (Run3_jEle_Et->at(j)<15000.) continue;
	    // Looking for matching L1Calo electron
	    float Dphi = fabs(TruthEles_phi->at(i)-Run3_jEle_phi->at(j));
	    if (Dphi > PI) Dphi = 2* PI - Dphi;   
	    float Deta = fabs(TruthEles_eta->at(i)-Run3_jEle_eta->at(j)); 
	    float DR = sqrt(pow(Dphi,2)+pow(Deta,2));
	    DRvals.push_back(DR);
	    L1CaloEt_vals.push_back(Run3_jEle_Et->at(j)/1000.);
	 
	    //cut on eta
	     if(fabs(TruthEles_eta->at(i))>etaMin && fabs(TruthEles_eta->at(i))<etaMax) TDEtratio = (Run3_jEle_Et->at(j)/1000)/(TruthEles_pt->at(i)/1000);  
	    //cut on et
	     if(TruthEles_pt->at(i)>etMax) TDEta = fabs(TruthEles_eta->at(i));
	    if(DR < 0.3){
	      //current - not satisfying DR < 0.3 singular condition, & eta doesn't work, cuts applied 
	       //  hL1Tret->Fill(TruthEles_pt->at(i)/1000, Run3_jEle_Et->at(j)/1000);
	      if(TDEta != 0 && TDEtratio != 0) hL1Treta->Fill(TDEta, TDEtratio);
	    }
	    //2d histos
	    //x-axis truth et, y-axis, l1calo et
	    //
	    //x-axis fabs eta, y-axis, ratio of L1Calo/truthet
	    //hL1Treta->Fill(
	  }

	  
	  //int DRmin_el = std::min_element(DRvals.begin(), DRvals.end() - DRvals.begin());
	  //float DRmin_el = *std::min_element(DRvals.begin(), DRvals.end() - DRvals.begin());

	  //ensures only 1 electron fills hL1 histos
	  float DRmin = 1.;
	  int DRmin_el = 1;
	  for (uint i = 0; i < DRvals.size(); i++){
	    if (DRvals[i] < DRmin){
	      DRmin = DRvals[i];
	      DRmin_el = i;
	    }
	  }
	  if(DRmin < 0.3){
	    hasMatched = 1;
	    if((Run3_jEle_iso->at(i)/Run3_jEle_Et->at(i))< 0.5){
	      hasMatchedIso = 1;
	    }
	    hL1EleEt->Fill(Run3_jEle_Et->at(DRmin_el)/1000.);
            hL1EleEta->Fill(Run3_jEle_eta->at(DRmin_el));
	    hL1ElePhi->Fill(Run3_jEle_phi->at(DRmin_el));
	  }
	  //Cut of 2.4 < |eta| < 5.0 for et histogram
	  if(fabs(TruthEles_eta->at(i))>etaMin && fabs(TruthEles_eta->at(i))<etaMax){
	    /* if(i < L1CaloEt_vals.size()){
		TDEtratio = L1CaloEt_vals[i]/TruthEles_pt->at(i)/1000.;
		}*/
	    hTrEleEt->Fill(TruthEles_pt->at(i)/1000.);
	    hTrEleEtFab->Fill(fabs(TruthEles_pt->at(i)/1000.));	
	    if(hasMatched) hTrmEleEt->Fill(TruthEles_pt->at(i)/1000.);
	    if(hasMatched) hTrmEleEtFab->Fill(fabs(TruthEles_pt->at(i)/1000.));
	    if(hasMatchedIso) hTrmEleEtFabIso->Fill(fabs(TruthEles_pt->at(i)/1000.));
	      //parse thru vector, divide each elements by i, add to histo
	   
	  }
	  //Cut of 20GeV on pt
	  if(TruthEles_pt->at(i)>etMax){
	    hTrEleEta->Fill(TruthEles_eta->at(i));
	    hTrEleEtaFab->Fill(fabs(TruthEles_eta->at(i)));
	    hTrEleEtaCheck->Fill(TruthEles_eta->at(i));
	    if(hasMatched) hTrmEleEta->Fill(TruthEles_eta->at(i));
	    if(hasMatched) hTrmEleEtaFab->Fill(fabs(TruthEles_eta->at(i)));
	    if(hasMatched) hTrmEleEtaCheck->Fill(TruthEles_eta->at(i));	   
	    if(hasMatchedIso) hTrmEleEtaFabIso->Fill(fabs(TruthEles_eta->at(i)));
	    if(hasMatchedIso) hTrmEleEtaCheckIso->Fill(TruthEles_eta->at(i));
	    // if(hasMatched) hL1Treta->Fill(fabs(TruthEles_eta->at(i)), TDEtratio);
	    if(hasMatched) hL1Tret->Fill(TruthEles_pt->at(i)/1000, L1CaloEt_vals[i]);
	  }
	  //Cut of 2.4 < |eta| < 5.0 for et histogram & Cut of 20GeV on pt
	  if(TruthEles_pt->at(i)>etMax && fabs(TruthEles_eta->at(i))>etaMin && fabs(TruthEles_eta->at(i))<etaMax){
	    hTrElePhi->Fill(TruthEles_phi->at(i));
	    hTrElePhiFab->Fill(fabs(TruthEles_phi->at(i)));	    	   
	    if(hasMatched) hTrmElePhi->Fill(TruthEles_phi->at(i));
	    if(hasMatched) hTrmElePhiFab->Fill(fabs(TruthEles_phi->at(i)));
	    if(hasMatchedIso) hTrmElePhiFabIso->Fill(fabs(TruthEles_phi->at(i)));
	  }


	      // Loop over L1Calo eFEX  electrons
	  for (int j = 0;j<nl1gel;j++){
            if (Run3_L1Ele_Et_New.at(j)<10000.) continue;
	    // Looking for matching L1Calo electron
	    float Dphi = fabs(TruthEles_phi->at(i)-Run3_L1Ele_phi->at(j));
	    if (Dphi > PI) Dphi = 2* PI - Dphi;   
	    float Deta = fabs(TruthEles_eta->at(i)-Run3_L1Ele_eta->at(j)); 
	    float DR = sqrt(pow(Dphi,2)+pow(Deta,2));
	    //// Matching criterion, can be refined further
            if (DR < 0.3){
	      // Matched L1Calo eFEX  electrons
              hL1EfexEleEt->Fill(Run3_L1Ele_Et_New.at(j)/1000.);
              hL1EfexEleEta->Fill(Run3_L1Ele_eta->at(j));
	      hL1EfexElePhi->Fill(Run3_L1Ele_phi->at(j));
	      // Matched Truth  electrons
	      hTrmEfexEleEt->Fill(TruthEles_pt->at(i)/1000.);
	      hTrmEfexEleEta->Fill(TruthEles_eta->at(i));
	      hTrmEfexElePhi->Fill(TruthEles_phi->at(i));
	    }
	  }   
	   // end loop over L1Calo electrons
	}
      }// end loop over truth electrons
      
      float isorange = 0;
      float hadrange = 0;

      // float isocuts[10] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
      float isocuts[5] = { 0.3, 0.4, 0.5, 0.6,  0.7};
      float cutrejrate[5] ={0};

      //rejection vars, loop over iso var
       for (int j = 0;j<nl1iel;j++){
	 if (Run3_jEle_Et->at(j)>20000. && fabs(Run3_jEle_eta->at(j)) > 2.4 && fabs(Run3_jEle_eta->at(j)) < 5.0) {
               hL1EleIso->Fill(Run3_jEle_iso->at(j));
	       hL1EleHad->Fill(Run3_jEle_Had_New.at(j));
	       isorange = Run3_jEle_iso->at(j)/Run3_jEle_Et->at(j);
	       // cout<<isorange<<endl;
	       // range of iso var of ele iso / ele et
	       hL1EleRange->Fill(isorange);
	       hadrange = Run3_jEle_Had_New.at(j)/Run3_jEle_Had_New.at(j);
	       hL1HadRange->Fill(hadrange);
	   
	       if (0 < isorange && isorange < isocuts[0] ){
		  effiso[0] += 1;
	  	  h[4]->Fill(Run3_jEle_eta->at(j));
		  h[9]->Fill(Run3_jEle_Et->at(j)/1000.);
	       }
	       if (0 < isorange && isorange < isocuts[1]){
		  effiso[1] += 1;
		  h[3]->Fill(Run3_jEle_eta->at(j));
		  h[8]->Fill(Run3_jEle_Et->at(j)/1000.);

	       }
	       if (0 < isorange && isorange < isocuts[2]){
		  effiso[2]++;
		  h[2]->Fill(Run3_jEle_eta->at(j));
		  h[7]->Fill(Run3_jEle_Et->at(j)/1000.);
	       }
	       if (0 < isorange && isorange < isocuts[3]){
		  effiso[3]++;
		  h[1]->Fill(Run3_jEle_eta->at(j));
		  h[6]->Fill(Run3_jEle_Et->at(j)/1000.);
	       }
	       if (0 < isorange && isorange < isocuts[4]){
		  effiso[4]++;
		  h[0]->Fill(Run3_jEle_eta->at(j));
		  h[5]->Fill(Run3_jEle_Et->at(j)/1000.);

	       }
	       /*   
		else if (isorange > isocuts[4] && isorange < isocuts[5]){
		  effiso6++;
	       }
	        else if (isorange > isocuts[5] && isorange < isocuts[6]){
		  effiso7++;
	       }
	        else if (isorange > isocuts[6] && isorange < isocuts[7]){
		  effiso8++;
	       }
	        else if (isorange > isocuts[7] && isorange < isocuts[8]){
		  effiso9++;
	       }
	        else if (isorange > isocuts[9] && isorange < isocuts[10]){
		  effiso10++;
		  } */
	 }
       }
       for (int j = 0;j<nl1el;j++){
	 if (Run3_jEle_Et->at(j)<10000.){
	    pJets++;
	      // Matched had jets electrons
             //  hL1EleEtHad->Fill(Run3_jEle_Et->at(j)/1000.);
             //  hL1EleEtaHad->Fill(Run3_jEle_eta->at(j));
	      // hL1ElePhiHad->Fill(Run3_jEle_phi->at(j));
	 }
      }   
      
   }// end event loop

   //Counters/Rates
   cout << elecPtCounter << endl;
   cout << emFlagCounter << endl;
   cout << elecPtRate << endl;
   cout << elecPt2Rate << endl;

   cout << "cuts rejections (increments of .2"<<endl;
   cout << effiso[0] << endl;
   for(int i=0; i<effiso.size(); i++){
     int j = effiso[i];
     cout << j/nentries << endl;
     cout << j << endl;
   }

   float rJets = float(pJets)/float(nentries);
   //output 

   cout<<nentries<<endl;
   cout<<pJets<<endl;
   cout<<rJets<<endl;
   
   //Write histograms
   hL1EleEt->Write();
   hL1EleEta->Write();
   hL1ElePhi->Write();
   hL1EfexEleEt->Write();
   hL1EfexEleEta->Write();
   hL1EfexElePhi->Write();
   hTrEleEt->Write();
   hTrEleEta->Write();
   hTrElePhi->Write();
   hTrmEleEta->Write();
   hTrmEleEt->Write();
   hTrmElePhi->Write();
   hTrmEfexEleEt->Write();
   hTrmEfexEleEta->Write();
   hTrmEfexElePhi->Write();
  //  hL1EleEtHad->Write();
   // hL1EleEtaHad->Write();
  //  hL1ElePhiHad->Write();
   hL1EleIso->Write();
   hL1EleHad->Write();
   hL1EleRange->Write();
   hL1HadRange->Write();
   hTrEleEtFab->Write();  
   hTrEleEtaFab->Write();  
   hTrElePhiFab->Write();  
   hTrEleEtaCheck->Write();
   hTrmEleEtFab->Write();
   hTrmEleEtaFab->Write();
   hTrmElePhiFab->Write();
   hTrmEleEtaCheck->Write();
   hTrmEleEtFabIso->Write();
   hTrmEleEtaFabIso->Write();
   hTrmElePhiFabIso->Write();
   hTrmEleEtaCheckIso->Write();
   // hL1EtaCut->Write();
   // hL1PtCut->Write();
   for(int i=0; i<10; i++){
     h[i]->Write();
   }  
   //2d plots
   hL1Tret->Write();
   hL1Treta->Write(); 
   // Efficiency plots

   TGraphAsymmErrors tgEleta; 
   tgEleta.Divide(hTrmEleEta,hTrEleEta, "cl=0.683 b(1,1) mode"); 

   can1->SetLeftMargin(0.16);
   can1->SetBottomMargin(0.16);

   tgEleta.GetXaxis()->SetTitle("Electron #eta");
   tgEleta.GetYaxis()->SetTitle("Efficiency");
   tgEleta.SetMarkerStyle(20);
   tgEleta.SetMarkerColor(46);
   tgEleta.GetXaxis()->SetTitleOffset(1.2);
   tgEleta.GetYaxis()->SetTitleOffset(1.2);

   can1->cd(); 
   tgEleta.Draw("AP"); 

   can1->Write(); 
   can1->Print("EffElEta.png","png");  

   //create plot of effeciency to find truth electrion in L1Calo trigger as a func of phi variable
  
   TGraphAsymmErrors tgElet; 
   tgElet.Divide(hTrmEleEt, hTrEleEt, "cl=0.683 b(1,1) mode"); 

   can1->SetLeftMargin(0.16);
   can1->SetBottomMargin(0.16);

   tgElet.GetXaxis()->SetTitle("Electron et");
   tgElet.GetYaxis()->SetTitle("Efficiency");
   tgElet.SetMarkerStyle(20);
   tgElet.SetMarkerColor(46);
   tgElet.GetXaxis()->SetTitleOffset(1.2);
   tgElet.GetYaxis()->SetTitleOffset(1.2);

   can1->cd(); 
   tgElet.Draw("AP"); 

   can1->Write(); 
   can1->Print("EffElet.png","png"); 

 
   TGraphAsymmErrors tgElphi; 
   tgElphi.Divide(hTrmElePhi, hTrElePhi, "cl=0.683 b(1,1) mode"); 

   can1->SetLeftMargin(0.16);
   can1->SetBottomMargin(0.16);

   tgElphi.GetXaxis()->SetTitle("Electron #phi");
   tgElphi.GetYaxis()->SetTitle("Efficiency");
   tgElphi.SetMarkerStyle(20);
   tgElphi.SetMarkerColor(46);
   tgElphi.GetXaxis()->SetTitleOffset(1.2);
   tgElphi.GetYaxis()->SetTitleOffset(1.2);

   can1->cd(); 
   tgElphi.Draw("AP"); 

   can1->Write(); 
   can1->Print("EffElphi.png","png"); 

   TGraphAsymmErrors tgEfexEleta; 
   tgEfexEleta.Divide(hTrmEfexEleEta, hTrEleEta, "cl=0.683 b(1,1) mode"); 

   can1->SetLeftMargin(0.16);
   can1->SetBottomMargin(0.16);

   tgEfexEleta.GetXaxis()->SetTitle("Electron #eta");
   tgEfexEleta.GetYaxis()->SetTitle("Efficiency");
   tgEfexEleta.SetMarkerStyle(20);
   tgEfexEleta.SetMarkerColor(46);
   tgEfexEleta.GetXaxis()->SetTitleOffset(1.2);
   tgEfexEleta.GetYaxis()->SetTitleOffset(1.2);

   can1->cd(); 
   tgEfexEleta.Draw("AP"); 

   can1->Write(); 
   can1->Print("EffEfexEleta.png","png"); 

   TGraphAsymmErrors tgEfexElet; 
   tgEfexElet.Divide(hTrmEfexEleEt,hTrEleEt, "cl=0.683 b(1,1) mode"); 

   can1->SetLeftMargin(0.16);
   can1->SetBottomMargin(0.16);

   tgEfexElet.GetXaxis()->SetTitle("Electron et");
   tgEfexElet.GetYaxis()->SetTitle("Efficiency");
   tgEfexElet.SetMarkerStyle(20);
   tgEfexElet.SetMarkerColor(46);
   tgEfexElet.GetXaxis()->SetTitleOffset(1.2);
   tgEfexElet.GetYaxis()->SetTitleOffset(1.2);

   can1->cd(); 
   tgEfexElet.Draw("AP"); 

   can1->Write(); 
   can1->Print("EffEfexElet.png","png"); 

 
   TGraphAsymmErrors tgEfexElphi; 
   tgEfexElphi.Divide(hTrmEfexElePhi, hTrElePhi, "cl=0.683 b(1,1) mode"); 

   can1->SetLeftMargin(0.16);
   can1->SetBottomMargin(0.16);

   tgEfexElphi.GetXaxis()->SetTitle("Electron #phi");
   tgEfexElphi.GetYaxis()->SetTitle("Efficiency");
   tgEfexElphi.SetMarkerStyle(20);
   tgEfexElphi.SetMarkerColor(46);
   tgEfexElphi.GetXaxis()->SetTitleOffset(1.2);
   tgEfexElphi.GetYaxis()->SetTitleOffset(1.2);

   can1->cd(); 
   tgEfexElphi.Draw("AP"); 

   can1->Write(); 
   can1->Print("EffEfexElphi.png","png"); 

   //Efficiencies wrt truth |eta| instead of eta
   TGraphAsymmErrors tgEleEtaFab20; 
   tgEleEtaFab20.Divide(hTrmEleEtaFab, hTrEleEtaFab, "cl=0.683 b(1,1) mode"); 

   can1->SetLeftMargin(0.16);
   can1->SetBottomMargin(0.16);

   tgEleEtaFab20.GetXaxis()->SetTitle("Electron truth #eta");
   tgEleEtaFab20.GetYaxis()->SetTitle("Efficiency");
   tgEleEtaFab20.SetMarkerStyle(20);
   tgEleEtaFab20.SetMarkerColor(46);
   tgEleEtaFab20.GetXaxis()->SetTitleOffset(1.2);
   tgEleEtaFab20.GetYaxis()->SetTitleOffset(1.2);

   can1->cd(); 
   tgEleEtaFab20.Draw("AP"); 

   can1->Write(); 
   can1->Print("tgEleEtaFab20.png","png"); 

   TGraphAsymmErrors tgEleEtFab; 
   tgEleEtFab.Divide(hTrmEleEtFab, hTrEleEtFab, "cl=0.683 b(1,1) mode"); 

   can1->SetLeftMargin(0.16);
   can1->SetBottomMargin(0.16);

   tgEleEtFab.GetXaxis()->SetTitle("Electron truth #et");
   tgEleEtFab.GetYaxis()->SetTitle("Efficiency");
   tgEleEtFab.SetMarkerStyle(20);
   tgEleEtFab.SetMarkerColor(46);
   tgEleEtFab.GetXaxis()->SetTitleOffset(1.2);
   tgEleEtFab.GetYaxis()->SetTitleOffset(1.2);

   can1->cd(); 
   tgEleEtFab.Draw("AP"); 

   can1->Write(); 
   can1->Print("tgEleEtFab.png","png"); 

   TGraphAsymmErrors tgElePhiFab; 
   tgElePhiFab.Divide(hTrmElePhiFab, hTrElePhiFab, "cl=0.683 b(1,1) mode"); 

   can1->SetLeftMargin(0.16);
   can1->SetBottomMargin(0.16);

   tgElePhiFab.GetXaxis()->SetTitle("Electron truth #phi");
   tgElePhiFab.GetYaxis()->SetTitle("Efficiency");
   tgElePhiFab.SetMarkerStyle(20);
   tgElePhiFab.SetMarkerColor(46);
   tgElePhiFab.GetXaxis()->SetTitleOffset(1.2);
   tgElePhiFab.GetYaxis()->SetTitleOffset(1.2);

   can1->cd(); 
   tgElePhiFab.Draw("AP"); 

   can1->Write(); 
   can1->Print("tgElePhiFab.png","png"); 

    TGraphAsymmErrors tgEleEtaCheck; 
   tgEleEtaCheck.Divide(hTrmEleEtaCheck, hTrEleEtaCheck, "cl=0.683 b(1,1) mode"); 

   can1->SetLeftMargin(0.16);
   can1->SetBottomMargin(0.16);

   tgEleEtaCheck.GetXaxis()->SetTitle("Electron truth #eta");
   tgEleEtaCheck.GetYaxis()->SetTitle("Efficiency");
   tgEleEtaCheck.SetMarkerStyle(20);
   tgEleEtaCheck.SetMarkerColor(46);
   tgEleEtaCheck.GetXaxis()->SetTitleOffset(1.2);
   tgEleEtaCheck.GetYaxis()->SetTitleOffset(1.2);

   can1->cd(); 
   tgEleEtaFab20.Draw("AP"); 

   can1->Write(); 
   can1->Print("tgEleEtaCheck.png","png"); 
   
   //efficiencies wrt truth |eta| instead of eta, + isolation cut 
   TGraphAsymmErrors tgEleEtaFabIso20; 
   tgEleEtaFabIso20.Divide(hTrmEleEtaFabIso, hTrEleEtaFab, "cl=0.683 b(1,1) mode"); 

   can1->SetLeftMargin(0.16);
   can1->SetBottomMargin(0.16);

   tgEleEtaFabIso20.GetXaxis()->SetTitle("Electron truth #eta");
   tgEleEtaFabIso20.GetYaxis()->SetTitle("Efficiency");
   tgEleEtaFabIso20.SetMarkerStyle(20);
   tgEleEtaFabIso20.SetMarkerColor(46);
   tgEleEtaFabIso20.GetXaxis()->SetTitleOffset(1.2);
   tgEleEtaFabIso20.GetYaxis()->SetTitleOffset(1.2);

   can1->cd(); 
   tgEleEtaFabIso20.Draw("AP"); 

   can1->Write(); 
   can1->Print("tgEleEtaFabIso20.png","png"); 

   TGraphAsymmErrors tgEleEtFabIso; 
   tgEleEtFabIso.Divide(hTrmEleEtFabIso, hTrEleEtFab, "cl=0.683 b(1,1) mode"); 

   can1->SetLeftMargin(0.16);
   can1->SetBottomMargin(0.16);

   tgEleEtFabIso.GetXaxis()->SetTitle("Electron truth #et");
   tgEleEtFabIso.GetYaxis()->SetTitle("Efficiency");
   tgEleEtFabIso.SetMarkerStyle(20);
   tgEleEtFabIso.SetMarkerColor(46);
   tgEleEtFabIso.GetXaxis()->SetTitleOffset(1.2);
   tgEleEtFabIso.GetYaxis()->SetTitleOffset(1.2);

   can1->cd(); 
   tgEleEtFabIso.Draw("AP"); 

   can1->Write(); 
   can1->Print("tgEleEtFabIso.png","png"); 

   TGraphAsymmErrors tgElePhiFabIso; 
   tgElePhiFabIso.Divide(hTrmElePhiFabIso, hTrElePhiFab, "cl=0.683 b(1,1) mode"); 

   can1->SetLeftMargin(0.16);
   can1->SetBottomMargin(0.16);

   tgElePhiFabIso.GetXaxis()->SetTitle("Electron truth #phi");
   tgElePhiFabIso.GetYaxis()->SetTitle("Efficiency");
   tgElePhiFabIso.SetMarkerStyle(20);
   tgElePhiFabIso.SetMarkerColor(46);
   tgElePhiFabIso.GetXaxis()->SetTitleOffset(1.2);
   tgElePhiFabIso.GetYaxis()->SetTitleOffset(1.2);

   can1->cd(); 
   tgElePhiFab.Draw("AP"); 

   can1->Write(); 
   can1->Print("tgElePhiFabIso.png","png"); 

    TGraphAsymmErrors tgEleEtaCheckIso; 
   tgEleEtaCheck.Divide(hTrmEleEtaCheckIso, hTrEleEtaCheck, "cl=0.683 b(1,1) mode"); 

   can1->SetLeftMargin(0.16);
   can1->SetBottomMargin(0.16);

   tgEleEtaCheckIso.GetXaxis()->SetTitle("Electron truth #eta");
   tgEleEtaCheckIso.GetYaxis()->SetTitle("Efficiency");
   tgEleEtaCheckIso.SetMarkerStyle(20);
   tgEleEtaCheckIso.SetMarkerColor(46);
   tgEleEtaCheckIso.GetXaxis()->SetTitleOffset(1.2);
   tgEleEtaCheckIso.GetYaxis()->SetTitleOffset(1.2);

   can1->cd(); 
   tgEleEtaFab20.Draw("AP"); 

   can1->Write(); 
   can1->Print("tgEleEtaCheck.png","png"); 
   /*
   TGraphAsymmErrors tgEleEtaCheck; 
   tgEleEtaCheck.Divide(hTrmEleEta, hTrEleEtaCheck, "cl=0.683 b(1,1) mode"); 

   can1->SetLeftMargin(0.16);
   can1->SetBottomMargin(0.16);

   tgEleEtaCheck.GetXaxis()->SetTitle("Electron truth #et");
   tgEleEtaCheck.GetYaxis()->SetTitle("Efficiency");
   tgEleEtaCheck.SetMarkerStyle(20);
   tgEleEtaCheck.SetMarkerColor(46);
   tgEleEtaCheck.GetXaxis()->SetTitleOffset(1.2);
   tgEleEtaCheck.GetYaxis()->SetTitleOffset(1.2);

   can1->cd(); 
   tgEleEtaCheck.Draw("AP"); 

   can1->Write(); 
   can1->Print("tgEleEtaCheck.png","png"); */
  

}// end Loop method

