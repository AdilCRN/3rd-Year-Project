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
  TH1F *hL1EleEt   = new TH1F("hL1EleEt",  "hL1EleEt",  20,10.,100.);
  TH1F *hL1EleEta  = new TH1F("hL1EleEta", "hL1EleEta", 20,-5.,5.);
  TH1F *hL1ElePhi  = new TH1F("hL1ElePhi", "hL1ElePhi", 20,-5.,5.);
  //eFEX electrons
  TH1F *hL1EfexEleEt   = new TH1F("hL1EfexEleEt",  "hL1EfexEleEt",  20,10.,100.);
  TH1F *hL1EfexEleEta  = new TH1F("hL1EfexEleEta", "hL1EfexEleEta", 20,-5.,5.);
  TH1F *hL1EfexElePhi  = new TH1F("hL1EfexElePhi", "hL1EfexElePhi", 20,-5.,5.);
  // All Truth
  TH1F *hTrEleEt   = new TH1F("hTrEleEt",  "hTrEleEt",  20,10.,100.);
  TH1F *hTrEleEta  = new TH1F("hTrEleEta", "hTrEleEta", 20,-5.,5.);
  TH1F *hTrElePhi  = new TH1F("hTrElePhi", "hTrElePhi", 20,-5.,5.);
  // Truth matched to L1Calo

  TH1F *hTrmEleEt  = new TH1F("hTrmEleEt", "hTrmEleEt", 20,10.,100.);
  TH1F *hTrmEleEta  = new TH1F("hTrmEleEta", "hTrmEleEta", 20,-5.,5.);
  TH1F *hTrmElePhi  = new TH1F("hTrmElePhi", "hTrmElePhi", 20,-5.,5.);

  // Truth matched to L1Calo eFEX
  TH1F *hTrmEfexEleEt  = new TH1F("hTrmEfexEleEt", "hTrmEfexEleEt", 20,10.,100.);
  TH1F *hTrmEfexEleEta  = new TH1F("hTrmEfexEleEta", "hTrmEfexEleEta", 20,-5.,5.);
  TH1F *hTrmEfexElePhi  = new TH1F("hTrmEfexElePhi", "hTrmEfexElePhi", 20,-5.,5.);



  // Hadronic Jet
//   TH1F *jEleEt   = new TH1F("hL1EleEtHad",  "hL1EleEtHad",  20,10.,100.);
//   TH1F *jEleEta  = new TH1F("hL1EleEtaHad", "hL1EleEtaHad", 20,-5.,5.);
//   TH1F *jElePhi  = new TH1F("hL1ElePhiHad", "hL1ElePhiHad", 20,-5.,5.);

  //
//   TH1F *hL1EleEt   = new TH1F("hjEleEt",  "hjEleEt",  20,10.,100.);
//   TH1F *hjEleEta  = new TH1F("hjEleEta", "hjEleEta", 20,-5.,5.);
//   TH1F *hjElePhi  = new TH1F("hjElePhi", "hjElePhi", 20,-5.,5.);
  //Jet rejection variables
  TH1F *hL1EleIso   = new TH1F("hL1EleIso",  "hL1EleIso",  20,-10000.,40000.);
  TH1F *hL1EleHad   = new TH1F("hL1EleHad",  "hL1EleHad",  20,-1.,1.);

  TH1F *hL1EleRange   = new TH1F("hL1EleRange",  "hL1EleRange",  20,0.,1.);
  TH1F *hL1HadRange   = new TH1F("hL1HadRange",  "hL1HadRange",  20,-2.,2.);

  //cuts histograms - unusable
  // TH1F *hL1EleEtaCut   = new TH1F("hL1EleEtaCut",  "hL1EleEtaCut",  20,-2.,2.);
  // TH1F *hL1ElePtCut   = new TH1F("hL1ElePtCut",  "hL1ElePtCut",  20,-2.,2.);

  //Creating histograms for varying pt and eta cuts
  /* TString cutNames[9];
   cutNames[0] = "EtaCut1", cutNames[1] = "EtaCut2", cutNames[2] = "EtaCut3", cutNames[3] = "EtaCut4", cutNames[4] = "EtaCut5";
     cutNames[5] = "PtCut1", cutNames[6] = "PtCut2", cutNames[7] = "PtCut3", cutNames[8] = "PtCut4", cutNames[9] = "PtCut5"; */
  
  //int effiso[5] ={0}; 
  std::valarray<int> effiso(0,5);

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
   
  /*
   int m=10;
   int j = 1;
   TH1F** h=new TH1F*[m];

   for (int i=0; i<m;i++){
      if(i<5){
       cutNames[i].Append("hL1EtaCut", j); 
       cout << cutNames[i];
     }
     else{
       j = 1;
       cutNames[i].Append("hL1PtCut", j);
       }
     h[i] = new TH1F("a", "a", 20., -5., 5.);
     j++;
     }*/

  int pJets = 0;

  //nentries = no of events
  //divide pjets/nentries

   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
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

           // Loop over truth electrons
      for (int i = 0;i<ntel;i++){
	
	if(TruthEles_motherID->at(i)==23 && TruthEles_pt->at(i)>10000.  && TruthEles_barcode->at(i)<200000  && fabs(TruthEles_eta->at(i))>2.0){

	  // Fill truth histograms 
	  hTrEleEt->Fill(TruthEles_pt->at(i)/1000.);
	  hTrEleEta->Fill(TruthEles_eta->at(i));
	  hTrElePhi->Fill(TruthEles_phi->at(i));

	  // Loop over L1Calo electrons
	  for (int j = 0;j<nl1el;j++){
            if (Run3_jEle_Et->at(j)<10000.) continue;
	    // Looking for matching L1Calo electron
	    float Dphi = fabs(TruthEles_phi->at(i)-Run3_jEle_phi->at(j));
	    if (Dphi > PI) Dphi = 2* PI - Dphi;   
	    float Deta = fabs(TruthEles_eta->at(i)-Run3_jEle_eta->at(j)); 
	    float DR = sqrt(pow(Dphi,2)+pow(Deta,2));
	    //// Matching criterion, can be refined further
            if (DR < 0.3){
	      // Matched L1Calo electrons
              hL1EleEt->Fill(Run3_jEle_Et->at(j)/1000.);
              hL1EleEta->Fill(Run3_jEle_eta->at(j));
	      hL1ElePhi->Fill(Run3_jEle_phi->at(j));
	      // Matched Truth  electrons
	      hTrmEleEt->Fill(TruthEles_pt->at(i)/1000.);
	      hTrmEleEta->Fill(TruthEles_eta->at(i));
	      hTrmElePhi->Fill(TruthEles_phi->at(i));
	    }
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

      //rejection vars
       for (int j = 0;j<nl1iel;j++){
	 if (Run3_jEle_Et->at(j)>20000. && fabs(Run3_jEle_eta->at(j)) > 2.5 && fabs(Run3_jEle_eta->at(j)) < 4.5) {
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

   cout << "cuts rejections (increments of .2"<<endl;
   cout << effiso[0] << endl;
   cout << "why wtf where" << endl;
   for(int i=0; i<effiso.size(); i++){
     int j = effiso[i];
     cout << j/nentries << endl;
     cout << j << endl;
   }
   cout << "pls" << endl;

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

   // hL1EtaCut->Write();
   // hL1PtCut->Write();

   for(int i=0; i<10; i++){
     h[i]->Write();
   }  
  
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

}// end Loop method

