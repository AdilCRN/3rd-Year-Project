#include <TF1.h>
#include <iostream>
#include <math.h>
#include <string>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>

using namespace std;

void DrawPlots(){


  int option = 0;// not used yet

  // declare Canvas instance
  TCanvas *can1 = new TCanvas("Canvas1","Histograms",600,600);
  TCanvas *can2 = new TCanvas("Canvas2","Histograms",600,600);


  // select the histogram inputfile
  TFile* hisfile   = new TFile("./histos.root", "OPEN");
  TFile* hisfile2   = new TFile("./histos2.root", "OPEN");

  //======================
  // User input
  //======================
  // select the truth and L1Calo histogram names
  string sTru("hTrEleEta");
  string sL1("hL1EleEta");
  string sTrm("hTrmEleEta");
  string xtitle("Electron #eta");
  string plottitle("EleEta");

  string sTruet("hTrEleEt");
  string sL1et("hL1EleEt");
  string sTrmet("hTrmEleEt");
  string xtitleet("Electron #et");
  string plottitleet("EleEt");

  string sTruphi("hTrElePhi");
  string sL1phi("hL1ElePhi");
  string sTrmphi("hTrmElePhi");
  string xtitlephi("Electron #phi");
  string plottitlephi("ElePhi");

  string hL1EC1("hL1EtaCut1");
  string xtitlecuteta("Electron #etacut (rad)");
  string plottitlecuteta("EtaCutRates");
  string hL1EC2("hL1EtaCut2");
  string hL1EC3("hL1EtaCut3");
  string hL1EC4("hL1EtaCut4");
  string hL1EC5("hL1EtaCut5");

  string hL1PC1("hL1PtCut1");
  string xtitlecutpt("Electron #ptcut (MeV)");
  string plottitlecutpt("PtCutRates");
  string hL1PC2("hL1PtCut2");
  string hL1PC3("hL1PtCut3");
  string hL1PC4("hL1PtCut4");
  string hL1PC5("hL1PtCut5");


  //======================

  // get the histograms
  TH1F* hTru =  (TH1F*)hisfile->Get(sTru.data());
  TH1F* hL1  =  (TH1F*)hisfile->Get(sL1.data());
  TH1F* hTrm =  (TH1F*)hisfile->Get(sTrm.data());

  TH1F* hTruet =  (TH1F*)hisfile->Get(sTruet.data());
  TH1F* hL1et  =  (TH1F*)hisfile->Get(sL1et.data());
  TH1F* hTrmet =  (TH1F*)hisfile->Get(sTrmet.data());

  TH1F* hTruphi =  (TH1F*)hisfile->Get(sTruphi.data());
  TH1F* hL1phi  =  (TH1F*)hisfile->Get(sL1phi.data());
  TH1F* hTrmphi =  (TH1F*)hisfile->Get(sTrmphi.data());

  //Cut histos for All Electrons, Cut1 - tightest,  Cut5 - loosest (cut)
  TH1F* hL1EtaCut1  =  (TH1F*)hisfile->Get(hL1EC1.data());
  TH1F* hL1EtaCut2  =  (TH1F*)hisfile->Get(hL1EC2.data());
  TH1F* hL1EtaCut3  =  (TH1F*)hisfile->Get(hL1EC3.data());
  TH1F* hL1EtaCut4  =  (TH1F*)hisfile->Get(hL1EC4.data());
  TH1F* hL1EtaCut5  =  (TH1F*)hisfile->Get(hL1EC5.data());

  TH1F* hL1PtCut1   =  (TH1F*)hisfile->Get(hL1PC1.data());
  TH1F* hL1PtCut2   =  (TH1F*)hisfile->Get(hL1PC2.data());
  TH1F* hL1PtCut3   =  (TH1F*)hisfile->Get(hL1PC3.data());
  TH1F* hL1PtCut4   =  (TH1F*)hisfile->Get(hL1PC4.data());
  TH1F* hL1PtCut5   =  (TH1F*)hisfile->Get(hL1PC5.data());

  //Stack of hists for Eta and Pt
  // THStack *hsEta = new THStack("hsEta","");
  // THStack *hsPt = new THStack("hsPt","");

  // set the line style
  hTru->SetLineWidth(2);
  hTru->SetLineColor(1);
  hTru->SetLineWidth(2);
  hL1->SetLineWidth(2);
  hL1->SetLineColor(38);
  hL1->SetLineWidth(2);
  hTrm->SetLineWidth(2);
  hTrm->SetLineColor(46);
  hTrm->SetLineWidth(2);

  hL1EtaCut1->SetLineColor(1);
  hL1EtaCut2->SetLineColor(38);
  hL1EtaCut3->SetLineColor(46);
  hL1EtaCut4->SetLineColor(29);
  hL1EtaCut5->SetLineColor(41);

  hL1PtCut1->SetLineColor(1);
  hL1PtCut2->SetLineColor(38);
  hL1PtCut3->SetLineColor(46);
  hL1PtCut4->SetLineColor(29);
  hL1PtCut5->SetLineColor(41);

  hTruet->SetLineWidth(2);
  hTruet->SetLineColor(1);
  hTruet->SetLineWidth(2);
  hL1et->SetLineWidth(2);
  hL1et->SetLineColor(38);
  hL1et->SetLineWidth(2);
  hTrmet->SetLineWidth(2);
  hTrmet->SetLineColor(46);
  hTrmet->SetLineWidth(2);

  hTruphi->SetLineWidth(2);
  hTruphi->SetLineColor(1);
  hTruphi->SetLineWidth(2);
  hL1phi->SetLineWidth(2);
  hL1phi->SetLineColor(38);
  hL1phi->SetLineWidth(2);
  hTrmphi->SetLineWidth(2);
  hTrmphi->SetLineColor(46);
  hTrmphi->SetLineWidth(2);
  
  /*add Eta/Pt histos to relevant Stacks
  hsEta->Add(hL1EtaCut1);
  hsEta->Add(hL1EtaCut2);
  hsEta->Add(hL1EtaCut3);
  hsEta->Add(hL1EtaCut4);
  hsEta->Add(hL1EtaCut5);

  hsPt->Add(hL1PtCut1);
  hsPt->Add(hL1PtCut2);
  hsPt->Add(hL1PtCut3);
  hsPt->Add(hL1PtCut4);
  hsPt->Add(hL1PtCut5);
  */

  // Draw the histograms

  can1->SetLeftMargin(0.16);
  can1->SetBottomMargin(0.16);
  can1->cd();

  hTru->Draw("");
  hTrm->Draw("hist same");
  hL1->Draw("hist same");

  string plot1 = plottitle + string(".png");
  can1->Print(plot1.data(),"png");

  //
  hTruet->Draw("");
  hTrmet->Draw("hist same");
  hL1et->Draw("hist same");

  string plot2 = plottitleet + string(".png");
  can1->Print(plot2.data(),"png");

  // 
  hTruphi->Draw("");
  hTrmphi->Draw("hist same");
  hL1phi->Draw("hist same");

  string plot3 = plottitlephi + string(".png");
  can1->Print(plot3.data(),"png");

  hL1EtaCut1->GetXaxis()->SetTitle(xtitlecuteta.data());
  hL1EtaCut1->Draw("");
  hL1EtaCut2->Draw("hist same");
  hL1EtaCut3->Draw("hist same");
  hL1EtaCut4->Draw("hist same");
  hL1EtaCut5->Draw("hist same");

  string plot7 = plottitlecuteta + string(".png");
  can1->Print(plot7.data(),"png");

  //add legend to etacut combined histo
  auto legend1 = new TLegend(10, 10, 20, 20);
   legend1->SetHeader("Cut Values","R"); // option "C" allows to center the header
   legend1->AddEntry(hL1PtCut1,"Cut: 0.0-0.7","l");
   legend1->AddEntry(hL1PtCut2,"Cut: 0.0-0.6","l");
   legend1->AddEntry(hL1PtCut3,"Cut: 0.0-0.5","l");
   legend1->AddEntry(hL1PtCut4,"Cut: 0.0-0.4","l");
   legend1->AddEntry(hL1PtCut5,"Cut: 0.0-0.3","l");

   legend1->Draw();

  hL1PtCut1->GetXaxis()->SetTitle(xtitlecutpt.data());
  hL1PtCut1->Draw("");
  hL1PtCut2->Draw("hist same");
  hL1PtCut3->Draw("hist same");
  hL1PtCut4->Draw("hist same");
  hL1PtCut5->Draw("hist same");

  string plot8 = plottitlecutpt + string(".png");
  can1->Print(plot8.data(),"png");
 
     //Draw Legend for Cut histos
  auto legend2 = new TLegend(10, 10, 10, 10);
   legend2->SetHeader("Cut Values","R"); // option "C" allows to center the header
   legend2->AddEntry(hL1PtCut1,"Cut: 0.0-0.7","l");
   legend2->AddEntry(hL1PtCut2,"Cut: 0.0-0.6","l");
   legend2->AddEntry(hL1PtCut3,"Cut: 0.0-0.5","l");
   legend2->AddEntry(hL1PtCut4,"Cut: 0.0-0.4","l");
   legend2->AddEntry(hL1PtCut5,"Cut: 0.0-0.3","l");

   legend2->Draw();


   /*
   //Drawing Stack using Canvas
   TCanvas *cs = new TCanvas("cs","cs",10,10,500,700);
   TText T;
   hsEta->Draw();
   T.DrawTextNDC(.5,.95, "EtaCuts");
   hsEta->GetXaxis()->SetTitle("Electron #etacut (rad)");

   TCanvas *cs2 = new TCanvas("cs","cs",10,10,500,700);
   TText T2;
   hsPt->Draw();
   T2.DrawTextNDC(.5,.95, "PtCuts");
   hsPt->GetXaxis()->SetTitle("Electron #Ptcut (MeV)");
   */

    // Draw efficiency plot
  TGraphAsymmErrors tgEle;
  tgEle.Divide(hTrm,hTru, "cl=0.683 b(1,1) mode");

  can2->SetLeftMargin(0.16);
  can2->SetBottomMargin(0.16);

  // Title
  tgEle.GetXaxis()->SetTitle(xtitle.data());
  tgEle.GetYaxis()->SetTitle("Efficiency");
  tgEle.SetMarkerStyle(20);
  tgEle.SetMarkerColor(46);
  tgEle.GetXaxis()->SetTitleOffset(1.2);
  tgEle.GetYaxis()->SetTitleOffset(1.2);

  can2->cd();
  tgEle.Draw("AP");

  string plot4 = plottitle + string("Eff.png");
  can2->Print(plot4.data(),"png");

 TGraphAsymmErrors tgEleEt;
  tgEleEt.Divide(hTrmet,hTruet, "cl=0.683 b(1,1) mode");

  can2->SetLeftMargin(0.16);
  can2->SetBottomMargin(0.16);

  // Title
  tgEleEt.GetXaxis()->SetTitle(xtitleet.data());
  tgEleEt.GetYaxis()->SetTitle("Efficiency");
  tgEleEt.SetMarkerStyle(20);
  tgEleEt.SetMarkerColor(46);
  tgEleEt.GetXaxis()->SetTitleOffset(1.2);
  tgEleEt.GetYaxis()->SetTitleOffset(1.2);

  can2->cd();
  tgEleEt.Draw("AP");

  string plot5 = plottitleet + string("Effet.png");
  can2->Print(plot5.data(),"png");

 TGraphAsymmErrors tgElePhi;
  tgElePhi.Divide(hTrmphi,hTruphi, "cl=0.683 b(1,1) mode");

  can2->SetLeftMargin(0.16);
  can2->SetBottomMargin(0.16);

  // Title
  tgElePhi.GetXaxis()->SetTitle(xtitlephi.data());
  tgElePhi.GetYaxis()->SetTitle("Efficiency");
  tgElePhi.SetMarkerStyle(20);
  tgElePhi.SetMarkerColor(46);
  tgElePhi.GetXaxis()->SetTitleOffset(1.2);
  tgElePhi.GetYaxis()->SetTitleOffset(1.2);

  can2->cd();
  tgEleEt.Draw("AP");

  string plot6 = plottitlephi + string("Effphi.png");
  can2->Print(plot6.data(),"png");



}
