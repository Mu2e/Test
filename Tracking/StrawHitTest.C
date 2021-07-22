#include <algorithm>
#include "TTree.h"
#include "TCut.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TLine.h"
#include "TBox.h"
#include "TColor.h"
#include "TStyle.h"
#include <iostream>
using namespace std;

Int_t myhpart(Int_t mcpdg,Int_t mcproc){
  enum hpart{Hadron=0,LowEe,Mudecay,CE,Other};
  if(mcproc==167)return CE;
  if(mcpdg==11&&(mcproc==14||mcproc==114))return Mudecay;
  if(abs(mcpdg)==11&&(mcproc!=14&&mcproc!=114&&mcproc<165))return LowEe;
  if(mcpdg>2000)return Hadron;
  return Other;
}

void StrawHitTest (TTree* hits, const char* page="bcan",unsigned nevents=1000 ) {

  TString spage(page);
//  gStyle->SetOptStat(0);
  TCut electron("abs(mcpdg)==11");
  TCut conv("mcpdg==11&&mcproc==167");
  TCut oele("abs(mcpdg)==11&&(mcproc!=12&&mcproc!=14&&mcproc!=114&&mcproc<166)");
  TCut tgtdio("mcpdg==11&&(mcproc==166)");
  TCut ootdio("mcpdg==11&&(mcproc==14||mcproc==114)");
  TCut bkg("mcpdg==&&mcproc==17");
  TCut pconv("mcpdg==11&&mcproc==11");
  TCut compt("mcpdg==11&&mcproc==12");
  TCut proton("mcpdg==2212");
  TCut deuteron("mcpdg==1000010020");
  TCut hadron("mcpdg>2000");
  TCut photon("mcpdg==22");
  TCut neutron("mcpdg==2112");
  TCut bkge("abs(mcpdg)==11&&mcppdg==22");
  TCut bkgo("abs(mcpdg)!=11&&mcpdg!=2212");
  TCut xtalk("mcxtalk");
  TCut direct("!mcxtalk");
  TCut opart =(!electron)&&(!hadron);

  TCut highEe = conv||tgtdio||ootdio;

  TCut pproton("mcppdg==2212");
  TCut pdeuteron("mcppdg==1000010020");
  TCut pneuton("mcppdg==2112");
  TCut pele("mcppdg==11");
  TCut ppos("mcppdg===11");
  TCut pgam("mcppdg==22");
  TCut nopar("mcppdg==0");

  TCut dioorigin("mcppdg==13&&mcproc==166");
  TCut ootmuonorigin("mcppdg==13&&mcproc<167");
  TCut norigin("mcppdg==2112");
  TCut porigin("mcppdg==22");
  TCut stpprotonorigin("mcppdg==2212");
  TCut stpdeuteronorigin("mcppdg==1000010020");
  TCut pprotonorigin("mcppdg==2212");
  TCut flashstraw("plane>33&&straw>=90");

  TCut hitsel("esel&&rsel&&tsel&&(!bkg)");

  TCut goodevt("mcom>100");
  TCut goodpeak("abs(tpeak-mct0-25)<30");
  if(spage =="hits"){
    TCanvas* hcan = new TCanvas("hcan","Hit Properties",800,800);
    hcan->Divide(2,2);
    TH2F* eve = new TH2F("eve","StrawHit Edep vs #Sigma MC Edep;MeV;MeV",100,0,0.006,100,0,0.006);
    eve->SetStats(0);
    TH2F* etot = new TH2F("etot","Electron TOT vs MC Transverse Drift Distance;True Drift Distance (mm);TOT (ns)",50,0,2.5,16,0,64);
    etot->SetStats(0);
    TH2F* tdiv = new TH2F("tdiv","Electron Longitudinal Position;True Longitudinal Position (mm);#Delta T Longitudinal Position",50,-600,600,50,-600,600);
    tdiv->SetStats(0);
    TH1F* tdres = new TH1F("tdres","Electron time Division Resolution",100,-400,400);
    hits->Project("tdres","shlen-mcshlen",highEe);
    hits->Project("etot","0.5*(totcal+tothv):abs(mcshd)",highEe);
    hcan->cd(1);
    etot->Draw("colorz");
    hcan->cd(2);
    tdres->Fit("gaus");
    hcan->cd(3);
    hits->Draw("shlen:mcshlen>>tdiv",highEe,"colorz",100000);
    hcan->cd(4);
    hits->Draw("edep:mcedep>>eve","","colorz",100000);
    
  } else if(spage =="tcan"){
    THStack* tstack = new THStack("tc","Reco Hit Time by Particle;Hit Time (ns);Hits/event/ns");
    TH1F* ctime = new TH1F("ctime","Conversion Electron Reco Hit Time",300,250,1750);
    TH1F* ptime = new TH1F("ptime","Hadron Reco Hit Time",300,250,1750);
    TH1F* tgtetime = new TH1F("tgtetime","Target DIO Electron Reco Hit Time",300,250,1750);
    TH1F* ootetime = new TH1F("ootetime","OOT DIO Electron Reco Hit Time",300,250,1750);
    TH1F* cetime = new TH1F("cetime","Compton Electron Reco Hit Time",300,250,1750);
    TH1F* oetime = new TH1F("oetime","Other Electron Reco Hit Time",300,250,1750);
    TH1F* otime = new TH1F("otime","Other Particle Reco Hit Time",300,250,1750);
    ctime->SetFillColor(kRed);
    ptime->SetFillColor(kBlack);
    tgtetime->SetFillColor(kCyan);
    ootetime->SetFillColor(kOrange);
    cetime->SetFillColor(kGreen);
    oetime->SetFillColor(kBlue);
    otime->SetFillColor(kOrange);

    double scale = 1.0/(ctime->GetBinWidth(1)*nevents);
    double factor = (ctime->GetBinWidth(1));
    hits->Project("otime","ctime",opart);
    otime->Scale(scale);
    tstack->Add(otime);
    hits->Project("ctime","ctime",conv);
    ctime->Scale(scale);
    tstack->Add(ctime);
    hits->Project("ptime","ctime",hadron);
    ptime->Scale(scale);
    tstack->Add(ptime);
    hits->Project("tgtetime","ctime",tgtdio);
    tgtetime->Scale(scale);
    tstack->Add(tgtetime);
    hits->Project("ootetime","ctime",ootdio);
    ootetime->Scale(scale);
    tstack->Add(ootetime);
    hits->Project("oetime","ctime",oele);
    oetime->Scale(scale);
    tstack->Add(oetime);
    hits->Project("cetime","ctime",compt);
    cetime->Scale(scale);
    tstack->Add(cetime);

    double total = otime->Integral() + ctime->Integral() + ptime->Integral() + tgtetime->Integral()+ ootetime->Integral() + cetime->Integral() + oetime->Integral();
 
    cout << "All other integral = " << otime->Integral() 
    << " CE integral = " << ctime->Integral()
    << " P integral = " << ptime->Integral()
    << " Tgt DIO e integral = " << tgtetime->Integral() 
    << " OOT DIO e integral = " << ootetime->Integral() 
    << " Compton e integral = " << cetime->Integral() 
    << " Other e integral = " << oetime->Integral() 
    << " Other particle integral = " << otime->Integral() 
    << " Total = " << total << endl;

    TLegend* tleg = new TLegend(.6,.7,.9,.9);
    char title[50];
    snprintf(title,50,"CE, #int=%4.0f",ctime->Integral()*factor);
    tleg->AddEntry(ctime,title,"F");
    snprintf(title,50,"Hadron, #int=%4.0f",ptime->Integral()*factor);
    tleg->AddEntry(ptime,title,"F");
    snprintf(title,50,"Target DIO e^{-}, #int=%4.0f",tgtetime->Integral()*factor);
    tleg->AddEntry(tgtetime,title,"F");
    snprintf(title,50,"OOT DIO e^{-}, #int=%4.0f",ootetime->Integral()*factor);
    tleg->AddEntry(ootetime,title,"F");
    snprintf(title,50,"Other e^{#pm}, #int=%4.0f",oetime->Integral()*factor);
    tleg->AddEntry(oetime,title,"F");
    snprintf(title,50,"Compton e^{-}, #int=%4.0f",cetime->Integral()*factor);
    tleg->AddEntry(cetime,title,"F");
    snprintf(title,50,"Other Particles, #int=%4.0f",otime->Integral()*factor);
    tleg->AddEntry(otime,title,"F");
    snprintf(title,50,"Total #int=%4.0f",total*factor);
    tleg->AddEntry((TObject*)0,title,"");

    TCanvas* tcan = new TCanvas("tcan","tcan",1200,800);
    tcan->Divide(2,1);
    tcan->cd(1);
    tstack->Draw("h");
    tleg->Draw();

    THStack* tstacks = new THStack("tcs","Reco Hit Time by Particle;Hit Time (ns);Hits/event/ns");
    TH1F* ctimes = new TH1F("ctimes","Conversion Electron Reco Hit Time",300,250,1750);
    TH1F* ptimes = new TH1F("ptimes","Hadron Reco Hit Time",300,250,1750);
    TH1F* tgtetimes = new TH1F("tgtetimes","Target DIO Electron Reco Hit Time",300,250,1750);
    TH1F* ootetimes = new TH1F("ootetimes","OOT DIO Electron Reco Hit Time",300,250,1750);
    TH1F* cetimes = new TH1F("cetimes","Compton Electron Reco Hit Time",300,250,1750);
    TH1F* oetimes = new TH1F("oetimes","Other Electron Reco Hit Time",300,250,1750);
    TH1F* otimes = new TH1F("otimes","Other Particle Reco Hit Time",300,250,1750);
    ctimes->SetFillColor(kRed);
    ptimes->SetFillColor(kBlack);
    tgtetimes->SetFillColor(kCyan);
    ootetimes->SetFillColor(kOrange);
    cetimes->SetFillColor(kGreen);
    oetimes->SetFillColor(kBlue);
    otimes->SetFillColor(kOrange);

    hits->Project("otimes","ctime",hitsel+opart);
    otimes->Scale(scale);
    tstacks->Add(otimes);
    hits->Project("ctimes","ctime",hitsel+conv);
    ctimes->Scale(scale);
    tstacks->Add(ctimes);
    hits->Project("ptimes","ctime",hitsel+hadron);
    ptimes->Scale(scale);
    tstacks->Add(ptimes);
    hits->Project("tgtetimes","ctime",hitsel+tgtdio);
    tgtetimes->Scale(scale);
    tstacks->Add(tgtetimes);
    hits->Project("ootetimes","ctime",hitsel+ootdio);
    ootetimes->Scale(scale);
    tstacks->Add(ootetimes);
    hits->Project("oetimes","ctime",hitsel+oele);
    oetimes->Scale(scale);
    tstacks->Add(oetimes);
    hits->Project("cetimes","ctime",hitsel+compt);
    cetimes->Scale(scale);
    tstacks->Add(cetimes);

    double totals = otimes->Integral() + ctimes->Integral() + ptimes->Integral() + tgtetimes->Integral()+ ootetimes->Integral() + cetimes->Integral() + oetimes->Integral();
 
    cout << "Selected All other integral = " << otimes->Integral() 
    << " CE integral = " << ctimes->Integral()
    << " P integral = " << ptimes->Integral()
    << " Tgt DIO e integral = " << tgtetimes->Integral() 
    << " OOT DIO e integral = " << ootetimes->Integral() 
    << " Compton e integral = " << cetimes->Integral() 
    << " Other e integral = " << oetimes->Integral() 
    << " Other particle integral = " << otimes->Integral() 
    << " Total = " << totals << endl;

    TLegend* tlegs = new TLegend(.6,.7,.9,.9);
    snprintf(title,50,"CE, #int=%4.0f",ctimes->Integral()*factor);
    tlegs->AddEntry(ctime,title,"F");
    snprintf(title,50,"Hadron, #int=%4.0f",ptimes->Integral()*factor);
    tlegs->AddEntry(ptime,title,"F");
    snprintf(title,50,"Target DIO e^{-}, #int=%4.0f",tgtetimes->Integral()*factor);
    tlegs->AddEntry(tgtetime,title,"F");
    snprintf(title,50,"OOT DIO e^{-}, #int=%4.0f",ootetimes->Integral()*factor);
    tlegs->AddEntry(ootetime,title,"F");
    snprintf(title,50,"Other e^{#pm}, #int=%4.0f",oetimes->Integral()*factor);
    tlegs->AddEntry(oetime,title,"F");
    snprintf(title,50,"Compton e^{-}, #int=%4.0f",cetimes->Integral()*factor);
    tlegs->AddEntry(cetime,title,"F");
    snprintf(title,50,"Other Particles, #int=%4.0f",otimes->Integral()*factor);
    tlegs->AddEntry(otime,title,"F");
    snprintf(title,50,"Total #int=%4.0f",totals*factor);
    tlegs->AddEntry((TObject*)0,title,"");


    tcan->cd(2);
    tstacks->Draw("h");
    tlegs->Draw();

  } else if(spage=="particle"){
//    gStyle->SetOptStat(111111);
    THStack* estack = new THStack("edep","Reco Hit Energy by Particle;Deposited Energy (KeV)");
    TH1F* econv = new TH1F("econv","Straw Hit Energy;Deposited Energy (KeV)",200,-1.0,6.0);
    estack->Add(econv);
    TH1F* emu = new TH1F("emu","Straw Hit Energy;Deposited Energy (KeV)",200,-1.0,6.0);
    estack->Add(emu);
    TH1F* egam = new TH1F("egam","Straw Hit Energy;Deposited Energy (KeV)",200,-1.0,6.0);
    estack->Add(egam);
    TH1F* ehad = new TH1F("ehad","Straw Hit Energy;Deposited Energy (KeV)",200,-1.0,6.0);
    estack->Add(ehad);
    econv->SetFillColor(kRed);
    emu->SetFillColor(kCyan);
    egam->SetFillColor(kGreen);
    ehad->SetFillColor(kMagenta);
//    econv->SetStats(0);
    emu->SetStats(0);
    egam->SetStats(0);
    ehad->SetStats(0);

    THStack* rstack = new THStack("rho","Reco Hit Radius by Particle;Transverse Radius (mm)");
    TH1F* rconv = new TH1F("rconv","Straw Hit Radius;Transverse Radius (mm)",100,360,800);
    rstack->Add(rconv);
    TH1F* rmu = new TH1F("rmu","Straw Hit Radius;Transverse Radius (mm)",100,360,800);
    rstack->Add(rmu);
    TH1F* rgam = new TH1F("rgam","Straw Hit Radius;Transverse Radius (mm)",100,360,800);
    rstack->Add(rgam);
    TH1F* rhad = new TH1F("rhad","Straw Hit Radius;Transverse Radius (mm)",100,360,800);
    rstack->Add(rhad);
//    TH1F* rx = new TH1F("rx","Straw Hit Radius;Transverse Radius (mm)",100,360,800);
    rconv->SetFillColor(kRed);
    rmu->SetFillColor(kCyan);
    rgam->SetFillColor(kGreen);
    rhad->SetFillColor(kMagenta);
    hits->Project("econv","edep*1000.0",conv+direct);
    hits->Project("emu","edep*1000.0",tgtdio+direct);
    hits->Project("egam","edep*1000.0",bkge+direct);
    hits->Project("ehad","edep*1000.0",hadron+direct);

    hits->Project("rconv","sqrt(shpos.data[1]^2+shpos.data[0]^2)",conv+direct);
    hits->Project("rmu","sqrt(shpos.data[1]^2+shpos.data[0]^2)",tgtdio+direct);
    hits->Project("rgam","sqrt(shpos.data[1]^2+shpos.data[0]^2)",bkge+direct);
    hits->Project("rhad","sqrt(shpos.data[1]^2+shpos.data[0]^2)",hadron+direct);
    
    TCanvas* bcan = new TCanvas("bcan","background",1000,800);
    bcan->Divide(1,2);
    bcan->cd(1);
    gPad->SetLogy();
    estack->Draw();
    econv->Draw("same");
    TLegend* leg2 = new TLegend(0.55,0.7,0.9,0.9);
    leg2->AddEntry(econv,"Conversion electron","f");
    leg2->AddEntry(emu,"#mu#rightarrowe#nu#nu","f");
    leg2->AddEntry(egam,"Low-E e","f");
    leg2->AddEntry(ehad,"Hadron","f");
    leg2->Draw();

    bcan->cd(2);
 //   gPad->SetLogy();
    rstack->Draw();
  } else if(spage=="sorigin"){

    THStack* sorigin = new THStack("sorigin","Selected Hit Time by Generator Particle;Hit Time (ns);Hits/event/ns");
    TH1F* dtime = new TH1F("dtime","DIO Reco Hit Time",300,250,1750);
    TH1F* gtime = new TH1F("gtime","Photon Reco Hit Time",300,250,1750);
    TH1F* pptime = new TH1F("pptime","Primary Proton Reco Hit Time",300,250,1750);
    TH1F* stptime = new TH1F("stptime","Stopping Target Proton Reco Hit Time",300,250,1750);
    TH1F* stdtime = new TH1F("stdtime","Stopping Target Deuteron Reco Hit Time",300,250,1750);
    TH1F* ntime = new TH1F("ntime","Neutron Reco Hit Time",300,250,1750);
    TH1F* mtime = new TH1F("mtime","OOT Muon Reco Hit Time",300,250,1750);
    TH1F* cetime = new TH1F("cetime","CE Reco Hit Time",300,250,1750);
    dtime->SetFillColor(kOrange);
    gtime->SetFillColor(kBlack);
    pptime->SetFillColor(kBlue);
    stptime->SetFillColor(kGreen);
    stdtime->SetFillColor(kMagenta);
    ntime->SetFillColor(kCyan);
    mtime->SetFillColor(kYellow);
    cetime->SetFillColor(kRed);

    double scale = 1.0/(dtime->GetBinWidth(1)*nevents);
    hits->Project("cetime","ctime",conv+hitsel);
    cetime->Scale(scale);
    sorigin->Add(cetime);
     hits->Project("dtime","ctime",dioorigin+hitsel);
    dtime->Scale(scale);
    sorigin->Add(dtime);
    hits->Project("stdtime","ctime",stpdeuteronorigin+hitsel);
    stdtime->Scale(scale);
    sorigin->Add(stdtime);
    hits->Project("mtime","ctime",ootmuonorigin+hitsel);
    mtime->Scale(scale);
    sorigin->Add(mtime);
    hits->Project("gtime","ctime",porigin+hitsel);
    gtime->Scale(scale);
    sorigin->Add(gtime);
    hits->Project("stptime","ctime",stpprotonorigin+hitsel);
    stptime->Scale(scale);
    sorigin->Add(stptime);
    hits->Project("pptime","ctime",pprotonorigin+hitsel);
    pptime->Scale(scale);
    sorigin->Add(pptime);
    hits->Project("ntime","ctime",norigin+hitsel);
    ntime->Scale(scale);
    sorigin->Add(ntime);


    double total = dtime->Integral() + pptime->Integral() + stptime->Integral() + gtime->Integral() + ntime->Integral() + mtime->Integral() + cetime->Integral();

    cout << "DIO integral = " << dtime->Integral()
      << " Primary Proton inegral = " << pptime->Integral()
      << " ST Proton inegral = " << stptime->Integral()
      << " Photon inegral = " << gtime->Integral()
      << " Neutron inegral = " << ntime->Integral()
      << " OOT muon inegral = " << mtime->Integral()
      << " CE inegral = " << cetime->Integral()
      << " total = " << total << endl;

    TCanvas* socan = new TCanvas("socan","socan",800,800);
    sorigin->Draw("h");
    TLegend* tleg = new TLegend(.6,.6,.9,.9);
    char title[50];
    double factor = dtime->GetBinWidth(1);
    snprintf(title,50,"Neutron, #int=%4.0f",ntime->Integral()*factor);
    tleg->AddEntry(ntime,title,"F");
    snprintf(title,50,"Primary Proton, #int=%4.0f",pptime->Integral()*factor);
    tleg->AddEntry(pptime,title,"F");
    snprintf(title,50,"Stopping Target Proton, #int=%4.0f",stptime->Integral()*factor);
    tleg->AddEntry(stptime,title,"F");
    snprintf(title,50,"Photon, #int=%4.0f",gtime->Integral()*factor);
    tleg->AddEntry(gtime,title,"F");
    snprintf(title,50,"OOT Muon, #int=%4.0f",mtime->Integral()*factor);
    tleg->AddEntry(mtime,title,"F");
    snprintf(title,50,"Stopping Target Deuteron, #int=%4.0f",stdtime->Integral()*factor);
    tleg->AddEntry(stdtime,title,"F");
    snprintf(title,50,"DIO, #int=%4.0f",dtime->Integral()*factor);
    tleg->AddEntry(dtime,title,"F");
    snprintf(title,50,"CE, #int=%4.0f",cetime->Integral()*factor);
    tleg->AddEntry(cetime,title,"F");
    snprintf(title,50,"Total #int=%4.0f",total*factor);
    tleg->AddEntry(dtime,title,"");
    tleg->Draw();
  } else if(spage=="hitsel"){
    TH2F* hsel = new TH2F("hsel","Hit Selection;Producing Particle;Cut efficiency (%)",5,-0.5,4.5,4,-0.5,3.5);
    TAxis* yax = hsel->GetYaxis();
    unsigned ibin(1);
    yax->SetBinLabel(ibin++,"Hit Time");
    yax->SetBinLabel(ibin++,"Hit Energy");
    yax->SetBinLabel(ibin++,"Hit Radius");
    yax->SetBinLabel(ibin++,"Bkg Hit");
    TAxis* xax = hsel->GetXaxis();
    ibin = 1;
    xax->SetBinLabel(ibin++,"Hadron");
    xax->SetBinLabel(ibin++,"LowEe");
    xax->SetBinLabel(ibin++,"#mu#rightarrowe#nu#nu");
    xax->SetBinLabel(ibin++,"Ce");
    xax->SetBinLabel(ibin++,"Other");
    // first, get normalization
    TH1F* myhp = new TH1F("myhp","Hit Rate;Producing Particle;Hits/event",5,-0.5,4.5);
    TH1F* myhpg = new TH1F("myhpg","Hit Rate;Producing Particle;Hits/event",5,-0.5,4.5);
    xax = myhp->GetXaxis();
    ibin = 1;
    xax->SetBinLabel(ibin++,"Hadron");
    xax->SetBinLabel(ibin++,"LowEe");
    xax->SetBinLabel(ibin++,"#mu#rightarrowe#nu#nu");
    xax->SetBinLabel(ibin++,"Ce");
    xax->SetBinLabel(ibin++,"Other");
    hits->Project("myhp","myhpart(mcpdg,mcproc)","tsel");
    hits->Project("myhpg","myhpart(mcpdg,mcproc)",hitsel);
    double pscale(1.0/nevents);
    myhp->Scale(pscale);
    myhpg->Scale(pscale);
    myhpg->SetFillColor(kGreen);
    // now loop over selections
    std::vector<TCut> selcuts = {"tsel","tsel&&esel","tsel&&esel&&rsel",hitsel};
    for(size_t icut=0;icut< selcuts.size();++icut){
      char val[100];
      cout << "Projecting cut " << selcuts[icut] << endl;
      snprintf(val,100,"%lu:myhpart(mcpdg,mcproc)",icut);
      hits->Project("+hsel",val,selcuts[icut]);
    }
// normalize by row
    for(int ibin=1;ibin <= myhp->GetXaxis()->GetNbins();++ibin){
      double norm = 100.0/hsel->GetBinContent(ibin,1);
      cout << "Normalization for " << hsel->GetYaxis()->GetBinLabel(ibin)  << " = " << norm << endl;
      for(int jbin=1;jbin <= hsel->GetYaxis()->GetNbins(); ++jbin) {
	double val =hsel->GetBinContent(ibin,jbin);
	cout << "value for ibin " << ibin <<" jbin " << jbin << " val " << val << endl;
	hsel->SetBinContent(ibin,jbin,val*norm);
      }
    }
    TCanvas* hscan = new TCanvas("hscan","hscan",750,750);
    hscan->Divide(1,2);
    hscan->cd(1);
    hsel->Draw("boxtext0");
    hscan->cd(2);
    TLegend* leg = new TLegend(0.6,0.7,0.8,0.9);
    leg->AddEntry(myhp,"All Hits","l");
    leg->AddEntry(myhpg,"Selected Hits","f");
    myhp->Draw("histtext0");
    myhpg->Draw("histtext90same");
    leg->Draw();
  }
  
}
