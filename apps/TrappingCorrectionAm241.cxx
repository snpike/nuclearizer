/* 
 * TrappingCorrectionAm241.cxx
 *
 *
 * Copyright (C) by Sean Pike.
 * All rights reserved.
 *
 *
 * This code implementation is the intellectual property of
 * Sean Pike.
 *
 * By copying, distributing or modifying the Program (or any work
 * based on the Program) you indicate your acceptance of this statement,
 * and all its terms.
 *
 */

// Standard
#include <iostream>
#include <string>
#include <sstream>
#include <csignal>
#include <cstdlib>
#include <map>
#include <vector>
#include <algorithm>
#include <ctime>
#include <cstdio>
using namespace std;

// ROOT
#include <TROOT.h>
#include <TEnv.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph2D.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TStopwatch.h>
#include <TProfile.h>

// MEGAlib
#include "MGlobal.h"
#include "MFile.h"
#include "MReadOutElementDoubleStrip.h"
#include "MFileReadOuts.h"
#include "MReadOutAssembly.h"
#include "MStripHit.h"
#include "MReadOutSequence.h"
#include "MSupervisor.h"
#include "MModuleLoaderMeasurementsHDF.h"
#include "MModuleEnergyCalibrationUniversal.h"
#include "MModuleEventFilter.h"
#include "MModuleStripPairingGreedy.h"
#include "MModuleStripPairingChiSquare.h"
#include "MModuleTACcut.h"
#include "MAssembly.h"


double g_MinCTD = -300;
double g_MaxCTD = 300;
int g_MinCounts = 1000;

////////////////////////////////////////////////////////////////////////////////


//! A standalone program based on MEGAlib and ROOT
class TrappingCorrectionAm241
{
public:
  //! Default constructor
  TrappingCorrectionAm241();
  //! Default destructor
  ~TrappingCorrectionAm241();
  
  //! Parse the command line
  bool ParseCommandLine(int argc, char** argv);
  //! Analyze what ever needs to be analyzed...
  bool Analyze();
  //! Interrupt the analysis
  void Interrupt() { m_Interrupt = true; }

  MStripHit* GetDominantStrip(vector<MStripHit*>& Strips, double& EnergyFraction);

private:
  //! True, if the analysis needs to be interrupted
  bool m_Interrupt;
  //! The input file name
  MString m_HVFileName;
  MString m_LVFileName;
  MString m_EcalFile;
  MString m_TACCalFile;
  MString m_TACCutFile;
  MString m_StripMapFile;
  //! output file names
  MString m_OutFile;
  //! option to do a pixel-by-pixel calibration (instead of detector-by-detector)
  bool m_PixelCorrect;
  bool m_GreedyPairing;
  bool m_ExcludeNN;

  double m_MinEnergy;
  double m_MaxEnergy;

};

////////////////////////////////////////////////////////////////////////////////


//! Default constructor
TrappingCorrectionAm241::TrappingCorrectionAm241() : m_Interrupt(false)
{
  gStyle->SetPalette(1, 0);
}


////////////////////////////////////////////////////////////////////////////////


//! Default destructor
TrappingCorrectionAm241::~TrappingCorrectionAm241()
{
  // Intentionally left blank
}


////////////////////////////////////////////////////////////////////////////////


//! Parse the command line
bool TrappingCorrectionAm241::ParseCommandLine(int argc, char** argv)
{
  ostringstream Usage;
  Usage<<endl;
  Usage<<"  Usage: TrappingCorrectionAm241 <options>"<<endl;
  Usage<<"    General options:"<<endl;
  Usage<<"         --HVfile:   HV illumination input file name (.hdf5 or .txt with list of hdf5s)"<<endl;
  Usage<<"         --LVfile:   LV illumination input file name (.hdf5 or .txt with list of hdf5s)"<<endl;
  Usage<<"         --emin:   minimum Event energy (default 40 keV)"<<endl;
  Usage<<"         --emax:   maximum Event energy (default 70 kev)"<<endl;
  Usage<<"         -e:   energy calibration file (.ecal)"<<endl;
  Usage<<"         --tcal:   TAC calibration file"<<endl;
  Usage<<"         --tcut:   TAC cut file"<<endl;
  Usage<<"         -p:   do pixel-by-pixel correction"<<endl;
  Usage<<"         -m:   strip map file name (.map)"<<endl;
  Usage<<"         -g:   greedy strip pairing (default is chi-square)"<<endl;
  Usage<<"         -n:   exclude nearest neighbors"<<endl;
  Usage<<"         -o:   outfile (default YYYYMMDDHHMMSS)"<<endl;
  Usage<<"         -h:   print this help"<<endl;
  Usage<<endl;

  string Option;

  // Check for help
  for (int i = 1; i < argc; i++) {
    Option = argv[i];
    if (Option == "-h" || Option == "--help" || Option == "?" || Option == "-?") {
      cout<<Usage.str()<<endl;
      return false;
    }
  }

  m_PixelCorrect = false;
  m_GreedyPairing = false;
  m_MinEnergy = 40;
  m_MaxEnergy = 70;
  
  time_t rawtime;
  tm* timeinfo;
  time(&rawtime);
  timeinfo = gmtime(&rawtime);

  char buffer [80];
  strftime(buffer,80,"%Y%m%d%H%M%S",timeinfo);

  m_OutFile = MString(buffer);

  // Now parse the command line options:
  for (int i = 1; i < argc; i++) {
    Option = argv[i];

    // First check if each option has sufficient arguments:
    // Single argument
    if ((Option == "--HVfile") || (Option == "--LVfile") || (Option == "-o") || (Option == "--emin") || (Option == "--emax") || (Option == "--tcal") || (Option == "--tcut") || (Option == "-m")) {
      if (!((argc > i+1) && (argv[i+1][0] != '-' || isalpha(argv[i+1][1]) == 0))){
        cout<<"Error: Option "<<argv[i][1]<<" needs a second argument!"<<endl;
        cout<<Usage.str()<<endl;
        return false;
      }
    } 

    // Then fulfill the options:
    if (Option == "--HVfile") {
      m_HVFileName = argv[++i];
      cout<<"Accepting file name: "<<m_HVFileName<<endl;
    } 

    if (Option == "--LVfile") {
      m_LVFileName = argv[++i];
      cout<<"Accepting file name: "<<m_LVFileName<<endl;
    } 

    if (Option == "-e") {
      m_EcalFile = argv[++i];
      cout<<"Accepting file name: "<<m_EcalFile<<endl;
    } 

    if (Option == "--emin") {
      m_MinEnergy = stod(argv[++i]);
    } 

    if (Option == "--emax") {
      m_MaxEnergy = stod(argv[++i]);
    } 

    if (Option == "--tcal") {
      m_TACCalFile = argv[++i];
      cout<<"Accepting file name: "<<m_TACCalFile<<endl;
    } 

    if (Option == "--tcut") {
      m_TACCutFile = argv[++i];
      cout<<"Accepting file name: "<<m_TACCutFile<<endl;
    }

    if (Option == "-m"){
      m_StripMapFile = argv[++i];
      cout<<"Accepting file name: "<<m_StripMapFile<<endl;
    } 

    if (Option == "-o"){
      m_OutFile = argv[++i];
      cout<<"Accepting file name: "<<m_OutFile<<endl;
    }

    if (Option == "-p"){
      m_PixelCorrect = true;
    }

    if (Option == "-g"){
      m_GreedyPairing = true;
    }

    if (Option == "-n"){
      m_ExcludeNN = true;
    }

  }

  return true;
}


////////////////////////////////////////////////////////////////////////////////


//! Do whatever analysis is necessary
bool TrappingCorrectionAm241::Analyze()
{
  //time code just to see
  TStopwatch watch;
  watch.Start();

  if (m_Interrupt == true) return false;

  // [CTD, HV Energy, LV Energy] for HV illumination and LV illumination
  map<unsigned int, vector<vector<double>>> Endpoints;

  vector<MString> FileNames;
  FileNames.push_back(m_HVFileName);
  FileNames.push_back(m_LVFileName);

  vector<MString> IllumSide;
  IllumSide.push_back(MString("HV"));
  IllumSide.push_back(MString("LV"));

  vector<map<int, TH1D*>> CTDHistograms;
  map<int, TH1D*> tempCTDHistHV;
  map<int, TH1D*> tempCTDHistLV;
  CTDHistograms.push_back(tempCTDHistHV);
  CTDHistograms.push_back(tempCTDHistLV);

  vector<map<int, TH1D*>> HVEnergyHistograms;
  map<int, TH1D*> tempHVHistHV;
  map<int, TH1D*> tempHVHistLV;
  HVEnergyHistograms.push_back(tempHVHistHV);
  HVEnergyHistograms.push_back(tempHVHistLV);

  vector<map<int, TH1D*>> LVEnergyHistograms;
  map<int, TH1D*> tempLVHistHV;
  map<int, TH1D*> tempLVHistLV;
  LVEnergyHistograms.push_back(tempLVHistHV);
  LVEnergyHistograms.push_back(tempLVHistLV);


  for (unsigned int s=0; s<2; ++s) {

    MString InputFile = FileNames[s];
    vector<MString> HDFNames;

    if ((InputFile.GetSubString(InputFile.Length() - 4)) == "hdf5") {
      HDFNames.push_back(InputFile);
    } else if ((InputFile.GetSubString(InputFile.Length() - 3)) == "txt") {
      MFile F;
      if (F.Open(InputFile)==false) {
        cout<<"Error: Failed to open input file."<<endl;
      } else {
        MString Line;
        while (F.ReadLine(Line)) {
          HDFNames.push_back(Line.Trim());
        }
      }
    }

    for (unsigned int f = 0; f<HDFNames.size(); ++f) {

      MString File = HDFNames[f];
      cout<<"Beginning analysis of file "<<File<<endl;

      MSupervisor* S = MSupervisor::GetSupervisor();
      
    	MModuleLoaderMeasurementsHDF* Loader;
    	MModuleTACcut* TACCalibrator;
    	MModuleEnergyCalibrationUniversal* EnergyCalibrator;
    	MModuleEventFilter* EventFilter;

      unsigned int MNumber = 0;
      cout<<"Creating HDF5 loader"<<endl;
      Loader = new MModuleLoaderMeasurementsHDF();
      Loader->SetFileNameStripMap(m_StripMapFile);
      Loader->SetFileName(File);
      Loader->SetLoadContinuationFiles(true);
      S->SetModule(Loader, MNumber);
      ++MNumber;

      cout<<"Creating TAC calibrator"<<endl;
      TACCalibrator = new MModuleTACcut();
      TACCalibrator->SetTACCalFileName(m_TACCalFile);
      TACCalibrator->SetTACCutFileName(m_TACCutFile);
      S->SetModule(TACCalibrator, MNumber);
      ++MNumber;
     
      cout<<"Creating energy calibrator"<<endl;
      EnergyCalibrator = new MModuleEnergyCalibrationUniversal();
      EnergyCalibrator->SetFileName(m_EcalFile);
      EnergyCalibrator->EnablePreampTempCorrection(false);
      S->SetModule(EnergyCalibrator, MNumber);
      ++MNumber;

      cout<<"Creating Event filter"<<endl;
      //! Only use events with 1 Strip Hit on each side to avoid strip pairing complications
      EventFilter = new MModuleEventFilter();
      EventFilter->SetMinimumLVStrips(1);
      EventFilter->SetMaximumLVStrips(3);
      EventFilter->SetMinimumHVStrips(1);
      EventFilter->SetMaximumHVStrips(3);
      EventFilter->SetMinimumHits(0);
      EventFilter->SetMaximumHits(100);
      EventFilter->SetMinimumTotalEnergy(m_MinEnergy);
      EventFilter->SetMaximumTotalEnergy(m_MaxEnergy*2);
      S->SetModule(EventFilter, MNumber);
      ++MNumber;
      
      cout<<"Creating strip pairing"<<endl;
      MModule* Pairing;
      if (m_GreedyPairing == true) {
        Pairing = new MModuleStripPairingGreedy();
      } else {
        Pairing = new MModuleStripPairingChiSquare();
      }
      S->SetModule(Pairing, MNumber);

      cout<<"Initializing Loader"<<endl;
      if (Loader->Initialize() == false) return false;
      cout<<"Initializing TAC calibrator"<<endl;
      if (TACCalibrator->Initialize() == false) return false;
      cout<<"Initializing Energy calibrator"<<endl;
      if (EnergyCalibrator->Initialize() == false) return false;
      cout<<"Initializing Event filter"<<endl;
      if (EventFilter->Initialize() == false) return false;
      cout<<"Initializing Pairing"<<endl;
      if (Pairing->Initialize() == false) return false;

      bool IsFinished = false;
      MReadOutAssembly* Event = new MReadOutAssembly();

      cout<<"Analyzing..."<<endl;
      while ((IsFinished == false) && (m_Interrupt == false)) {
        Event->Clear();

        if (Loader->IsReady()) {

          Loader->AnalyzeEvent(Event);
          TACCalibrator->AnalyzeEvent(Event);
          EnergyCalibrator->AnalyzeEvent(Event);
          bool Unfiltered = EventFilter->AnalyzeEvent(Event);

          if (Unfiltered == true) {
            
            Pairing->AnalyzeEvent(Event);

            if ((Event->HasAnalysisProgress(MAssembly::c_StripPairing) == true) && (Unfiltered==true) && (Event->GetNHits()==1)) {
              
              for (unsigned int h = 0; h < Event->GetNHits(); ++h) {
                double HVEnergy = 0.0;
                double LVEnergy = 0.0;
                double HVEnergyResolution = 0.0;
                double LVEnergyResolution = 0.0;
                vector<MStripHit*> HVStrips;
                vector<MStripHit*> LVStrips;

                MHit* H = Event->GetHit(h);
                
                int DetID = H->GetStripHit(0)->GetDetectorID();

                for (unsigned int sh = 0; sh < H->GetNStripHits(); ++sh) {
                  MStripHit* SH = H->GetStripHit(sh);

                  if ((m_ExcludeNN==false) || ((m_ExcludeNN==true) && (SH->IsNearestNeighbor()==false))) {
                    if (SH->IsLowVoltageStrip()==true) {
                      LVEnergy += SH->GetEnergy();
                      LVEnergyResolution += (SH->GetEnergyResolution())*(SH->GetEnergyResolution());
                      LVStrips.push_back(SH);
                    } else {
                      HVEnergy += SH->GetEnergy();
                      HVEnergyResolution += (SH->GetEnergyResolution())*(SH->GetEnergyResolution());
                      HVStrips.push_back(SH);
                    }
                  }
                }

                if ((HVStrips.size()>0) && (LVStrips.size()>0)) {
                  
                  double HVEnergyFraction = 0;
                  double LVEnergyFraction = 0;
                  MStripHit* HVSH = GetDominantStrip(HVStrips, HVEnergyFraction); 
                  MStripHit* LVSH = GetDominantStrip(LVStrips, LVEnergyFraction);
                  
                  if ((LVSH->HasCalibratedTiming()==true) && (HVSH->HasCalibratedTiming()==true)) {
                    
                    double CTD = LVSH->GetTiming() - HVSH->GetTiming();
                    
                    int PixelID = (10000*DetID) + (100*LVSH->GetStripID()) + (HVSH->GetStripID());

                    if (Endpoints.find(PixelID)==Endpoints.end()) {
                      vector<double> tempHVvec;
                      vector<double> tempLVvec;
                      vector<vector<double>> tempvec;
                      Endpoints[PixelID] = tempvec;
                      Endpoints[PixelID].push_back(tempHVvec);
                      Endpoints[PixelID].push_back(tempLVvec);
                    }

                    TH1D* CTDHist = CTDHistograms[s][PixelID];
                    TH1D* HVHist = HVEnergyHistograms[s][PixelID];
                    TH1D* LVHist = LVEnergyHistograms[s][PixelID];
                    
                    if (CTDHist == nullptr) {
                      char name[64]; sprintf(name,"CTD: PixelID %d %s Illumination",PixelID,IllumSide[s].Data());
                      CTDHist = new TH1D(name, name, (g_MaxCTD - g_MinCTD)/2, g_MinCTD, g_MaxCTD);
                      CTDHist->SetXTitle("CTD (ns)");
                      CTDHist->SetYTitle("Hits");
                      CTDHistograms[s][PixelID] = CTDHist;
                    }

                    if (HVHist == nullptr) {
                      char name[64]; sprintf(name,"HV energy: PixelID %d %s Illumination",PixelID,IllumSide[s].Data());
                      HVHist = new TH1D(name, name, (m_MaxEnergy - m_MinEnergy)*2, m_MinEnergy, m_MaxEnergy);
                      HVHist->SetXTitle("HV Energy (keV)");
                      HVHist->SetYTitle("Hits");
                      HVEnergyHistograms[s][PixelID] = HVHist;
                    }

                    if (LVHist == nullptr) {
                      char name[64]; sprintf(name,"LV energy: PixelID %d %s Illumination",PixelID,IllumSide[s].Data());
                      LVHist = new TH1D(name, name, (m_MaxEnergy - m_MinEnergy)*2, m_MinEnergy, m_MaxEnergy);
                      LVHist->SetXTitle("LV Energy (keV)");
                      LVHist->SetYTitle("Hits");
                      LVEnergyHistograms[s][PixelID] = LVHist;
                    }

                    CTDHist->Fill(CTD);
                    HVHist->Fill(HVEnergy);
                    LVHist->Fill(LVEnergy);
                  }
                }
              }
            }
          }
        }
        IsFinished = Loader->IsFinished();
      }
    }
  }

  map<unsigned int, vector<vector<double>>> FullDetEndpoints;

  vector<map<int, TH1D*>> FullDetCTDHistograms;
  map<int, TH1D*> FullDettempCTDHistHV;
  map<int, TH1D*> FullDettempCTDHistLV;
  FullDetCTDHistograms.push_back(FullDettempCTDHistHV);
  FullDetCTDHistograms.push_back(FullDettempCTDHistLV);

  vector<map<int, TH1D*>> FullDetHVEnergyHistograms;
  map<int, TH1D*> FullDettempHVHistHV;
  map<int, TH1D*> FullDettempHVHistLV;
  FullDetHVEnergyHistograms.push_back(FullDettempHVHistHV);
  FullDetHVEnergyHistograms.push_back(FullDettempHVHistLV);

  vector<map<int, TH1D*>> FullDetLVEnergyHistograms;
  map<int, TH1D*> FullDettempLVHistHV;
  map<int, TH1D*> FullDettempLVHistLV;
  FullDetLVEnergyHistograms.push_back(FullDettempLVHistHV);
  FullDetLVEnergyHistograms.push_back(FullDettempLVHistLV);


  for (unsigned int s=0; s<2; ++s) {

    for (auto H: CTDHistograms[s]) {
      
      int PixelID = H.first;
      int DetID = (PixelID-(PixelID%10000))/10000;
      int LVSHID = (PixelID-(DetID*10000) - (PixelID%100))/100;
      int HVSHID = (PixelID-(DetID*10000) - (LVSHID *100));

      if (FullDetEndpoints.find(DetID)==FullDetEndpoints.end()) {
        vector<double> tempHVvec;
        vector<double> tempLVvec;
        vector<vector<double>> tempvec;
        FullDetEndpoints[DetID] = tempvec;
        FullDetEndpoints[DetID].push_back(tempHVvec);
        FullDetEndpoints[DetID].push_back(tempLVvec);
      }

      TH1D* CTDHist = FullDetCTDHistograms[s][DetID];
      TH1D* HVHist = FullDetHVEnergyHistograms[s][DetID];
      TH1D* LVHist = FullDetLVEnergyHistograms[s][DetID];
      
      if (CTDHist == nullptr) {
        char name[64]; sprintf(name,"CTD: DetID %d %s Illumination",DetID,IllumSide[s].Data());
        CTDHist = new TH1D(name, name, (g_MaxCTD - g_MinCTD)/2, g_MinCTD, g_MaxCTD);
        CTDHist->SetXTitle("CTD (ns)");
        CTDHist->SetYTitle("Hits");
        FullDetCTDHistograms[s][DetID] = CTDHist;
      }

      if (HVHist == nullptr) {
        char name[64]; sprintf(name,"HV energy: DetID %d %s Illumination",DetID,IllumSide[s].Data());
        HVHist = new TH1D(name, name, (m_MaxEnergy - m_MinEnergy)*2, m_MinEnergy, m_MaxEnergy);
        HVHist->SetXTitle("HV Energy (keV)");
        HVHist->SetYTitle("Hits");
        FullDetHVEnergyHistograms[s][DetID] = HVHist;
      }

      if (LVHist == nullptr) {
        char name[64]; sprintf(name,"LV energy: DetID %d %s Illumination",DetID,IllumSide[s].Data());
        LVHist = new TH1D(name, name, (m_MaxEnergy - m_MinEnergy)*2, m_MinEnergy, m_MaxEnergy);
        LVHist->SetXTitle("LV Energy (keV)");
        LVHist->SetYTitle("Hits");
        FullDetLVEnergyHistograms[s][DetID] = LVHist;
      }

      CTDHist->Add(CTDHist, H.second);
      HVHist->Add(HVHist, HVEnergyHistograms[s][ID]);
      LVHist->Add(LVHist, LVEnergyHistograms[s][ID]);

      if (m_PixelCorrect==true) {
        if ((H.second->Integral() > g_MinCounts) && (HVEnergyHistograms[s][ID]->Integral() > g_MinCounts) && (LVEnergyHistograms[s][ID]->Integral() > g_MinCounts)) {

          TFitResultPtr CTDFit = H.second->Fit("gaus", "SQ");
          TFitResultPtr HVFit = HVEnergyHistograms[s][ID]->Fit("gaus", "SQ");
          TFitResultPtr LVFit = LVEnergyHistograms[s][ID]->Fit("gaus", "SQ");

          if ((!(CTDFit->IsEmpty())) && (!(HVFit->IsEmpty())) && (!(LVFit->IsEmpty()))) {
            Endpoints[ID][s].push_back(CTDFit->Parameter(1));
            Endpoints[ID][s].push_back(HVFit->Parameter(1));
            Endpoints[ID][s].push_back(LVFit->Parameter(1));
          } else {
            cout<<"Fits failed for Pixel "<<ID<<endl;
          }
        } else {
          cout<<"Fewer than "<<g_MinCounts<<" counts in Pixel ID "<<ID<<endl;
        }
      }
    }

    for (auto H: FullDetCTDHistograms[s]) {

      int DetID = H.first;

      if ((H.second->Integral() > g_MinCounts) && (FullDetHVEnergyHistograms[s][DetID]->Integral() > g_MinCounts) && (FullDetLVEnergyHistograms[s][DetID]->Integral() > g_MinCounts)) {

        TFitResultPtr CTDFit = H.second->Fit("gaus", "SQ");
        TFitResultPtr HVFit = FullDetHVEnergyHistograms[s][DetID]->Fit("gaus", "SQ");
        TFitResultPtr LVFit = FullDetLVEnergyHistograms[s][DetID]->Fit("gaus", "SQ");

        if ((!(CTDFit->IsEmpty())) && (!(HVFit->IsEmpty())) && (!(LVFit->IsEmpty()))) {
          FullDetEndpoints[DetID][s].push_back(CTDFit->Parameter(1));
          FullDetEndpoints[DetID][s].push_back(HVFit->Parameter(1));
          FullDetEndpoints[DetID][s].push_back(LVFit->Parameter(1));

          cout<<"Results of Det "<<DetID<<" CTD fit:"<<endl;
          CTDFit->Print();
          cout<<"Results of Det "<<DetID<<" HV energy fit:"<<endl;
          HVFit->Print();
          cout<<"Results of Det "<<DetID<<" LV energy fit:"<<endl;
          LVFit->Print();

          TFile CTDFile(m_OutFile+MString("_Det")+DetID+MString("_CTDHist_") + IllumSide[s]+ MString("Illum.root"),"recreate");

          TCanvas* CTDCanvas = new TCanvas();
          CTDCanvas->SetLogz();
          CTDCanvas->cd();
          H.second->Draw("colz");

          H.second->Write();
          CTDFile.Close();

          TFile HVHistFile(m_OutFile+MString("_Det")+DetID+MString("_HVEnergyHist_") + IllumSide[s]+ MString("Illum.root"),"recreate");

          TCanvas* HVHistCanvas = new TCanvas();
          HVHistCanvas->SetLogz();
          HVHistCanvas->cd();
          FullDetHVEnergyHistograms[s][DetID]->Draw("colz");

          FullDetHVEnergyHistograms[s][DetID]->Write();
          HVHistFile.Close();

          TFile LVHistFile(m_OutFile+MString("_Det")+DetID+MString("_LVEnergyHist_") + IllumSide[s]+ MString("Illum.root"),"recreate");

          TCanvas* LVHistCanvas = new TCanvas();
          LVHistCanvas->SetLogz();
          LVHistCanvas->cd();
          FullDetLVEnergyHistograms[s][DetID]->Draw("colz");

          FullDetLVEnergyHistograms[s][DetID]->Write();
          LVHistFile.Close();

        } else {
          cout<<"Fits failed for Det "<<DetID<<endl;
        }
      } else {
        cout<<"Fewer than "<<g_MinCounts<<" counts in Det ID "<<ID<<endl;
      }
    }
  }

  //setup output file
  ofstream OutputCalFile;
  OutputCalFile.open(m_OutFile+MString("_parameters.txt"));
  OutputCalFile<<"Det ID"<<"HV Strip ID"<<"LV Strip ID"<<'\t'<<"HV Slope"<<'\t'<<"LV Slope"<<endl<<endl;

  map<int, TH2D*> DeltaHV;
  map<int, TH2D*> DeltaLV;

  for (auto E: Endpoints) {
    int ID = E.first;
    if ((E.second[0].size() > 0) && (E.second[1].size() > 0)) {
      double HVSlope = (E.second[0][1] - E.second[1][1])/(E.second[0][0] - E.second[1][0]);
      double LVSlope = (E.second[0][2] - E.second[1][2])/(E.second[0][0] - E.second[1][0]);
      if (m_PixelCorrect==true) {
        int DetID = (ID-(ID%10000))/10000;
        int LVSHID = (ID-(DetID*10000) - (ID%100))/100;
        int HVSHID = (ID-(DetID*10000) - (LVSHID *100));
        OutputCalFile<<to_string(DetID)<<'\t'<<to_string(HVSHID)<<'\t'<<to_string(LVSHID)<<'\t'<<to_string(HVSlope)<<'\t'<<to_string(LVSlope)<<endl<<endl;

        TH2D* TempHV = DeltaHV[DetID];
        TH2D* TempLV = DeltaLV[DetID];
        if (TempHV == nullptr) {

          char HVname[64]; sprintf(HVname,"Delta HV Map: Det %d",DetID);
          TempHV = new TH2D(HVname, HVname, 64, 0, 64, 64, 0, 64);
          TempHV->SetXTitle("HV Strip");
          TempHV->SetYTitle("LV Strip");
          TempHV->SetZTitle("Delta HV Energy");
          
          char LVname[64]; sprintf(LVname,"Delta LV Map: Det %d",DetID);
          TempLV = new TH2D(LVname, LVname, 64, 0, 64, 64, 0, 64);
          TempLV->SetXTitle("HV Strip");
          TempLV->SetYTitle("LV Strip");
          TempLV->SetZTitle("Delta LV Energy");

          DeltaHV[DetID] = TempHV;
          DeltaLV[DetID] = TempLV;
        }

        DeltaHV[DetID]->SetBinContent(HVSHID, LVSHID, E.second[0][1] - E.second[1][1]);
        DeltaLV[DetID]->SetBinContent(HVSHID, LVSHID, E.second[0][2] - E.second[1][2]);

      } else {
        for (int hv=0; hv<64; ++hv) {
          for (int lv=0; lv<64; ++lv) {
            OutputCalFile<<to_string(ID)<<'\t'<<to_string(hv)<<'\t'<<to_string(lv)<<'\t'<<to_string(HVSlope)<<'\t'<<to_string(LVSlope)<<endl<<endl;
          }
        }
      }
    }
  }

  OutputCalFile.close();

  if (m_PixelCorrect==true) {
    for (auto H: DeltaHV) {
      
      int DetID = H.first;

      TFile HVFile(m_OutFile+MString("_Det")+DetID+MString("_DeltaHVMap.root"),"recreate");
      TCanvas* HVCanvas = new TCanvas();
      HVCanvas->cd();
      H.second->Draw("colz");
      H.second->Write();
      HVFile.Close();

      TFile LVFile(m_OutFile+MString("_Det")+DetID+MString("_DeltaLVMap.root"),"recreate");
      TCanvas* LVCanvas = new TCanvas();
      LVCanvas->cd();
      DeltaLV[DetID]->Draw("colz");
      DeltaLV[DetID]->Write();
      LVFile.Close();
    }
  }

  watch.Stop();
  cout<<"total time (s): "<<watch.CpuTime()<<endl;
 
  return true;
}


////////////////////////////////////////////////////////////////////////////////


TrappingCorrectionAm241* g_Prg = 0;
int g_NInterruptCatches = 1;

MStripHit* TrappingCorrectionAm241::GetDominantStrip(vector<MStripHit*>& Strips, double& EnergyFraction)
{
  double MaxEnergy = -numeric_limits<double>::max(); // AZ: When both energies are zero (which shouldn't happen) we still pick one
  double TotalEnergy = 0.0;
  MStripHit* MaxStrip = nullptr;

  // Iterate through strip hits and get the strip with highest energy
  for (const auto SH : Strips) {
    double Energy = SH->GetEnergy();
    TotalEnergy += Energy;
    if (Energy > MaxEnergy) {
      MaxStrip = SH;
      MaxEnergy = Energy;
    }
  }
  if (TotalEnergy == 0) {
    EnergyFraction = 0;
  } else {
    EnergyFraction = MaxEnergy/TotalEnergy;
  }
  return MaxStrip;
}


////////////////////////////////////////////////////////////////////////////////


//! Called when an interrupt signal is flagged
//! All catched signals lead to a well defined exit of the program
void CatchSignal(int a)
{
  if (g_Prg != 0 && g_NInterruptCatches-- > 0) {
    cout<<"Catched signal Ctrl-C (ID="<<a<<"):"<<endl;
    g_Prg->Interrupt();
  } else {
    abort();
  }
}


////////////////////////////////////////////////////////////////////////////////


//! Main program
int main(int argc, char** argv)
{
  // Catch a user interupt for graceful shutdown
  signal(SIGINT, CatchSignal);

  // Initialize global MEGALIB variables, especially mgui, etc.
  MGlobal::Initialize("Standalone", "a standalone example program");

  TApplication TrappingCorrectionApp("TrappingCorrectionApp", 0, 0);

  g_Prg = new TrappingCorrectionAm241();

  if (g_Prg->ParseCommandLine(argc, argv) == false) {
    cerr<<"Error during parsing of command line!"<<endl;
    return -1;
  } 
  if (g_Prg->Analyze() == false) {
    cerr<<"Error during analysis!"<<endl;
    return -2;
  } 

  TrappingCorrectionApp.Run();

  cout<<"Program exited normally!"<<endl;

  return 0;
}


////////////////////////////////////////////////////////////////////////////////
