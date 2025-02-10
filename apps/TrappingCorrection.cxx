/* 
 * TrappingCorrection.cxx
 *
 *
 * Copyright (C) by Andreas Zoglauer.
 * All rights reserved.
 *
 *
 * This code implementation is the intellectual property of
 * Andreas Zoglauer.
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
#include <TStopwatch.h>
#include <TProfile.h>
#include <Minuit2/FCNBase.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/MnApplication.h>
#include <Minuit2/MnStrategy.h>
using namespace ROOT::Minuit2;

// MEGAlib
#include "MGlobal.h"
#include "MFile.h"
#include "MReadOutElementDoubleStrip.h"
#include "MFileReadOuts.h"
#include "MReadOutAssembly.h"
#include "MStripHit.h"
#include "MReadOutSequence.h"
#include "MSupervisor.h"
#include "MModuleLoaderMeasurementsROA.h"
#include "MModuleEnergyCalibrationUniversal.h"
#include "MModuleEventFilter.h"
#include "MModuleStripPairingGreedy.h"
#include "MModuleStripPairingChiSquare.h"
#include "MModuleTACcut.h"
#include "MAssembly.h"


////////////////////////////////////////////////////////////////////////////////


class SymmetryFCN : public FCNBase
{
public:
	//! Operator which returns the symmetry of m_CTDHistogram given the parameters passed
	double operator()(vector<double> const &v) const override;
	double Up() const override { return 1; }

	void AddCTD(double CTD){ m_CTDs.push_back(CTD); }
	void AddHVEnergy(double HVEnergy){ m_HVEnergies.push_back(HVEnergy); }
	void AddLVEnergy(double LVEnergy){ m_LVEnergies.push_back(LVEnergy); }

private:

	//! The measured CTD and HV/LV energies
	vector<double> m_CTDs;
	vector<double> m_HVEnergies;
	vector<double> m_LVEnergies;

};


////////////////////////////////////////////////////////////////////////////////


//! A standalone program based on MEGAlib and ROOT
class TrappingCorrection
{
public:
  //! Default constructor
  TrappingCorrection();
  //! Default destructor
  ~TrappingCorrection();
  
  //! Parse the command line
  bool ParseCommandLine(int argc, char** argv);
  //! Analyze what eveer needs to be analyzed...
  bool Analyze();
	//!load cross talk correction
	vector<vector<vector<float> > > LoadCrossTalk();
  //! Interrupt the analysis
  void Interrupt() { m_Interrupt = true; }

	void dummy_func() { return; }

	MStripHit* GetDominantStrip(vector<MStripHit*>& Strips, double& EnergyFraction);

private:
  //! True, if the analysis needs to be interrupted
  bool m_Interrupt;
  //! The input file name
  MString m_FileName;
  MString m_EcalFile;
  MString m_TACFile;
	//! output file names
	MString m_OutFile;
	//! energy E0
	// float m_E0;
	//! option to correct charge loss or not
	// bool m_CorrectCL;
	//! option to do a pixel-by-pixel calibration (instead of detector-by-detector)
	bool m_PixelCorrect;
	bool m_GreedyPairing;
	bool m_CardCageOverride;

};

////////////////////////////////////////////////////////////////////////////////


double SymmetryFCN::operator()(vector<double> const &v) const
{
	double HVSlope = v[0];
	double HVIntercept = v[1];
	double LVSlope = v[2];
	double LVIntercept = v[3];

	char name[64]; sprintf(name,"name");
	TH2D* CorrectedHistogram = new TH2D(name, name, 51, -250, 250, 51, 0.95, 1.05);

	for (unsigned int i = 0; i < m_CTDs.size(); ++i) {
		// Correct the HV and LV energies by dividing by the CCE. DeltaCCE is defined as a linear function with units percentage energy lost.
		double CorrectedHVEnergy = m_HVEnergies[i]/(1 - (HVSlope*m_CTDs[i] + HVIntercept)/100);
		double CorrectedLVEnergy = m_LVEnergies[i]/(1 - (LVSlope*m_CTDs[i] + LVIntercept)/100);
		CorrectedHistogram->Fill(m_CTDs[i], CorrectedHVEnergy/CorrectedLVEnergy);
		// m_CorrectedHistogramReflected->Fill(-m_CTDs[i], CorrectedHVEnergy/CorrectedLVEnergy)
	}

	vector<vector<double>> BinValues;
	vector<vector<double>> ReflectedBinValues;
	double Symmetry = 0;

	for (unsigned int y = 0; y < CorrectedHistogram->GetNbinsY(); ++y) {
		
		vector<double> XValues;

		for (unsigned int x = 0; x < CorrectedHistogram->GetNbinsX(); ++x) {
			XValues.push_back(CorrectedHistogram->GetBinContent(x,y));
		}

		// BinValues.push_back(XValues);
		vector<double> ReflectedXValues = XValues;
		reverse(ReflectedXValues.begin(), ReflectedXValues.end());

		for (unsigned int x = 0; x < XValues.size(); ++x) {
			Symmetry += XValues[x] * ReflectedXValues[x];
		}
		// ReflectedBinValues.push_back(ReflectedXValues);
	}

	Symmetry /= 2;

	return -2*log(Symmetry);

}


////////////////////////////////////////////////////////////////////////////////


//! Default constructor
TrappingCorrection::TrappingCorrection() : m_Interrupt(false)
{
  gStyle->SetPalette(1, 0);
}


////////////////////////////////////////////////////////////////////////////////


//! Default destructor
TrappingCorrection::~TrappingCorrection()
{
  // Intentionally left blank
}


////////////////////////////////////////////////////////////////////////////////


//! Parse the command line
bool TrappingCorrection::ParseCommandLine(int argc, char** argv)
{
  ostringstream Usage;
  Usage<<endl;
  Usage<<"  Usage: TrappingCorrection <options>"<<endl;
  Usage<<"    General options:"<<endl;
  Usage<<"         -i:   input file name (.roa)"<<endl;
  Usage<<"         -e:   energy calibration file (.ecal)"<<endl;
  Usage<<"         -t:   TAC calibration file"<<endl;
	Usage<<"         -p:   do pixel-by-pixel correction"<<endl;
	Usage<<"         -g:   greedy strip pairing (default is chi-square)"<<endl;
	Usage<<"         -c:   Card cage data (i.e. no TAC calibration required)"<<endl;
	Usage<<"         -o:   outfile"<<endl;
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

  // Now parse the command line options:
  for (int i = 1; i < argc; i++) {
    Option = argv[i];

    // First check if each option has sufficient arguments:
    // Single argument
    if (Option == "-i" || Option == "-e" || Option == "-t" || Option == "-o") {
      if (!((argc > i+1) && 
            (argv[i+1][0] != '-' || isalpha(argv[i+1][1]) == 0))){
        cout<<"Error: Option "<<argv[i][1]<<" needs a second argument!"<<endl;
        cout<<Usage.str()<<endl;
        return false;
      }
    } 
    // Multiple arguments template
    /*
    else if (Option == "-??") {
      if (!((argc > i+2) && 
            (argv[i+1][0] != '-' || isalpha(argv[i+1][1]) == 0) && 
            (argv[i+2][0] != '-' || isalpha(argv[i+2][1]) == 0))){
        cout<<"Error: Option "<<argv[i][1]<<" needs two arguments!"<<endl;
        cout<<Usage.str()<<endl;
        return false;
      }
    }
    */

    // Then fulfill the options:
    if (Option == "-i") {
      m_FileName = argv[++i];
      cout<<"Accepting file name: "<<m_FileName<<endl;
    } 

    if (Option == "-e") {
      m_EcalFile = argv[++i];
      cout<<"Accepting file name: "<<m_EcalFile<<endl;
    } 

    if (Option == "-t") {
      m_TACFile = argv[++i];
      cout<<"Accepting file name: "<<m_TACFile<<endl;
    } 

		if (Option == "-o"){
			m_OutFile = argv[++i];
			cout<<"Accepting file name: "<<m_OutFile<<endl;
		}

		if (Option == "-c"){
			m_CardCageOverride = true;
		}

		if (Option == "-p"){
			m_PixelCorrect = true;
		}

		if (Option == "-g"){
			m_GreedyPairing = true;
		}

	}

  return true;
}


////////////////////////////////////////////////////////////////////////////////


//! Do whatever analysis is necessary
bool TrappingCorrection::Analyze()
{

/*	TH2D* h2 = new TH2D("fracmap","fracmap",356*2,-356,356,35,356-30,356+5);
//	for (int i=-662; i<662; i++){
//		for (int j=662-30; j<662+5; j++){
	for (int i=0; i<356*2; i++){
		for (int j=0; j<35; j++){
			float c = ((float)i-356)/((float)j+356-30);
			cout<<c<<endl;
			h2->SetBinContent(i,j,c);
		}
	}
	TCanvas *ctemp = new TCanvas();
	h2->Draw("colz");
	ctemp->Print("frac_map.pdf");
*/	
	//time code just to see
	TStopwatch watch;
	watch.Start();

	if (m_Interrupt == true) return false;

  MSupervisor* S = MSupervisor::GetSupervisor();
  
  cout<<"Creating ROA loader"<<endl;
	MModuleLoaderMeasurementsROA* Loader = new MModuleLoaderMeasurementsROA();
  Loader->SetFileName(m_FileName);
  S->SetModule(Loader, 0);

	cout<<"Creating TAC calibrator"<<endl;
  MModuleTACcut* TACCalibrator = new MModuleTACcut();
  TACCalibrator->SetTACCalFileName(m_TACFile);
  S->SetModule(TACCalibrator, 1);
 
 	cout<<"Creating energy calibrator"<<endl;
  MModuleEnergyCalibrationUniversal* EnergyCalibrator = new MModuleEnergyCalibrationUniversal();
  EnergyCalibrator->SetFileName(m_EcalFile);
	EnergyCalibrator->EnablePreampTempCorrection(false);
  S->SetModule(EnergyCalibrator, 2);

  cout<<"Creating Event filter"<<endl;
  //! Only use events with 1 Strip Hit on each side to avoid strip pairing complications
  MModuleEventFilter* EventFilter = new MModuleEventFilter();
  EventFilter->SetMinimumLVStrips(1);
  EventFilter->SetMaximumLVStrips(1);
  EventFilter->SetMinimumHVStrips(1);
  EventFilter->SetMaximumHVStrips(1);
  S->SetModule(EventFilter, 3);
  
  cout<<"Creating strip pairing"<<endl;
  MModule* Pairing;
  if (m_GreedyPairing == true){
  	// Pairing = dynamic_cast<MModuleStripPairingChiSquare*>(Pairing);
  	Pairing = new MModuleStripPairingGreedy();
  }
  else {
  	// Pairing = dynamic_cast<MModuleStripPairingChiSquare*>(Pairing);
  	Pairing = new MModuleStripPairingChiSquare();
  }
  S->SetModule(Pairing, 4);

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

  map<int, TH2D*> Histograms;
  map<int, SymmetryFCN*> FCNs;

  bool IsFinished = false;
  MReadOutAssembly* Event = new MReadOutAssembly();

  cout<<"Analyzing..."<<endl;
  while (IsFinished == false && m_Interrupt == false) {
    Event->Clear();

		if (Loader->IsReady() ){
	    Loader->AnalyzeEvent(Event);
	    TACCalibrator->AnalyzeEvent(Event);
	    EnergyCalibrator->AnalyzeEvent(Event);
	    bool Unfiltered = EventFilter->AnalyzeEvent(Event);
	    Pairing->AnalyzeEvent(Event);

			if ((Event->HasAnalysisProgress(MAssembly::c_StripPairing) == true) && (Unfiltered==true)) {
				for (unsigned int h = 0; h < Event->GetNHits(); ++h) {
					double HVEnergy = 0.0;
					double LVEnergy = 0.0;
					vector<MStripHit*> HVStrips;
		      vector<MStripHit*> LVStrips;

		      MHit* H = Event->GetHit(h);
		      
		      int DetID = H->GetStripHit(0)->GetDetectorID();
		      TH2D* Hist = Histograms[DetID];
		      SymmetryFCN* FCN = FCNs[DetID]; 
					
					if (Hist == nullptr) {
						char name[64]; sprintf(name,"%d",DetID);
						Hist = new TH2D(name, name, 51, -250, 250, 51, 0.95, 1.05);
						Histograms[DetID] = Hist;
					}

					if (FCN == nullptr) {
						FCN = new SymmetryFCN();
						FCNs[DetID] = FCN;
					}

		      for (unsigned int sh = 0; sh < H->GetNStripHits(); ++sh) {
		        MStripHit* SH = H->GetStripHit(sh);

						if (SH->IsLowVoltageStrip()==true) {
							LVEnergy += SH->GetEnergy();
							LVStrips.push_back(SH);
						} else {
							HVEnergy += SH->GetEnergy();
							HVStrips.push_back(SH);
		      	}
		      }

		      double HVEnergyFraction = 0;
		      double LVEnergyFraction = 0;
		      MStripHit* HVSH = GetDominantStrip(HVStrips, HVEnergyFraction); 
		      MStripHit* LVSH = GetDominantStrip(LVStrips, LVEnergyFraction); 
					double EnergyFraction = HVEnergy/LVEnergy;
					double CTD = HVSH->GetTiming() - LVSH->GetTiming();
					Hist->Fill(CTD, EnergyFraction);
					FCN->AddCTD(CTD);
					FCN->AddHVEnergy(HVEnergy);
					FCN->AddLVEnergy(LVEnergy);
				}
			}
		}
    IsFinished = Loader->IsFinished();
	}

	//setup output file
	ofstream logFitStats;
	logFitStats.open("TrappingCorrection.log");
	logFitStats<<"Det"<<'\t'<<"HV Slope"<<'\t'<<"HV Intercept"<<'\t'<<"LV Slope"<<'\t'<<"LV Intercept"<<endl<<endl;

	for (auto H: Histograms) {
//		if (H.first == 0){
		TCanvas* C = new TCanvas();
		C->SetLogz();
		C->cd();
		H.second->Draw("colz");

		int det = H.first;
		TFile f(m_OutFile+MString("_Det")+det+MString("_Hist.root"),"new");
		H.second->Write();
		f.Close();

	}

	for (auto F: FCNs) {

		MnUserParameters* InitialState = new MnUserParameters();
		InitialState->Add("HVSlope", -1e-3, 1e-4, -10, 0);
		InitialState->Add("HVIntercept", 0.5, 0.1, -10, 10);
		InitialState->Add("LVSlope", 1e-3, 1e-4, 0, 10);
		InitialState->Add("LVIntercept", 0.5, 0.1, -10, 10);

		MnMigrad migrad(*F.second, *InitialState);
		 // Minimize
  	FunctionMinimum min = migrad();
	  // output
	  cout<<min<<endl;
	}

	watch.Stop();
	cout<<"total time (s): "<<watch.CpuTime()<<endl;
 
  return true;
}


////////////////////////////////////////////////////////////////////////////////


TrappingCorrection* g_Prg = 0;
int g_NInterruptCatches = 1;

MStripHit* TrappingCorrection::GetDominantStrip(vector<MStripHit*>& Strips, double& EnergyFraction)
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

  g_Prg = new TrappingCorrection();

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
