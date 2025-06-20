/*
 * MModuleTACcut.cxx
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


////////////////////////////////////////////////////////////////////////////////
//
// MModuleTACcut
//
////////////////////////////////////////////////////////////////////////////////


// Include the header:
#include "MModuleTACcut.h"

// Standard libs:

// ROOT libs:
#include "TGClient.h"

// MEGAlib libs:
#include "MModule.h"
#include "MGUIOptionsTACcut.h"
#include "MGUIExpoTACcut.h"


////////////////////////////////////////////////////////////////////////////////


#ifdef ___CLING___
ClassImp(MModuleTACcut)
#endif


////////////////////////////////////////////////////////////////////////////////


MModuleTACcut::MModuleTACcut() : MModule()
{
  // Construct an instance of MModuleTACcut

  // Set all module relevant information

  // Set the module name --- has to be unique
  m_Name = "TAC Calibration";

  // Set the XML tag --- has to be unique --- no spaces allowed
  m_XmlTag = "XmlTagTACcut";

  // Set all modules, which have to be done before this module
  AddPreceedingModuleType(MAssembly::c_EventLoader);
  //AddPreceedingModuleType(MAssembly::c_EnergyCalibration);

  // Set all types this modules handles
  AddModuleType(MAssembly::c_TACcut);

  // Set all modules, which can follow this module
  AddSucceedingModuleType(MAssembly::c_StripPairing);
  AddSucceedingModuleType(MAssembly::c_DepthCorrection);

  // Set if this module has an options GUI
  // Overwrite ShowOptionsGUI() with the call to the GUI!
  m_HasOptionsGUI = true;
  // If true, you have to derive a class from MGUIOptions (use MGUIOptionsTACcut)
  // and implement all your GUI options

  // Can the program be run multi-threaded
  m_AllowMultiThreading = true;

  // Can we use multiple instances of this class
  m_AllowMultipleInstances = false;

  m_DisableTime = 1396;
  m_FlagToEnDelay = 104;
}


////////////////////////////////////////////////////////////////////////////////


MModuleTACcut::~MModuleTACcut()
{
  // Delete this instance of MModuleTACcut
}


////////////////////////////////////////////////////////////////////////////////


bool MModuleTACcut::Initialize()
{
  // Initialize the module 

  if (LoadTACCalFile(m_TACCalFile) == false) {
    cout<<m_XmlTag<<": Error: TAC Calibration file could not be loaded."<<endl;
    return false;
  }

  if (LoadTACCutFile(m_TACCutFile) == false) {
    cout<<m_XmlTag<<": Error: TAC Calibration file could not be loaded."<<endl;
    return false;
  }

  // Some sanity checks:
  if (m_LVTACCal.size() == 0) {
    cout<<m_XmlTag<<": Error: The low voltage TAC calibration data set is empty"<<endl;
    return false;
  }
  if (m_HVTACCal.size() == 0) {
    cout<<m_XmlTag<<": Error: The high voltage TAC calibration data set is empty"<<endl;
    return false;
  }
  if (m_LVTACCut.size() == 0) {
    cout<<m_XmlTag<<": Error: The low voltage TAC cut data set is empty"<<endl;
    return false;
  }
  if (m_HVTACCal.size() == 0) {
    cout<<m_XmlTag<<": Error: The high voltage TAC cut data set is empty"<<endl;
    return false;
  }

  return MModule::Initialize();
}


////////////////////////////////////////////////////////////////////////////////


void MModuleTACcut::CreateExpos()
{
  // Create all expos

  if (HasExpos() == true) return;

  // Set the histogram display
  m_ExpoTACcut = new MGUIExpoTACcut(this);
  m_ExpoTACcut->SetTACHistogramArrangement(m_DetectorIDs);
  for (unsigned int i = 0; i < m_DetectorIDs.size(); ++i) {
    unsigned int DetID = m_DetectorIDs[i];
    m_ExpoTACcut->SetTACHistogramParameters(DetID, 120, 0, 7000);
  }
  m_Expos.push_back(m_ExpoTACcut);
}

////////////////////////////////////////////////////////////////////////////////


bool MModuleTACcut::AnalyzeEvent(MReadOutAssembly* Event) 
{
  // Start with sanity checks:
  for (unsigned int i = 0; i < Event->GetNStripHits(); ++i) {
    MStripHit* SH = Event->GetStripHit(i);

    int DetID = SH->GetDetectorID();
    int StripID = SH->GetStripID();

    if (DetID >= m_LVTACCal.size()) {
      cout<<m_XmlTag<<": Error: DetID "<<DetID<<" is not in LVTACCal (max det ID: "<<m_LVTACCal.size()-1<<") - skipping event"<<endl;
      return false;
    }
    if (StripID >= m_LVTACCal[DetID].size()) {
      cout<<m_XmlTag<<": Error: StripID "<<StripID<<" is not in LVTACCal (max strip ID: "<<m_LVTACCal[DetID].size()-1<<") - skipping event"<<endl;
      return false;
    }
    if (DetID >= m_HVTACCal.size()) {
      cout<<m_XmlTag<<": Error: DetID "<<DetID<<" is not in HVTACCal (max det ID: "<<m_HVTACCal.size()-1<<") - skipping event"<<endl;
      return false;
    }
    if (StripID >= m_HVTACCal[DetID].size()) {
      cout<<m_XmlTag<<": Error: StripID "<<StripID<<" is not in HVTACCal (max strip ID: "<<m_HVTACCal[DetID].size()-1<<") - skipping event"<<endl;
      return false;
    }
  }

  double MaxTAC = -numeric_limits<double>::max();
  for (unsigned int i = 0; i < Event->GetNStripHits(); ++i) {
    MStripHit* SH = Event->GetStripHit(i);
    double TAC_timing = SH->GetTAC();
    double ns_timing;
    int DetID = SH->GetDetectorID();
    int StripID = SH->GetStripID();
    if (SH->IsLowVoltageStrip() == true) {
      ns_timing = TAC_timing*m_LVTACCal[DetID][StripID][0] + m_LVTACCal[DetID][StripID][1];
    } else {
      ns_timing = TAC_timing*m_HVTACCal[DetID][StripID][0] + m_HVTACCal[DetID][StripID][1];
    }
    SH->SetTiming(ns_timing);
    if (ns_timing > MaxTAC) {
      MaxTAC = ns_timing;
    }
  }

  for (unsigned int i = 0; i < Event->GetNStripHits();) {
    MStripHit* SH = Event->GetStripHit(i);
    double SHTiming = SH->GetTiming();
    int DetID = SH->GetDetectorID();
    int StripID = SH->GetStripID();
    double ShapingOffset;
    double CoincidenceWindow;
    if (SH->IsLowVoltageStrip() == true) {
      ShapingOffset = m_LVTACCut[DetID][StripID][0];
      CoincidenceWindow = m_LVTACCut[DetID][StripID][1];
    } else {
      ShapingOffset = m_HVTACCut[DetID][StripID][0];
      CoincidenceWindow = m_HVTACCut[DetID][StripID][1];
    }
    double TotalOffset = ShapingOffset + m_DisableTime + m_FlagToEnDelay;
    if ((SHTiming < TotalOffset + CoincidenceWindow) && (SHTiming > TotalOffset) && (SHTiming > MaxTAC - CoincidenceWindow)){
      if (HasExpos()==true) {
        m_ExpoTACcut->AddTAC(DetID, SHTiming);
      }
      ++i;
    } else {
      Event->RemoveStripHit(i);
      delete SH;
    }
  }

  Event->SetAnalysisProgress(MAssembly::c_TACcut);

  return true;
}


////////////////////////////////////////////////////////////////////////////////


void MModuleTACcut::Finalize()
{
  MModule::Finalize();
}


////////////////////////////////////////////////////////////////////////////////


void MModuleTACcut::ShowOptionsGUI()
{
  //! Show the options GUI --- has to be overwritten!

  MGUIOptionsTACcut* Options = new MGUIOptionsTACcut(this);
  Options->Create();
  gClient->WaitForUnmap(Options);
}


////////////////////////////////////////////////////////////////////////////////


bool MModuleTACcut::ReadXmlConfiguration(MXmlNode* Node)
{
  //! Read the configuration data from an XML node

  MXmlNode* TACCalFileNameNode = Node->GetNode("TACCalFileName");
  if (TACCalFileNameNode != nullptr) {
    SetTACCalFileName(TACCalFileNameNode->GetValue());
  }

  MXmlNode* TACCutFileNameNode = Node->GetNode("TACCutFileName");
  if (TACCutFileNameNode != nullptr) {
    SetTACCutFileName(TACCutFileNameNode->GetValue());
  }

  MXmlNode* DisableTimeNode = Node->GetNode("DisableTime");
  if (DisableTimeNode != nullptr) {
    m_DisableTime = DisableTimeNode->GetValueAsDouble();
  }

  MXmlNode* FlagToEnDelayNode = Node->GetNode("FlagToEnDelay");
  if (FlagToEnDelayNode != nullptr) {
    m_FlagToEnDelay = FlagToEnDelayNode->GetValueAsDouble();
  }

  return true;
}


////////////////////////////////////////////////////////////////////////////////


MXmlNode* MModuleTACcut::CreateXmlConfiguration() 
{
  //! Create an XML node tree from the configuration

  MXmlNode* Node = new MXmlNode(0, m_XmlTag);
  
  new MXmlNode(Node, "TACCalFileName", m_TACCalFile);
  new MXmlNode(Node, "TACCutFileName", m_TACCutFile);
  new MXmlNode(Node, "DisableTime", m_DisableTime);
  new MXmlNode(Node, "FlagToEnDelay", m_FlagToEnDelay);

  return Node;
}


////////////////////////////////////////////////////////////////////////////////


bool MModuleTACcut::LoadTACCalFile(MString FName)
{
  // Read in the TAC Calibration file, which should contain for each strip:
  //  DetID, Side (h or l for high or low voltage), TAC cal, TAC cal error, TAC cal offset, TAC offset error
  // OR:
  // ReadOutID, Detector, Side, Strip, TAC cal, TAC cal error, TAC offset, TAC offset error
  MFile F;
  if (F.Open(FName) == false) {
    cout << "MModuleTACcut: failed to open TAC Calibration file." << endl;
    return false;
  } else {
    MString Line;
    while (F.ReadLine(Line)) {
      if (!Line.BeginsWith("#")) {
        std::vector<MString> Tokens = Line.Tokenize(",");
        if ((Tokens.size() == 7) || (Tokens.size() == 8)) {
          int IndexOffset = Tokens.size() % 7;
          int DetID = Tokens[0+IndexOffset].ToInt();
          int StripID = Tokens[2+IndexOffset].ToInt();
          double TACCal = Tokens[3+IndexOffset].ToDouble();
          double TACCalError = Tokens[4+IndexOffset].ToDouble();
          double Offset = Tokens[5+IndexOffset].ToDouble();
          double OffsetError = Tokens[6+IndexOffset].ToDouble();
          vector<double> CalValues;
          CalValues.push_back(TACCal); CalValues.push_back(Offset); CalValues.push_back(TACCalError); CalValues.push_back(OffsetError);
          
          if (find(m_DetectorIDs.begin(), m_DetectorIDs.end(), DetID) == m_DetectorIDs.end()) {
            unordered_map<int, vector<double>> TempMapHV;
            m_HVTACCal[DetID] = TempMapHV;
            unordered_map<int, vector<double>> TempMapLV;
            m_LVTACCal[DetID] = TempMapLV;
            m_DetectorIDs.push_back(DetID);
          }
          
          if ((Tokens[1+IndexOffset] == "l") || (Tokens[1+IndexOffset] == "0")) {
            m_LVTACCal[DetID][StripID] = CalValues;
          } else if ((Tokens[1+IndexOffset] == "h") || (Tokens[1+IndexOffset] == "1")) {
            m_HVTACCal[DetID][StripID] = CalValues;
          }
        }
      }
    }
    F.Close();
    sort(m_DetectorIDs.begin(), m_DetectorIDs.end());
  }

  return true;
}


////////////////////////////////////////////////////////////////////////////////


bool MModuleTACcut::LoadTACCutFile(MString FName)
{
  // Read in the TAC Cut file, which should contain for each strip:
  //  DetID, h or l for high or low voltage, shaping offset, coincidence window
  MFile F;
  if (F.Open(FName) == false) {
    cout << "MModuleTACcut: failed to open TAC Cut file." << endl;
    return false;
  } else {
    MString Line;
    while (F.ReadLine(Line)) {
      if (!Line.BeginsWith("#")) {
        std::vector<MString> Tokens = Line.Tokenize(",");
        if (Tokens.size() == 5) {
          int DetID = Tokens[0].ToInt();
          int StripID = Tokens[2].ToInt();
          double ShapingOffset = Tokens[3].ToDouble();
          double CoincidenceWindow = Tokens[4].ToDouble();
          vector<double> CutParams;
          CutParams.push_back(ShapingOffset); CutParams.push_back(CoincidenceWindow);
          
          if (find(m_DetectorIDs.begin(), m_DetectorIDs.end(), DetID) == m_DetectorIDs.end()) {
            unordered_map<int, vector<double>> TempMapHV;
            m_HVTACCut[DetID] = TempMapHV;
            unordered_map<int, vector<double>> TempMapLV;
            m_LVTACCut[DetID] = TempMapLV;
            m_DetectorIDs.push_back(DetID);
          }
          
          if (Tokens[1] == "l") {
            m_LVTACCut[DetID][StripID] = CutParams;
          } else if (Tokens[1] == "h") {
            m_HVTACCut[DetID][StripID] = CutParams;
          }
        }
      }
    }
    F.Close();
    sort(m_DetectorIDs.begin(), m_DetectorIDs.end());
  }

  return true;

}

// MModuleTACcut.cxx: the end...
////////////////////////////////////////////////////////////////////////////////


