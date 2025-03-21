/*
 * MModuleTACcut.h
 *
 * Copyright (C) by Andreas Zoglauer.
 * All rights reserved.
 *
 * Please see the source-file for the copyright-notice.
 *
 */


#ifndef __MModuleTACcut__
#define __MModuleTACcut__


////////////////////////////////////////////////////////////////////////////////


// Standard libs:
#include <algorithm>

// ROOT libs:
#include "TGClient.h"
#include "TH1.h"

// MEGAlib libs:
#include "MGlobal.h"
#include "MModule.h"
#include "MGUIExpoTACcut.h"


// Forward declarations:


////////////////////////////////////////////////////////////////////////////////


class MModuleTACcut : public MModule
{
  // public interface:
 public:
  //! Default constructor
  MModuleTACcut();
  //! Default destructor
  virtual ~MModuleTACcut();
  
  //! Create a new object of this class 
  virtual MModuleTACcut* Clone() { return new MModuleTACcut(); }

  //! Initialize the module
  virtual bool Initialize();

  //! Create expos 
  virtual void CreateExpos();

  //! Finalize the module
  virtual void Finalize();

  //! Main data analysis routine, which updates the event to a new level 
  virtual bool AnalyzeEvent(MReadOutAssembly* Event);

  //! Show the options GUI
  virtual void ShowOptionsGUI();


  //! Read the configuration data from an XML node
  virtual bool ReadXmlConfiguration(MXmlNode* Node);
  //! Create an XML node tree from the configuration
  virtual MXmlNode* CreateXmlConfiguration();

  ///////////// Creating functions that will update and get the min/max TAC values //////////////////////////

 //! Set the minimum TAC value!
  void SetMinimumTAC(unsigned int MinimumTAC) { m_MinimumTAC = MinimumTAC; }
  //! Get the minimum TAC value!
  unsigned int GetMinimumTAC() const { return m_MinimumTAC; }

  //! Set the maximum TAC value!
  void SetMaximumTAC(unsigned int MaximumTAC) { m_MaximumTAC = MaximumTAC; }
  //! Get the maximum TAC value!
  unsigned int GetMaximumTAC() const { return m_MaximumTAC; }

  //! Set filename for TAC Calibration
  void SetTACCalFileName( const MString& FileName) {m_TACCalFile = FileName;}
  //! Get filename for TAC Calibration
  MString GetTACCalFileName() const {return m_TACCalFile;}

  //! Set the shaping offset
  void SetShapingOffset(double ShapingOffset) { m_ShapingOffset = ShapingOffset; }
  //! Get the shaping offset
  unsigned int GetShapingOffset() const { return m_ShapingOffset; }

  //! Set the disable time
  void SetDisableTime(double DisableTime) { m_DisableTime = DisableTime; }
  //! Get the disable time
  unsigned int GetDisableTime() const { return m_DisableTime; }

  //! Set the shaping flag_to_en_delay
  void SetFlagToEnDelay(double FlagToEnDelay) { m_FlagToEnDelay = FlagToEnDelay; }
  //! Get the shaping flag_to_en_delay
  unsigned int GetFlagToEnDelay() const { return m_FlagToEnDelay; }

  //! Set the shaping coincidence window
  void SetCoincidenceWindow(double CoincidenceWindow) { m_CoincidenceWindow = CoincidenceWindow; }
  //! Get the shaping coincidence window
  unsigned int GetCoincidenceWindow() const { return m_CoincidenceWindow; }

  //! Load the TAC calibration file
  bool LoadTACCalFile(MString FName);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 


  // protected methods:
 protected:

  
  // private methods:
 private:



  // protected members:
 protected:

  // private members:
 private:

// declare min and max TAC variables here
unsigned int m_MinimumTAC, m_MaximumTAC;
double m_ShapingOffset, m_DisableTime, m_FlagToEnDelay, m_CoincidenceWindow;
MString m_TACCalFile;
unordered_map<int, unordered_map<int, vector<double>>> m_HVTACCal;
unordered_map<int, unordered_map<int, vector<double>>> m_LVTACCal;

vector<unsigned int> m_DetectorIDs;

MGUIExpoTACcut* m_ExpoTACcut;


#ifdef ___CLING___
 public:
  ClassDef(MModuleTACcut, 0) // no description
#endif

};

#endif


////////////////////////////////////////////////////////////////////////////////
