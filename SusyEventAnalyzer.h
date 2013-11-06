// -*- C++ -*-
//
// Package:    SusyNtuplizer
// Class:      SusyEventAnalyzer.h
//
/*

 Description: an analyzer for susy::Event

 Implementation:

*/
//
// Original Author:  Dongwook Jang
// $Id: SusyEventAnalyzer.h,v 1.6 2013/06/21 20:02:11 weinberg Exp $
//

#ifndef SusyEventAnalyzer_h
#define SusyEventAnalyzer_h

#include <TChain.h>
#include <TString.h>
#include <TPRegexp.h>
#include <TArrayI.h>

#include <map>
#include <set>
#include <fstream>

#include "SusyEvent.h"

class SusyEventAnalyzer {
public:
  SusyEventAnalyzer(TTree&);
  virtual ~SusyEventAnalyzer();

  /* Main analyzer function to be defined */
  virtual void Run();

  /* parameter configuration functions */
  void IncludeAJson(TString const&);
  void SetOutputDirectory(TString const& v) { outputDirectoryName = v; }
  void SetOutput(TString const& v) { outputName = v; }
  void SetLogFile(TString const& v) { logFileName = v;}
  void SetPrintInterval(int v) { printInterval = v; }
  void SetPrintLevel(int v) { printLevel = v; }
  void SetProcessNEvents(int v) { processNEvents = v; }
  void AddHltName(TString const& v) { hltNames.push_back(v + "_v*"); }
  void CopyEvents(bool v) { copyEvents = v; }

  void SetNScaledEvents(Float_t lumi, Float_t xSec) { lumiCalc = lumi; nScaledEvents = lumi * xSec; }
  void SetDatasetName  (TString const& v)           { datasetName = v; }
  void SetCutComplement(TString const& v)           { cutComplement = v; SetObjectCuts(cutComplement); }
  void SetJecUncert    (Float_t jec)                { jecSystematic = jec; }

  void    SetObjectCuts       (TString);
  Float_t SetSusyXSec         (Float_t);
  Float_t SetSignalEventWeight(TString, Float_t, Float_t);

  void createHistogram(const char*, const char*, const char*, const char*, Int_t, Double_t, Double_t);

protected:
  bool IsGoodLumi(UInt_t, UInt_t) const;
  bool PassTriggers() const;

  /* container of all event data */
  susy::Event event;
  /* input tree */
  TTree *fTree;
  /* directory of the output file */
  TString outputDirectoryName;
  /* suffix of the output file */
  TString outputName;
  /* log file name */
  TString logFileName;
  /* verbosity - 0 => no printout, 1 => print function control flow, 2 => print event processing flow, 3 => print event dump */
  int printLevel;
  /* print frequency */
  unsigned printInterval;
  /* maximum number of events */
  int processNEvents;
  /* HLT path names */
  std::vector<TString> hltNames;
  /* switch for saving skims */
  bool copyEvents;
  /* good lumi list */
  std::map<unsigned, std::set<unsigned> > goodLumiList;
  mutable std::pair<unsigned, unsigned> currentLumi;
  mutable bool currentLumiIsGood;

  std::map<TString, TH1F*> hName;  // Map for histograms
  Float_t lumiCalc;                // Luminosity used for scaling output
  Float_t nScaledEvents;           // Number of events to be used for scaling output histograms
  TString datasetName;
  TString cutComplement;
  Float_t jecSystematic;           // -1.0, 0.0, +1.0 for JEC down, nominal, JEC up

  // Cut variables
  Float_t muon1_ptCut,   muon2_ptCut,   electron1_etCut, electron2_etCut;
  Float_t photon1_etCut, photon2_etCut, jet_ptCut,       met_etCut;

  // Flat ntuple variables
  Int_t   runNo,      lumiNo;
  ULong_t eventNo;
  Float_t squark_m,   neutralino_m;
  Int_t   vertices_n, muons_n, electrons_n, photons_n, jets_n;
  Float_t evtWt;
  Float_t met_et,     met_phi, st;

  std::vector<float>* muon_e;
  std::vector<float>* muon_px;
  std::vector<float>* muon_py;
  std::vector<float>* muon_pz;
  std::vector<float>* muon_charge;
  std::vector<float>* electron_e;
  std::vector<float>* electron_px;
  std::vector<float>* electron_py;
  std::vector<float>* electron_pz;
  std::vector<float>* electron_charge;
  std::vector<float>* photon_e;
  std::vector<float>* photon_px;
  std::vector<float>* photon_py;
  std::vector<float>* photon_pz;
  std::vector<float>* jet_e;
  std::vector<float>* jet_px;
  std::vector<float>* jet_py;
  std::vector<float>* jet_pz;
  std::vector<bool>*  jet_bTag;
};

SusyEventAnalyzer::SusyEventAnalyzer(TTree& tree) :
  event(),
  fTree(&tree),
  outputDirectoryName("."),
  outputName("analysis"),
  logFileName(outputName + ".log"),
  printLevel(0),
  printInterval(1000),
  processNEvents(-1),
  hltNames(),
  copyEvents(false),
  goodLumiList(),
  currentLumi(0, 0),
  currentLumiIsGood(true),
  lumiCalc(1.0),
  nScaledEvents(-1.0),
  datasetName(""),
  cutComplement(""),
  jecSystematic(0.0)
{
  event.setInput(tree);

  muon1_ptCut     = 15.0;
  muon2_ptCut     = 15.0;
  electron1_etCut = 15.0;
  electron2_etCut = 15.0;
  photon1_etCut   = 15.0;
  photon2_etCut   = 15.0;
  jet_ptCut       = 30.0;
  met_etCut       = 15.0;
}

SusyEventAnalyzer::~SusyEventAnalyzer()
{
}

void
SusyEventAnalyzer::IncludeAJson(TString const& _fileName)
{
  if(_fileName == "") return;

  std::ifstream inputFile(_fileName);
  if(!inputFile.is_open()){
    std::cerr << "Cannot open JSON file " << _fileName << std::endl;
    return;
  }

  std::string line;
  TString jsonText;
  while(true){
    std::getline(inputFile, line);
    if(!inputFile.good()) break;
    jsonText += line;
  }
  inputFile.close();

  TPRegexp runBlockPat("\"([0-9]+)\":[ ]*\\[((?:\\[[0-9]+,[ ]*[0-9]+\\](?:,[ ]*|))+)\\]");
  TPRegexp lumiBlockPat("\\[([0-9]+),[ ]*([0-9]+)\\]");

  TArrayI positions(2);
  positions[1] = 0;
  while(runBlockPat.Match(jsonText, "g", positions[1], 10, &positions) == 3){
    TString runBlock(jsonText(positions[0], positions[1] - positions[0]));
    TString lumiPart(jsonText(positions[4], positions[5] - positions[4]));

    unsigned run(TString(jsonText(positions[2], positions[3] - positions[2])).Atoi());
    std::set<unsigned>& lumis(goodLumiList[run]);

    TArrayI lumiPos(2);
    lumiPos[1] = 0;
    while(lumiBlockPat.Match(lumiPart, "g", lumiPos[1], 10, &lumiPos) == 3){
      TString lumiBlock(lumiPart(lumiPos[0], lumiPos[1] - lumiPos[0]));
      int begin(TString(lumiPart(lumiPos[2], lumiPos[3] - lumiPos[2])).Atoi());
      int end(TString(lumiPart(lumiPos[4], lumiPos[5] - lumiPos[4])).Atoi());
      for(int lumi(begin); lumi <= end; ++lumi)
        lumis.insert(lumi);
    }
  }
}



void
SusyEventAnalyzer::SetObjectCuts(TString cuts)
{
  if (cuts == "singleMu")
    muon1_ptCut = 30.0;
  else if (cuts == "singleElectron")
    electron1_etCut = 30.0;
  else if (cuts == "doublePhoton") {
    photon1_etCut = 40.0;
    photon2_etCut = 25.0;
  }
  else if (cuts == "photonHad")
    photon1_etCut = 75.0;
}



bool
SusyEventAnalyzer::IsGoodLumi(UInt_t run, UInt_t lumi) const
{
  if(goodLumiList.size() == 0) return true;
  if(run == currentLumi.first && lumi == currentLumi.second) return currentLumiIsGood;
  currentLumi.first = run;
  currentLumi.second = lumi;
  currentLumiIsGood = false;

  std::map<unsigned, std::set<unsigned> >::const_iterator rItr(goodLumiList.find(run));
  if(rItr != goodLumiList.end()){
    std::set<unsigned>::const_iterator lItr(rItr->second.find(lumi));
    if(lItr != rItr->second.end()) currentLumiIsGood = true;
  }

  return currentLumiIsGood;
}

bool
SusyEventAnalyzer::PassTriggers() const
{
  unsigned nT(hltNames.size());
  if(nT == 0) return true;

  for(unsigned iT(0); iT != nT; ++iT)
    if(event.hltMap.pass(hltNames[iT])) return true;

  return false;
}

#endif

