// Original Author:  Dongwook Jang
// $Id: ana.C,v 1.1 2013/06/21 19:39:12 weinberg Exp $

void ana(TString outputName="stealthAnalysis"){

  // REMOVE THE LINE BELOW IF NOT RUNNING IN CMSSW ENVIRONMENT
  gSystem->Load("libCondFormatsJetMETObjects.so");

  gSystem->Load("libSusyEvent.so");

  gSystem->AddIncludePath("-I" + TString(gSystem->Getenv("CMSSW_RELEASE_BASE")) + "/src");

  // Analysis macro
  gROOT->LoadMacro("SusyEventAnalyzer.cc+");

  // chain of inputs
  TChain chain("susyTree");
//   chain.Add("susyEvents.root");
  chain.Add("dcap:///pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/stealthT2_gg/susyEvents_165_1_LBp.root");
//   chain.Add("dcap:///pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/stealthT2_gg/susyEvents_17_2_bXu.root");
//   chain.Add("dcap:///pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/stealthT2_gg/susyEvents_198_1_OJk.root");
//   chain.Add("dcap:///pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/stealthT2_gg/susyEvents_230_1_pAO.root");
//   chain.Add("dcap:///pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/stealthT2_gg/susyEvents_269_1_nZi.root");
//   chain.Add("dcap:///pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/stealthT2_gg/susyEvents_284_1_NgZ.root");
//   chain.Add("dcap:///pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/stealthT2_gg/susyEvents_322_1_a5Q.root");
//   chain.Add("dcap:///pnfs/cms/WAX/11/store/user/lpcsusystealth/ntuples/stealthT2_gg/susyEvents_68_1_Frc.root");

  // Disabling unused branches will speed up the processing significantly, but risks inconsistencies if wrong trees are turned off.
  // Make sure you know what you are doing.
  // Especially be careful when producing a skim - disabled branches will not be copied.
  // You can either enable all branches, or have two Event objects (one for event selection with reduced branch configuration, and the
  // other for event copying with branches fully enabled).
  chain.SetBranchStatus("*", 0);
  chain.SetBranchStatus("runNumber", 1);
  chain.SetBranchStatus("luminosityBlockNumber", 1);
  chain.SetBranchStatus("eventNumber", 1);
  chain.SetBranchStatus("isRealData", 1);
  chain.SetBranchStatus("metFilterBit", 1);
  chain.SetBranchStatus("rho", 1);
  chain.SetBranchStatus("rho25", 1);
  chain.SetBranchStatus("hlt*", 1);
  chain.SetBranchStatus("vertices*", 1);
  chain.SetBranchStatus("tracks*", 1);
  chain.SetBranchStatus("photons_photons*", 1);
  chain.SetBranchStatus("muons_muons*", 1);
  chain.SetBranchStatus("electrons_gsfElectrons*", 1);
  chain.SetBranchStatus("pfJets_ak5*", 1);
  chain.SetBranchStatus("met_pfType01CorrectedMet*", 1);
  chain.SetBranchStatus("gridParams*", 1);

  if(chain.LoadTree(0) != 0){
    cerr << "Error with input chain. Do the files exist?" << endl;
    return;
  }

  SusyEventAnalyzer sea(chain);

  sea.SetOutput(outputName);
  sea.SetLogFile("cout"); // set to a full path to a file to output log to a file
  sea.SetPrintInterval(10000);
  sea.SetPrintLevel(0);
//   sea.AddHltName("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50");
  sea.CopyEvents(false);
  sea.SetProcessNEvents(-1);
  sea.SetNScaledEvents(19.712, -1.0);  // (lumi, x-sec). If no given x-sec, set to -1.0;
  sea.SetSignalDatasetName("stealthT2_gg");  // Set this if running over signal

  TStopwatch ts;

  ts.Start();

  sea.Run();

  ts.Stop();

  std::cout << "RealTime : " << ts.RealTime()/60.0 << " minutes" << std::endl;
  std::cout << "CPUTime  : " << ts.CpuTime()/60.0 << " minutes" << std::endl;

}
