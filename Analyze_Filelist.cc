#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "SusyEvent.h"
#include "SusyEventAnalyzer.h"
#include "SusyNtuplizer_LinkDef.h"

#include <stdio.h>
#include <TROOT.h>
#include <TTree.h>
#include <TChain.h>
#include <TSystem.h>



using namespace std;



//====================
const int NPARAM = 7;   //--> NUMBER OF PARAMETERS
//====================



int
main(int argc, char* argv[])
{
  gSystem->Load("libSusyEvent.so");

  //==================== Check number of input parameters
  if (argc < NPARAM || argc > NPARAM + 1) {
    cout << "int main() >> ERROR : arcg = " << argc     <<
            " different from "              << NPARAM+1 << endl;

    cout << "--- Correct Usage : exe inputlist [json.txt]" << endl;
  }
  else {
    TString infile  = argv     [1];
    TString ds      = argv     [argc - 5];
    TString physics = argv     [argc - 4];
    Float_t lumi    = atof(argv[argc - 3]);
    Float_t xSec    = atof(argv[argc - 2]);
    Float_t jec     = atof(argv[argc - 1]);

    char listName[500];
    char pName   [500];
    char treeName[30];

    cout << "Loading Data ............" << endl <<
            "Creating TChain ........." << endl;

    sprintf(listName, argv[1]);
    strcpy (treeName, "susyTree");

    TChain* chain = new TChain(treeName);

    cout << "Selected tree " << treeName << endl;

    ifstream is(listName);

    if (is.good()) {
      cout << "Reading list : " << listName << " ......." << endl;

      while (is.getline(pName, 500, '\n')) {
        if (pName[0] == '#')
          continue;

        chain->Add(pName);
      }

      SusyEventAnalyzer analysis(chain);

      if (argc == NPARAM + 1) {
        analysis.SetUseJson  (true);
        analysis.IncludeAJson(argv[2]);
      }

      TString dataset = physics + "_" + ds;   // dataset name
      TString outdir = "";

      analysis.SetOutdir        (outdir); // For submitting to condor we specify no path, so that the output is written to the worker node.
      analysis.SetDataset       (dataset);
      analysis.SetPrintInterval (1e3);     // print frequency
      analysis.SetPrintLevel    (0);       // print level for event contents
      analysis.SetProcessNEvents(-1);      // number of events to be processed
      analysis.SetUseTrigger    (true);
      analysis.SetFilter        (false);   // filter events passing final cuts
      analysis.SetNScaledEvents (lumi, xSec);
      analysis.SetJecUncert     (jec);

//       analysis.AddHltName("HLT_Photon36_CaloId10_Iso50_Photon22_CaloId10_Iso50");
//       analysis.AddHltName("HLT_Photon36_CaloId10_Iso50_Photon22_R9Id85");
//       analysis.AddHltName("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50");
//       analysis.AddHltName("HLT_Photon36_R9Id85_Photon22_CaloId10_Iso50");
//       analysis.AddHltName("HLT_Photon36_R9Id85_Photon22_R9Id85");

      analysis.AddHltName("HLT_HT750");
      analysis.AddHltName("HLT_PFNoPUHT650");

      cout << "Starting Loop" << endl;

      analysis.Loop();
    }
    else
      cout << "int main() >> ERROR : file " << listName << " not read" << endl;

    is.close();
  }
}
