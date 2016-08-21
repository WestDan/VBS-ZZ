#ifndef MyAnalysis_MyxAODAnalysis_H
#define MyAnalysis_MyxAODAnalysis_H

#include <EventLoop/Algorithm.h>

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "AnalysisVar.h"


#include "PATInterfaces/SystematicRegistry.h"
#include "PATInterfaces/SystematicCode.h"
#include "PATInterfaces/CorrectionCode.h"
#include "PATInterfaces/SystematicsUtil.h"
#include "PATInterfaces/SystematicVariation.h" 


#include <TTree.h>
#include <TH1.h>
#include <TLorentzVector.h>
#include <time.h>
#include <TMacro.h>
#include "GoodRunsLists/GoodRunsListSelectionTool.h"
using namespace std;

// muon calibration and smearing tool
class AsgElectronLikelihoodTool;
class OverlapRemovalTool;

namespace ana{
     class QuickAna;
}

class MyxAODAnalysis : public EL::Algorithm, public AnalysisVar  // no definition even in AnalysisVar.h ??? --Zhang
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:

  time_t start, end; //!

  xAOD::TEvent *m_event;  //!
  //xAOD::EventInfo* eventInfo; //!
  // count number of events
  bool isMC;
  int m_eventCounter; //!
  int m_4leptons; // events with 4 leptons
  int m_2pairs;	// evnets with 4 leptons and 2 pairs
  int m_mll1;  // 4 lep + 2 pairs + 1 M.Z
  int m_mll2;	// 4lep + 2pairs + 2M.Z
  int m_numjets; // 4lep + 2pairs + 2M.Z + 2 jets

  int m_truth_xAOD;
  int m_truth_4leptons;
  int m_truth_leppt;
  int m_truth_lepeta;
  int m_truth_leadingleppt;
  int m_truth_2pairs;
  int m_truth_mll1;
  int m_truth_mll2;
  int m_truth_numjets;

  // weight of truth events
  double w_truth_xAOD;  
  double w_truth_4leptons;
  double w_truth_leppt;
  double w_truth_lepeta;
  double w_truth_leadingleppt;
  double w_truth_2pairs;
  double w_truth_mll1;
  double w_truth_mll2;
  double w_truth_numjets;
  

  double m_sumOfWeights;  // sumofWeights
  double m_filter; //!
  int m_fidu1; //!
  int m_fidu2; //!
  int m_fidu3; //!
  int m_1fwdE; //!

  float sumOfWeights;
  float sumOfWeightsSquared;

  unsigned int nInputElectrons; //!
  unsigned int nOverlapElectrons; //!

  string set; // treename
  string m_FileName; //!
  //bool passMuon, passEle, passFwdEle;
  
  ana::QuickAna *quickAna; //!
  std::vector<CP::SystematicSet> m_sysList; //!
  std::vector<CP::SystematicSet> pileup_sysList; //!
  CP::SystematicSet NomSet;//!

#ifndef __CINT__
//  CP::MuonCalibrationAndSmearingTool *m_muonCalibrationAndSmearingTool; //! 
//  CP::MuonSelectionTool *m_muonSelectionTool; //!
  AsgElectronLikelihoodTool* myLikelihood; //!
  OverlapRemovalTool *orTool;//!
#endif // not __CINT__

  GoodRunsListSelectionTool *m_grl; //!
  

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  // Tree *myTree; //!
  // TH1 *myHist; //!
  //std::vector<string>* m_setSysList; //!

  // this is a standard constructor
  MyxAODAnalysis ();
  MyxAODAnalysis (string );

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  // this is needed to distribute the algorithm to the workers
  ClassDef(MyxAODAnalysis, 1);     // a function? where does it come from ???  --Zhang

  void ClearFlags(MapType2_Int& map);
  void ClearFlags(MapType3_Int& map);
  void ClearWeight(MapType2_Double& map);
  void ClearVariables(MapType2_Int& map);
  void ClearVariables(MapType2_Float& map);
  void ClearVariables(MapType2_VFloat& map);
  void ClearVariables(MapType2_VInt& map);
  void ClearVariables(MapType2_Double& map);
  void ClearVariables(MapType2_Double2D& map);
  void ClearVariables(MapType2_VDouble& map);
  void ClearVariables(MapType2_V2DDouble& map);
  void InitStrVec(vector<string>& out, string in, string de=",");
  void InitObjSTEP(MapType_VString& STEP, string obj, string steps);
  void InitSetting(MapType2_Int& setmap, string setname, string settings);
  void InitVarVec(vector<string>& map, MapType_Double& list, string cut, string settings);
  void CreateCountingMap();

  void InitHistVar(string varlist, int nbin, double xmin, double xmax, string cutstep="");
  void InitTruthHistVar(string varlist, int nbin, double xmin, double xmax, string cutstep="");
  void InitHist2DVar(string varlist, int nxbin, double xmin, double xmax,
                                     int nybin, double ymin, double ymax, string cutstep="");
  void InitHist2DVVar(string varlist, int nxbin, double xmin, double xmax,
                                      int nybin, double ymin, double ymax, string cutstep="");
  void InitTreeVar(string varlist, string type);
  void InitVVar(string varlist, int nbin, double xmin, double xmax, string cutstep="xAOD");
  void AddVarIntoTree(TTree *tree);
  void CreateTreeVarHistoMap(TFile *file);
  void CreateTruthHistoMap(TFile* file);
  void FillHistograms(string sysname);
  void FillTruthHistograms(string sysname);
};

#endif
