#include <iostream>
#include <fstream>
#include "TSystem.h"
using namespace std;

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>

#include "xAODMuon/Muon.h"  // only muon needed??? --Zhang

#include <QuickAna/QuickAna.h>

#include "MyAnalysis/MyxAODAnalysis.h"
#include "MyAnalysis/AnalysisVar.h"

#include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"

#include "AssociationUtils/OverlapRemovalTool.h"
#include "PileupReweighting/PileupReweightingTool.h"

#include "xAODCutFlow/CutBookkeeper.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"

//ort::inputAccessor_t selectAcc("selected");
//ort::inputDecorator_t selectDec("selected");
//ort::outputAccessor_t overlapAcc("overlaps");
//ort::objLinkAccessor_t objLinkAcc("overlapObject");

EL::StatusCode MyxAODAnalysis :: initialize ()
{
    // Here you do everything that you need to do after the first input
    // file has been connected and before the first event is precessed,
    // e.g. create additional histograms based on which variables are
    // available in the input files. You can also create all of your
    // histograms and trees in here, but be aware that this method
    // doesn't get called if no events are processed. So any method
    // you create here won't be available in the output if you have no 
    // input events.   --Zhang

  
  time(&start);


  SETNAME.push_back("physics");
  // what's the other possible value of SETNAME;

  InitSetting(SETTING,"physics","docorr=1,dosys=0,doweight=1,dovbf=1,doopt=0");
  // SETTING["physics"]["docorr"] = 1
  // SETTING["physics"]["dosys"] = 0
  // SETTING["physics"]["doweight"] = 1
  // SETTING["physics"]["dovbf"] = 1
  // SETTING["physics"]["doopt"] = 0

  InitStrVec(CHN,"eeee,eemm,mmmm",",");
  InitStrVec(STEP_cut,"xAOD,JetClean,Trigger,4Leptons,2Pairs,Mll1,Mll2,NumJets", ",");
  InitStrVec(TRUTH_STEP_cut,"xAOD,4Leptons,LepPt,LepEta,LeadingLepPt,2Pairs,Mll1,Mll2,NumJets", ",");
  InitStrVec(STEP_cut_tree,"M2Lep,MetCut,DLepR,DMetZ,FracPt,BJetVeto");

  // --> original simple step cut for obj mu
  InitObjSTEP(STEP_obj,"mu","All,Tool,CB,Eta,Pt,Z0,D0,TrkIso,OverLap,Pt15,Medium");
  // --> original simple step cut for obj ele
  InitObjSTEP(STEP_obj,"ele","All,ID,LHID,Pt,Eta,ObjQ,Z0,D0,TrkIso,OverLap,Pt15,Medium");
  InitObjSTEP(STEP_obj,"jet","All,Pt,Eta,JVT,OverLap,Clean");
  

  // Varibale dependency: 
  // <--    Pt1..Pt4 ====> 4Leptons	-->
  // <--    Mll1,Mll2 ====> xAOD	-->(choose the two pairs with the most smallest Mass difference to Z)
  // <--    leptons number ====> xAOD  	-->
  // <--    ZZ ====> 4Leptons	  	-->
  InitTruthHistVar("truth_lep_num,truth_ele_num,truth_muon_num", 7, 0, 7, "xAOD,4Leptons,LepPt,LepEta,LeadingLepPt,2Pairs,Mll1,Mll2,NumJets");
  InitTruthHistVar("truth_Pt1,truth_Pt2,truth_Pt3,truth_Pt4", 24, 0, 240.e3, "4Leptons,LepPt,LepEta,LeadingLepPt,2Pairs,Mll1,Mll2,NumJets");
  InitTruthHistVar("truth_Eta1,truth_Eta2,truth_Eta3,truth_Eta4", 30, -3.0, 3.0, "4Leptons,LepPt,LepEta,LeadingLepPt,2Pairs,Mll1,Mll2,NumJets");
  InitTruthHistVar("truth_Phi1,truth_Phi2,truth_Phi3,truth_Phi4", 35, -3.5, 3.5, "4Leptons,LepPt,LepEta,LeadingLepPt,2Pairs,Mll1,Mll2,NumJets");
  InitTruthHistVar("truth_Z1,truth_Z2", 25, 66.e3, 116.e3, "xAOD,4Leptons,LepPt,LepEta,LeadingLepPt,2Pairs,Mll1,Mll2,NumJets");
  InitTruthHistVar("truth_ZZ", 50, 0, 500.e3, "4Leptons,LepPt,LepEta,LeadingLepPt,2Pairs,Mll1,Mll2,NumJets");
  InitTruthHistVar("truth_NumJets", 5, 0, 5, "xAOD,4Leptons,LepPt,LepEta,LeadingLepPt,2Pairs,Mll1,Mll2,NumJets");
  InitTruthHistVar("truth_leadingJet_Pt,truth_subleadingJet_Pt", 50, 0, 500.e3, "xAOD,4Leptons,LepPt,LepEta,LeadingLepPt,2Pairs,Mll1,Mll2,NumJets");
  InitTruthHistVar("truth_leadingJet_Eta,truth_subleadingJet_Eta", 40, -4.0, 4.0, "xAOD,4Leptons,LepPt,LepEta,LeadingLepPt,2Pairs,Mll1,Mll2,NumJets");
  InitTruthHistVar("truth_leadingJet_Phi,truth_subleadingJet_Phi", 35, -3.5, 3.5, "xAOD,4Leptons,LepPt,LepEta,LeadingLepPt,2Pairs,Mll1,Mll2,NumJets");
  InitTruthHistVar("truth_Mjj", 50, 0, 1000.e3, "xAOD,4Leptons,LepPt,LepEta,LeadingLepPt,2Pairs,Mll1,Mll2,NumJets");
  InitTruthHistVar("truth_Delta_Jet_Eta", 40, 0, 10, "xAOD,4Leptons,LepPt,LepEta,LeadingLepPt,2Pairs,Mll1,Mll2,NumJets");
  InitTruthHistVar("truth_Centrality", 40, -10, 10, "xAOD,4Leptons,LepPt,LepEta,LeadingLepPt,2Pairs,Mll1,Mll2,NumJets");

  InitHistVar("Pt1,Pt2,Pt3,Pt4", 24, 0, 240.e3, "All");
  InitHistVar("Eta1,Eta2,Eta3,Eta4", 20, -3.0, 3.0, "All");
  InitHistVar("Phi1,Phi2,Phi3,Phi4", 20, -3.5, 3.5, "All");
  InitHistVar("NumJets", 10, 0, 10, "All");
  InitHistVar("LeadingJetPt,SubLeadingJetPt", 30, 0, 300.e3, "All");
  InitHistVar("LeadingJetEta,SubLeadingJetEta", 20, -4.0, 4.0, "All");
  InitHistVar("LeadingJetPhi,SubLeadingJetPhi", 20, -3.5, 3.5, "All");
  InitHistVar("DiJetMass", 50, 0, 1000.e3, "All");
  InitHistVar("DeltaJetEta", 40, 0, 10, "All");
  InitHistVar("DiLepMass1", 50, 40, 140.e3, "Mll1,Mll2,NumJets");
  InitHistVar("DiLepMass2", 50, 40, 140.e3, "Mll2,NumJets");
  InitHistVar("DiLepPt1", 30, 0, 300.e3, "Mll1,Mll2,NumJets");
  InitHistVar("DiLepPt2", 30, 0, 300.e3, "Mll2,NumJets");
  InitHistVar("DiLepRap1", 15, 0, 3, "Mll1,Mll2,NumJets");
  InitHistVar("DiLepRap2", 15, 0, 3, "Mll2,NumJets");
  InitHistVar("4LepMass",30, 0, 600.e3,"Mll2,NumJets");
  InitHistVar("4LepPt", 30, 0, 300.e3, "Mll2,NumJets");
  InitHistVar("Centrality", 40, -10, 10, "Mll2,NumJets");

  m_event = wk()->xaodEvent();
  Info("initialize()", "Number of events = %lli", m_event->getEntries() );
  m_FileName = wk()->inputFile()->GetName();
  cout<<endl<<endl<<endl<<"File name: "<<m_FileName<<endl;

  const xAOD::EventInfo* eventInfo = 0;
  if( !m_event->retrieve( eventInfo, "EventInfo").isSuccess() )
  {
    Error("execute ()", "Failed to retrieve EventInfo. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  isMC = eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION );
  // eventType :
  // IS_SIMULATION | true: simulation, false: data
  // IS_TESTBEAM   | true: testbeam,   false: full detector
  // IS_CALIBRATION| true: calibration,false: physics

  uint64_t nEventsProcessed=0;
  sumOfWeights = 0.;
  sumOfWeightsSquared = 0.;

  if(isMC) {

    TTree *MetaData = dynamic_cast<TTree*>(wk()->inputFile()->Get("MetaData"));

    if (!MetaData) {
      Error("fileExecute()", "MetaData not found! Exiting.");
      return EL::StatusCode::FAILURE;
    }

    MetaData->LoadTree(0);
    bool m_isDerivation = !MetaData->GetBranch("StreamAOD");
    if(m_isDerivation ){

      const xAOD::CutBookkeeperContainer* incompleteCBC = nullptr;
      if(!m_event->retrieveMetaInput(incompleteCBC, "IncompleteCutBookkeepers").isSuccess()){
        Error("initializeEvent()","Failed to retrieve IncompleteCutBookkeepers from MetaData! Exiting.");
        return EL::StatusCode::FAILURE;
      }
      if ( incompleteCBC->size() != 0 ) {
        Error("initializeEvent()","Found incomplete Bookkeepers! Check file for corruption.");
        return EL::StatusCode::FAILURE;
      }

      const xAOD::CutBookkeeperContainer* completeCBC = 0;
      if(!m_event->retrieveMetaInput(completeCBC, "CutBookkeepers").isSuccess()){
        Error("initializeEvent()","Failed to retrieve CutBookkeepers from MetaData! Exiting.");
        return EL::StatusCode::FAILURE;
      }
      int minCycle = 10000;

      for ( auto cbk : *completeCBC ) {
        if ( ! cbk->name().empty()  && minCycle > cbk->cycle() ){ minCycle = cbk->cycle(); }
      }

      const xAOD::CutBookkeeper* allEventsCBK=0;
      for ( auto cbk :  *completeCBC ) {
        if ( minCycle == cbk->cycle() && cbk->name() == "AllExecutedEvents" ){
          allEventsCBK = cbk;
          break;
        }
      }

      nEventsProcessed  = allEventsCBK->nAcceptedEvents();
      sumOfWeights        = allEventsCBK->sumOfEventWeights();
      sumOfWeightsSquared = allEventsCBK->sumOfEventWeightsSquared();
    }
  }
  cout<<" nEventsProcessed: "<< nEventsProcessed <<" sumOfWeights: "<<sumOfWeights<<" sumOfWeightsSquared: " << sumOfWeightsSquared<<endl;
  TFile *f = new TFile ("SumOfWeights.root","RECREATE");
  TH1F *h1 = new TH1F("sumOfWeights", "sumOfWeights", 0,10,10);
  h1->Fill(5, sumOfWeights);
  h1->Write();
  f->Close();


  //orTool = new OverlapRemovalTool("OverlapRemovalTool");
  //orTool->setProperty("InputLabel", "selected");
  //orTool->setProperty("LinkOverlapObjects", true);
  //orTool->initialize();

  m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
  const char* GRLFilePath = "/afs/cern.ch/user/w/weibin/ZZ/Test/ExtraTools/grl/data15_13TeV.periodAllYear_DetStatus-v73-pro19-08_DQDefects-00-01-02_PHYS_StandardGRL_All_Good_25ns.xml";
  const char* fullGRLFilePath = gSystem->ExpandPathName (GRLFilePath);
  std::vector<std::string> vecStringGRL;
  vecStringGRL.push_back(fullGRLFilePath);
  m_grl->setProperty( "GoodRunsListVec", vecStringGRL);
  m_grl->setProperty("PassThrough", false); // if true (default) will ignore result of GRL and will just pass all events
  m_grl->initialize();

  std::vector<std::string> prwFiles;
  //string PileConfigFile = "/atlas/data19/wenguo/tools_and_files/prw/MC15_ttbar25ns_prw.root";
  string PileConfigFile = "/afs/cern.ch/user/w/weibin/ZZ/Test/ExtraTools/prw/MC15b.prw.ttbar.root";        
  prwFiles.push_back(PileConfigFile);
  std::vector<std::string> lumicalcFiles;
  string LumiCalcFile3 = "/afs/cern.ch/user/w/weibin/ZZ/Test/ExtraTools/prw/ilumicalc_histograms_None_276262-284484.root";
  lumicalcFiles.push_back(LumiCalcFile3);


  quickAna =new ana::QuickAna("quickana");
  if(!isMC) quickAna->isDataFlag="true";
  quickAna->eventinfoDef = "default";
//  if(!isMC) quickAna->triggerDef = "HLT_mu20_iloose_L1MU15 HLT_mu50 HLT_e24_lhmedium_L1EM20VH HLT_e60_lhmedium HLT_e120_lhloose";
//  else quickAna->triggerDef = "HLT_mu20_iloose_L1MU15 HLT_mu50 HLT_e24_lhmedium_L1EM18VH HLT_e60_lhmedium HLT_e120_lhloose";
  if(!isMC) quickAna->triggerDef = "HLT_mu20_iloose_L1MU15_OR_HLT_mu50 HLT_e24_lhmedium_L1EM20VH HLT_e60_lhmedium HLT_e120_lhloose";
  else quickAna->triggerDef = "HLT_mu20_iloose_L1MU15_OR_HLT_mu50 HLT_e24_lhmedium_L1EM18VH HLT_e60_lhmedium HLT_e120_lhloose";
  quickAna->muonDef="hzhinv_loose hzhinv_medium";
  quickAna->electronDef="hzhinv_loose hzhinv_medium";
  quickAna->jetKine = "pt > 20e3 && eta > -4.5 && eta < 4.5";
  quickAna->jetDef="antikt04_HZZ";
  quickAna->tauDef="none";
  quickAna->photonDef="none";
  quickAna->metDef="metZHinv";
  quickAna->orDef = "default";
  quickAna->muDataFiles=lumicalcFiles;
  quickAna->muMcFiles=prwFiles;
  quickAna->initialize().ignore();

  m_sysList= CP::make_systematics_vector (quickAna->recommendedSystematics());
  CP::SystematicCode::enableFailure();

  
  SYSNAME.push_back("NOCORR");
  for (auto sysListItr : m_sysList){
    cout << "name sys is : " << sysListItr.name() << endl;
    if(sysListItr.name()=="") {
      SYSNAME.push_back("NOMINAL");
    }
    else SYSNAME.push_back(sysListItr.name());
  }
  // according to the output file: 
  // we can know that 
  // SYSNAME = { "NOCORR",
  //		 "NOMINAL"
  //		 "EG_..."
  //		 "EL_..."
  //		 "JET_..."
  //		 "MET_..."
  //		 "MUON_..."
  //	       }

  TFile *file1 = wk()->getOutputFile ("tree_output");
  if(SETTING["physics"]["dobkgwz"]) { // ???
    CHN.clear();
    STEP_cut.clear();
    InitStrVec(CHN,"mmm,mme,eee,eem",",");    
    InitStrVec(STEP_cut,"xAOD,V3rdLep,OS");

    tree_bkg = new TTree("tree_bkg", "tree_bkg");
    tree_bkg->SetDirectory (file1);
    AddVarIntoTree(tree_bkg);
  }else if(SETTING["physics"]["dobkgz"]) {
    STEP_cut.clear();
    InitStrVec(STEP_cut,"xAOD,TwoMoreLep,OS,V3rdLep,M2Lep,MetCut,FracPt,NoJet");  
    tree_bkg = new TTree("tree_bkg", "tree_bkg");
    tree_bkg->SetDirectory (file1);
    AddVarIntoTree(tree_bkg);
  }else {
    for(int i=0; i<(int)SYSNAME.size(); i++) {

      string sys=SYSNAME[i];

      if(sys!="NOCORR" && (!SETTING["physics"]["docorr"])) continue;

      if(sys!="NOMINAL" && SETTING["physics"]["docorr"] && (!SETTING["physics"]["dosys"])) continue;

      if(sys=="NOCORR" && SETTING["physics"]["docorr"] && SETTING["physics"]["dosys"]) continue;

      string tree_name = "tree_" + sys;

      Tree[sys] = new TTree(tree_name.c_str(), "output tree");
      Tree[sys]->SetDirectory (file1);
      AddVarIntoTree(Tree[sys]);
    }
  }

  CreateCountingMap();

  TFile *file2 = wk()->getOutputFile ("hist_output");

  CreateTreeVarHistoMap(file2);

  TFile *file_truth = wk()->getOutputFile("truth_hist_output");
  CreateTruthHistoMap(file_truth);

  m_eventCounter = 0;
  m_4leptons = 0;
  m_2pairs = 0;
  m_mll1 = 0;
  m_mll2 = 0;
  m_numjets = 0;

  m_truth_xAOD = 0;
  m_truth_4leptons = 0;
  m_truth_leppt = 0;
  m_truth_lepeta = 0;
  m_truth_leadingleppt = 0;
  m_truth_2pairs = 0;
  m_truth_mll1 = 0;
  m_truth_mll2 = 0;
  m_truth_numjets = 0;

  w_truth_xAOD = 0.0;
  w_truth_4leptons = 0.0;
  w_truth_leppt = 0.0;
  w_truth_lepeta = 0.0;
  w_truth_leadingleppt = 0.0;
  w_truth_2pairs = 0.0;
  w_truth_mll1 = 0.0;
  w_truth_mll2 = 0.0;
  w_truth_numjets = 0.0;

  m_sumOfWeights = 0.0;
  m_filter = 0;
  m_fidu1 = 0;
  m_fidu2 = 0;
  m_fidu3 = 0;
  m_1fwdE = 0;

  nInputElectrons = 0;
  nOverlapElectrons = 0;

  return EL::StatusCode::SUCCESS;

}


void MyxAODAnalysis :: ClearFlags(MapType2_Int& map) {

    MapType2_Int::iterator it;
    for(it=map.begin(); it!=map.end(); it++) {
        string chn=(*it).first;
        MapType_Int::iterator it2;
        for(it2=(*it).second.begin(); it2!=(*it).second.end(); it2++) {
            string cut=(*it2).first;
            map[chn][cut] = 0;
        }
    }
} // isn't there any method to initialize a map container object to specified value ??? --Zhang

void MyxAODAnalysis :: ClearFlags(MapType3_Int& map) {

    MapType3_Int::iterator it;
    for(it=map.begin(); it!=map.end(); it++) {
      string chn=(*it).first;
      MapType2_Int::iterator it2;
      for(it2=(*it).second.begin(); it2!=(*it).second.end(); it2++) {
        string cut=(*it2).first;
        MapType_Int::iterator it3;
        for(it3=(*it2).second.begin(); it3!=(*it2).second.end(); it3++) {
          string var=(*it3).first;
          map[chn][cut][var] = 0;
        }
      }
    }
}


void MyxAODAnalysis :: ClearWeight(MapType2_Double& map) {
    MapType2_Double::iterator it;

    MapType_Double::iterator iit;

    for(it=map.begin(); it!=map.end(); it++)
      for(iit=(*it).second.begin(); iit!=(*it).second.end(); iit ++)
        map[(*it).first][(*iit).first]=1.0;
}



void MyxAODAnalysis :: CreateCountingMap() {

    MapType_VString::iterator it;
    for(it=STEP_obj.begin(); it!=STEP_obj.end(); it++) {
      string obj=(*it).first;
      for(int i=0; i<(int)(*it).second.size(); i++) {
        string cut=STEP_obj[obj][i];
        COUNT ini={0.,0.,0.};
	// COUNT { float num, wnum, err;};

        for(int j=0; j<(int)SYSNAME.size(); j++) {
          CNT_obj[SYSNAME[j]][obj][cut]=ini;
        }

      }
    }

    for(int i=0; i<(int)CHN.size(); i++) {
      string chn=CHN[i];

      for(int j=0; j<(int)STEP_cut.size(); j++) {
        string cut=STEP_cut[j];
        Evt_Weight[chn][cut]=1.0;
        FLAG_cut_temp[chn][cut]["default"]=0;
        FLAG_cut[chn][cut]["default"]=0; 

        COUNT ini={0.,0.,0.};

          for(int k=0; k<(int)SYSNAME.size(); k++) {
             CNT_cut[SYSNAME[k]][chn][cut]["default"]=ini;
          }
      }
    }
}


void MyxAODAnalysis :: InitStrVec(vector<string>& out, string in, string de) {
    int pos=0, pos_pre=0;
    while(true) {
        pos=in.find(de,pos_pre);
        if(pos==-1) {out.push_back(in.substr(pos_pre,in.size()-pos_pre)); break;}
        else  out.push_back(in.substr(pos_pre,pos-pos_pre)); // this depend on that no space between literal word and delimiter.
        pos_pre=pos+1;
    }
}



void MyxAODAnalysis :: InitObjSTEP(MapType_VString& STEP, string obj, string steps) {
    vector<string> str, objstr;
    InitStrVec(str, steps, ",");

    for(int i=0; i<(int)str.size(); i++) {
      string objs = str[i];
      objstr.push_back(objs);
    } // is there any difference between str and objstr ??? --Zhang

    STEP[obj]=objstr;
}

void MyxAODAnalysis :: ClearVariables(MapType2_Int& map) {

    MapType2_Int::iterator it;
    for(it=map.begin(); it!=map.end(); it++) {
      string varname = (*it).first;
      map[varname]["Value"]=-9999;
    }
}

void MyxAODAnalysis :: ClearVariables(MapType2_Float& map) {

    MapType2_Float::iterator it;
    for(it=map.begin(); it!=map.end(); it++) {
      string varname = (*it).first;
      map[varname]["Value"]=-9999.;
    }
}

void MyxAODAnalysis :: ClearVariables(MapType2_VFloat& map) {

    MapType2_VFloat::iterator it;
    for(it=map.begin(); it!=map.end(); it++) {
      string varname = (*it).first;
      map[varname]["Value"].clear();
    }
}

void MyxAODAnalysis :: ClearVariables(MapType2_VInt& map) {

    MapType2_VInt::iterator it;
    for(it=map.begin(); it!=map.end(); it++) {
      string varname = (*it).first;
      map[varname]["Value"].clear();
    }
}



void MyxAODAnalysis :: ClearVariables(MapType2_Double& map) {

    MapType2_Double::iterator it;
    for(it=map.begin(); it!=map.end(); it++) {
      string varname = (*it).first;
      map[varname]["Value"]=-9999.;
    }
}

void MyxAODAnalysis :: ClearVariables(MapType2_Double2D& map) {

    MapType2_Double2D::iterator it;
    pair<double, double> init (-9999.0, -9999.0);
    for(it=map.begin(); it!=map.end(); it++) {
      string varname = (*it).first;
      map[varname]["Value"] = init;
    }
}



void MyxAODAnalysis :: ClearVariables(MapType2_VDouble& map) {

    MapType2_VDouble::iterator it;
    for(it=map.begin(); it!=map.end(); it++) {
      string varname = (*it).first;
      map[varname]["Value"].clear();
      map[varname]["Weight"].clear();
    }
}

void MyxAODAnalysis :: ClearVariables(MapType2_V2DDouble& map) {

    MapType2_V2DDouble::iterator it;
    for(it=map.begin(); it!=map.end(); it++) {
      string varname = (*it).first;
      map[varname]["Value"].clear();
      map[varname]["NBins"].clear();
      map[varname]["Min"].clear();
      map[varname]["Max"].clear();
      map[varname]["Weight"].clear();
    }
}



void MyxAODAnalysis :: InitTreeVar(string varlist, string type) {
    vector<string> variables;
    InitStrVec(variables, varlist, ",");

    if(type=="I") {
      for(int i=0; i<(int)variables.size(); i++) {
        TreeIntVar[variables[i]]["Value"]=-9999;
      }
    }
    if(type=="U") {
      for(int i=0; i<(int)variables.size(); i++) {
        TreeUloVar[variables[i]]["Value"]=0;
      }
    }
    if(type=="F") {
      for(int i=0; i<(int)variables.size(); i++) {
        TreeFltVar[variables[i]]["Value"]=-9999.0;
      }
    }
    if(type=="S") {
      for(int i=0; i<(int)variables.size(); i++) {
        TreeStrVar[variables[i]]["Value"]=" ";
      }
    }
    if(type=="vector<F>") {
      for(int i=0; i<(int)variables.size(); i++) {
        TreeFltVVar[variables[i]]["Value"].clear();
      }
    }

    if(type=="vector<I>") {
      for(int i=0; i<(int)variables.size(); i++) {
        TreeIntVVar[variables[i]]["Value"].clear();
      }
    }


}


void MyxAODAnalysis :: InitHistVar(string varlist, int nbin, double xmin, double xmax, string cutstep) {
    vector<string> variables;
    InitStrVec(variables, varlist, ","); 
    // actually, this sentense is superfluous, because rarely, in 
    // reality, will we use the same nbin, xmin and xmax for 
    // different variables


    vector<string> cuts;
    if(cutstep=="") cuts.clear();
    else if(cutstep=="All") cuts = STEP_cut;
    else if(cutstep.compare(0,1,"-")==0) {
      cutstep = cutstep.erase(0,1);
      for(int i=0; i<(int)STEP_cut.size(); i++) {
        if(cutstep != STEP_cut[i]) cuts.push_back(STEP_cut[i]);
        else if(cutstep == STEP_cut[i]) {cuts.push_back(STEP_cut[i]); break;}
      }
    }else {
      InitStrVec(cuts, cutstep, ",");
    }

    for(int i=0; i<(int)variables.size(); i++) {

      HistVar[variables[i]]["Value"]=-9999.0;
      HistVar[variables[i]]["NBins"]=nbin;
      HistVar[variables[i]]["Xmin"]=xmin;
      HistVar[variables[i]]["Xmax"]=xmax;
      HistVar[variables[i]]["Weight"]=1.0;
   
      for(int j=0; j<(int)cuts.size(); j++) {
        HistVar[variables[i]][cuts[j]]=1;
      }
    }
}


void MyxAODAnalysis :: InitTruthHistVar(string varlist, int nbin, double xmin, double xmax, string cutstep) {
    vector<string> variables;
    InitStrVec(variables, varlist, ","); 

    vector<string> cuts;
    if(cutstep=="") cuts.clear();
    else if(cutstep=="All") cuts = TRUTH_STEP_cut;
    else if(cutstep.compare(0,1,"-")==0) {
      cutstep = cutstep.erase(0,1);
      for(int i=0; i<(int)TRUTH_STEP_cut.size(); i++) {
        if(cutstep != TRUTH_STEP_cut[i]) cuts.push_back(TRUTH_STEP_cut[i]);
        else if(cutstep == TRUTH_STEP_cut[i]) {cuts.push_back(TRUTH_STEP_cut[i]); break;}
      }
    }else {
      InitStrVec(cuts, cutstep, ",");
    }

    for(int i=0; i<(int)variables.size(); i++) {

      TruthHistVar[variables[i]]["Value"]=-9999.0;
      TruthHistVar[variables[i]]["NBins"]=nbin;
      TruthHistVar[variables[i]]["Xmin"]=xmin;
      TruthHistVar[variables[i]]["Xmax"]=xmax;
      TruthHistVar[variables[i]]["Weight"]=1.0;
   
      for(int j=0; j<(int)cuts.size(); j++) {
        TruthHistVar[variables[i]][cuts[j]]=1;
      }
    }
}

void MyxAODAnalysis :: InitHist2DVar(string varlist, int nxbin, double xmin, double xmax,int nybin, double ymin, double ymax, string cutstep) {
    vector<string> variables;
    InitStrVec(variables, varlist, ",");

    vector<string> cuts;
    if(cutstep=="") cuts.clear();
    else if(cutstep=="All") cuts = STEP_cut;
    else if(cutstep.compare(0,1,"-")==0) {
      cutstep = cutstep.erase(0,1);
      for(int i=0; i<(int)STEP_cut.size(); i++) {
        if(cutstep != STEP_cut[i]) cuts.push_back(STEP_cut[i]);
        else if(cutstep == STEP_cut[i]) {cuts.push_back(STEP_cut[i]); break;}
      }
    }else {
      InitStrVec(cuts, cutstep, ",");
    }

    pair<double, double> init (-9999.0, -9999.0);
    for(int i=0; i<(int)variables.size(); i++) {

      Hist2DVar[variables[i]]["Value"]=init;
      Hist2DVar[variables[i]]["NBins"]=make_pair(nxbin, nybin);
      Hist2DVar[variables[i]]["Min"]=make_pair(xmin, ymin);
      Hist2DVar[variables[i]]["Max"]=make_pair(xmax, ymax);

      for(int j=0; j<(int)cuts.size(); j++) {
        Hist2DVar[variables[i]][cuts[j]]=make_pair(1,1);
      }
    }
}


void MyxAODAnalysis :: InitVVar(string varlist, int nbin, double xmin, double xmax, string cutstep) {
    vector<string> variables;
    InitStrVec(variables, varlist, ",");

    vector<string> cuts;
    if(cutstep=="") cuts.clear();
    else if(cutstep=="All") cuts = STEP_cut;
    else if(cutstep.compare(0,1,"-")==0) {
      cutstep = cutstep.erase(0,1);
      for(int i=0; i<(int)STEP_cut.size(); i++) {
        if(cutstep != STEP_cut[i]) cuts.push_back(STEP_cut[i]);
        else if(cutstep == STEP_cut[i]) {cuts.push_back(STEP_cut[i]); break;}
      }
    }else {
      InitStrVec(cuts, cutstep, ",");
    }

    for(int i=0; i<(int)variables.size(); i++) {

      VVar[variables[i]]["Value"].clear();
      VVar[variables[i]]["NBins"].push_back(nbin);
      VVar[variables[i]]["Xmin"].push_back(xmin);
      VVar[variables[i]]["Xmax"].push_back(xmax);
      VVar[variables[i]]["Weight"].clear();
  
      for(int j=0; j<(int)cuts.size(); j++) {
        VVar[variables[i]][cuts[j]].push_back(1);
      }
    }
}

void MyxAODAnalysis :: InitHist2DVVar(string varlist, int nxbin, double xmin, double xmax,int nybin, double ymin, double ymax, string cutstep) {
    vector<string> variables;
    InitStrVec(variables, varlist, ",");

    vector<string> cuts;
    if(cutstep=="") cuts.clear();
    else if(cutstep=="All") cuts = STEP_cut;
    else if(cutstep.compare(0,1,"-")==0) {
      cutstep = cutstep.erase(0,1);
      for(int i=0; i<(int)STEP_cut.size(); i++) {
        if(cutstep != STEP_cut[i]) cuts.push_back(STEP_cut[i]);
        else if(cutstep == STEP_cut[i]) {cuts.push_back(STEP_cut[i]); break;}
      }
    }else {
      InitStrVec(cuts, cutstep, ",");
    }

    for(int i=0; i<(int)variables.size(); i++) {

      V2DVar[variables[i]]["Value"].clear();
      V2DVar[variables[i]]["NBins"].push_back(make_pair(nxbin, nybin));
      V2DVar[variables[i]]["Min"].push_back(make_pair(xmin, ymin));
      V2DVar[variables[i]]["Max"].push_back(make_pair(xmax, ymax));

      for(int j=0; j<(int)cuts.size(); j++) {
        V2DVar[variables[i]][cuts[j]].push_back(make_pair(1,1));
      }
    }

}



void MyxAODAnalysis :: AddVarIntoTree(TTree *tree) {
    
    char buf[1000];
    getcwd(buf, sizeof(buf));
    setenv("LOCAL", buf, 1);
    FILE* fp;
    char result_buf[1000];
    fp = popen("find $LOCAL -name MiniTree.txt", "r");
    fgets(result_buf, sizeof(result_buf), fp);

    string varfile(result_buf);
    size_t len = varfile.size();
    varfile.erase(len-1);
    ifstream file;
    file.open(varfile.c_str(), ios::out);

    if (file.is_open())  {
      char line[256];
      while (!file.eof() )  {
        string varname, type;

        file.getline (line,100);
        string sline(line);

        if(sline.find("Int_t")!=string::npos) {
          type = "I";
          varname = sline.substr(9);
          type = varname+"/"+type;
          InitTreeVar(varname,"I");
          tree->Branch(varname.c_str(),&TreeIntVar[varname]["Value"], type.c_str());
        }
        if(sline.find("Ulo_t")!=string::npos) {
          type = "U";
          varname = sline.substr(9);
          type = varname+"/"+type;
          InitTreeVar(varname, "U");
          tree->Branch(varname.c_str(),&TreeUloVar[varname]["Value"]);
        }
        if(sline.find("Float_t")!=string::npos) {
          type = "F";
          varname = sline.substr(9);
          type = varname+"/"+type;
          InitTreeVar(varname, "F");
          tree->Branch(varname.c_str(),&TreeFltVar[varname]["Value"], type.c_str());
        }
        if(sline.find("Str")!=string::npos) {
          type = "S";
          varname = sline.substr(9);
          type = varname+"/"+type;
          InitTreeVar(varname, "S");
          tree->Branch(varname.c_str(),&TreeStrVar[varname]["Value"]);
        }
        if(sline.find("Vector_I")!=string::npos) {
          varname = sline.substr(9);
          InitTreeVar(varname, "Vector_I");
          tree->Branch(varname.c_str(), &TreeIntVVar[varname]["Value"]);
        }
        if(sline.find("Vector_F")!=string::npos) {
          varname = sline.substr(9);
          InitTreeVar(varname, "Vector_F");
          tree->Branch(varname.c_str(), &TreeFltVVar[varname]["Value"]);
        }
      }
    }

}

void MyxAODAnalysis :: CreateTreeVarHistoMap(TFile* file) {
    
  for(int k=0; k<(int)SYSNAME.size(); k++) {
    string sysname = SYSNAME[k];

    if(sysname!="NOCORR" && (!SETTING["physics"]["docorr"])) continue;

    if(sysname!="NOMINAL" && SETTING["physics"]["docorr"] && (!SETTING["physics"]["dosys"])) continue;

    if(sysname=="NOCORR" && SETTING["physics"]["docorr"] && SETTING["physics"]["dosys"]) continue;

    MapType2_Double::iterator it;
    for(it=HistVar.begin(); it!=HistVar.end(); it++) {
      string varname=(*it).first;
      int nbin = int((*it).second["NBins"]);
      double xlow = double((*it).second["Xmin"]);
      double xhigh = double((*it).second["Xmax"]);

      if(nbin == 0) continue;
      if(xlow > xhigh) 
	cout << "Error, in create histogram: " << varname << ", Xmin should less than Xmax." << endl;

      for(int i=0; i<(int)CHN.size(); i++) {
        string chn=CHN[i];

        for(int j=0; j<(int)STEP_cut.size(); j++) {
          string cut=STEP_cut[j];
          if(int((*it).second[cut])!=1) continue;

          string histo_name = sysname + "_" + chn + "_" + cut + "_" + varname;
          TH1D *histo_pointer = new TH1D(histo_name.c_str(),histo_name.c_str(),nbin,xlow,xhigh);
          histo[sysname][chn][cut][varname]=histo_pointer;
          histo[sysname][chn][cut][varname]->SetDirectory(file);
        }
      }
    }

    MapType2_Double2D::iterator it2D;
    for(it2D=Hist2DVar.begin(); it2D!=Hist2DVar.end(); it2D++) {
      string varname=(*it2D).first;

      int nxbin = int((*it2D).second["NBins"].first);
      double xlow = double((*it2D).second["Min"].first);
      double xhigh = double((*it2D).second["Max"].first);
      int nybin = int((*it2D).second["NBins"].second);
      double ylow = double((*it2D).second["Min"].second);
      double yhigh = double((*it2D).second["Max"].second);

      if(nxbin == 0) continue;

      for(int i=0; i<(int)CHN.size(); i++) {
        string chn=CHN[i];

        for(int j=0; j<(int)STEP_cut.size(); j++) {
          string cut=STEP_cut[j];

          pair<double, double> ncuts = (*it2D).second[cut];
          if(int(ncuts.first)!=1) continue;

          string histo_name = sysname + "_" + chn + "_" + cut + "_" + varname;
          TH2F *histo_pointer = new TH2F(histo_name.c_str(),histo_name.c_str(),nxbin,xlow,xhigh,nybin,ylow,yhigh);
          histo_2D[sysname][chn][cut][varname]=histo_pointer;
          histo_2D[sysname][chn][cut][varname]->SetDirectory(file);
        }
      }
    }


    MapType2_VDouble::iterator itv;
    for(itv=VVar.begin(); itv!=VVar.end(); itv++) {
      string varname=(*itv).first;
      int nbin = int((*itv).second["NBins"][0]);
      double xlow = double((*itv).second["Xmin"][0]);
      double xhigh = double((*itv).second["Xmax"][0]);

      if(nbin == 0) continue;

      for(int i=0; i<(int)CHN.size(); i++) {
        string chn=CHN[i];

        for(int j=0; j<(int)STEP_cut.size(); j++) {
          string cut=STEP_cut[j];
          if(int((*itv).second[cut].size())!=1) continue;
          if(int((*itv).second[cut][0])!=1) continue;

          string histo_name = sysname + "_" + chn + "_" + cut + "_" + varname;
          TH1D *histo_pointer = new TH1D(histo_name.c_str(),histo_name.c_str(),nbin,xlow,xhigh);
          histo[sysname][chn][cut][varname]=histo_pointer;
          histo[sysname][chn][cut][varname]->SetDirectory(file);
        }
      }
    }

    MapType2_V2DDouble::iterator it2Dv;
    for(it2Dv=V2DVar.begin(); it2Dv!=V2DVar.end(); it2Dv++) {
        string varname=(*it2Dv).first;

        int nxbin = int((*it2Dv).second["NBins"][0].first);
        double xlow = double((*it2Dv).second["Min"][0].first);
        double xhigh = double((*it2Dv).second["Max"][0].first);
        int nybin = int((*it2Dv).second["NBins"][0].second);
        double ylow = double((*it2Dv).second["Min"][0].second);
        double yhigh = double((*it2Dv).second["Max"][0].second);

        if(nxbin == 0) continue;

        for(int i=0; i<(int)CHN.size(); i++) {
          string chn=CHN[i];

          for(int j=0; j<(int)STEP_cut.size(); j++) {
            string cut=STEP_cut[j];

            if(int((*it2Dv).second[cut].size())!=1) continue;
            pair<double, double> ncuts = (*it2Dv).second[cut][0];
            if(int(ncuts.first)!=1) continue;

            string histo_name = sysname + "_" + chn + "_" + cut + "_" + varname;
            TH2F *histo_pointer = new TH2F(histo_name.c_str(),histo_name.c_str(),nxbin,xlow,xhigh,nybin,ylow,yhigh);
            histo_v2D[sysname][chn][cut][varname]=histo_pointer;
            histo_v2D[sysname][chn][cut][varname]->SetDirectory(file);
          }
        }
      }
  }
}

void MyxAODAnalysis :: CreateTruthHistoMap(TFile* truth_file) {
    
  for(int k=0; k<(int)SYSNAME.size(); k++) {
    string sysname = SYSNAME[k];

    if(sysname!="NOCORR" && (!SETTING["physics"]["docorr"])) continue;

    if(sysname!="NOMINAL" && SETTING["physics"]["docorr"] && (!SETTING["physics"]["dosys"])) continue;

    if(sysname=="NOCORR" && SETTING["physics"]["docorr"] && SETTING["physics"]["dosys"]) continue;

    MapType2_Double::iterator it;
    
    for(it=TruthHistVar.begin(); it!=TruthHistVar.end(); ++it)
    {
      string varname = (*it).first;
      int nbin = int((*it).second["NBins"]);
      double xlow = double((*it).second["Xmin"]);
      double xhigh = double((*it).second["Xmax"]);

      if(nbin == 0) continue;
      if(xlow > xhigh) 
	cout << "Error, in create histogram: " << varname << ", Xmin should less than Xmax." << endl;

      for(int i=0; i<(int)CHN.size(); i++)
      {
	string chn = CHN[i];

	for(int j=0; j<(int)TRUTH_STEP_cut.size(); j++)
	{
	  string cut = TRUTH_STEP_cut[j];
          if(int((*it).second[cut])!=1) continue;

	  string histo_name = sysname + "_" + chn + "_" + cut + "_" + varname;
	  TH1D * histo_pointer = new TH1D(histo_name.c_str(), histo_name.c_str(), nbin, xlow, xhigh);
	  histo[sysname][chn][cut][varname] = histo_pointer;
	  histo[sysname][chn][cut][varname]->SetDirectory(truth_file);
	}
      }
    }
  }

  for (int i=0; i<(int)CHN.size(); i++)
  {
      string chn = CHN[i];
      for(int j=0; j<(int)TRUTH_STEP_cut.size(); j++)
      {
	  string cut = TRUTH_STEP_cut[j];
	  FLAG_Truth_cut[chn][cut]["default"] = 0;
      }
  }
}

void MyxAODAnalysis :: FillHistograms(string sysname) {

  //  another shortcoming with this method is that it only fill 
  //  histograms after each cut separately, without progressive 
  //  relationship.

  MapType2_Double::iterator it;
  for(it=HistVar.begin(); it!=HistVar.end(); it++) {
    string varname=(*it).first;

    if((*it).second["NBins"] == 0) continue;

    for(int i=0; i<(int)CHN.size(); i++) {
      string chn=CHN[i];

      for(int j=0; j<(int)STEP_cut.size(); j++) {
        string cut=STEP_cut[j];
        if(int((*it).second[cut])==1 && FLAG_cut[chn][cut]["default"]) 
          histo[sysname][chn][cut][varname]->Fill(HistVar[varname]["Value"], Evt_Weight[chn][cut]);
      } 
    }
  } 

  MapType2_Double2D::iterator it2D;
  for(it2D=Hist2DVar.begin(); it2D!=Hist2DVar.end(); it2D++) {
    string varname=(*it2D).first;

    pair<double, double> nbins = (*it2D).second["NBins"];
    if(nbins.first == 0) continue;

    for(int i=0; i<(int)CHN.size(); i++) {
      string chn=CHN[i];

      for(int j=0; j<(int)STEP_cut.size(); j++) {
        string cut=STEP_cut[j];
        pair<double, double> ncuts = (*it2D).second[cut];
        if(int(ncuts.first)==1 && FLAG_cut[chn][cut]["default"]) {
          pair<double, double> values = (*it2D).second["Value"];
          histo_2D[sysname][chn][cut][varname]->Fill(values.first, values.second);
        }
      }
    }
  }


  MapType2_VDouble::iterator itv;
  for(itv=VVar.begin(); itv!=VVar.end(); itv++) {
    string varname=(*itv).first;

    if((*itv).second["NBins"][0] == 0) continue;

    for(int i=0; i<(int)CHN.size(); i++) {
      string chn=CHN[i];

      for(int j=0; j<(int)STEP_cut.size(); j++) {
        string cut=STEP_cut[j];
        if(int((*itv).second[cut].size())!=1) continue;
        if(int((*itv).second[cut][0])==1 && FLAG_cut[chn][cut]["default"]) {
          for(int k=0; k<(int)VVar[varname]["Value"].size(); k++)
            if(VVar[varname]["Weight"].size()>0)
              histo[sysname][chn][cut][varname]->Fill(VVar[varname]["Value"][k], VVar[varname]["Weight"][k]);
            else histo[sysname][chn][cut][varname]->Fill(VVar[varname]["Value"][k]);
        }
      }
    }
  }

  MapType2_V2DDouble::iterator it2Dv;
  for(it2Dv=V2DVar.begin(); it2Dv!=V2DVar.end(); it2Dv++) {
    string varname=(*it2Dv).first;

    pair<double, double> nbins = (*it2Dv).second["NBins"][0];
    if(nbins.first == 0) continue;

    for(int i=0; i<(int)CHN.size(); i++) {
      string chn=CHN[i];

      for(int j=0; j<(int)STEP_cut.size(); j++) {
        string cut=STEP_cut[j];

        if(int((*it2Dv).second[cut].size())!=1) continue;

        if(int((*it2Dv).second[cut][0].first)==1 && FLAG_cut[chn][cut]["default"]) {
          for(int k=0; k<(int)V2DVar[varname]["Value"].size(); k++) {
            pair<double, double> values = V2DVar[varname]["Value"][k];
            if(V2DVar[varname]["Weight"].size()>0)
              histo_v2D[sysname][chn][cut][varname]->Fill(values.first, values.second, V2DVar[varname]["Weight"][k].first);
            else histo_v2D[sysname][chn][cut][varname]->Fill(values.first, values.second);
          }
        }
      }
    }
  }

}

void MyxAODAnalysis :: FillTruthHistograms(string sysname) {

  MapType2_Double::iterator it;
  for(it=TruthHistVar.begin(); it!=TruthHistVar.end(); it++) {
    string varname=(*it).first;

    if((*it).second["NBins"] == 0) continue;

    for(int i=0; i<(int)CHN.size(); i++) {
      string chn=CHN[i];

      for(int j=0; j<(int)TRUTH_STEP_cut.size(); j++) {
        string cut=TRUTH_STEP_cut[j];
        if(int((*it).second[cut])==1 && FLAG_Truth_cut[chn][cut]["default"]) 
          histo[sysname][chn][cut][varname]->Fill(TruthHistVar[varname]["Value"], Evt_Weight[chn][cut]);
      } 
    }
  } 

}

void MyxAODAnalysis :: InitSetting(MapType2_Int& setmap, string setname, string settings) {
    vector<string> vsettings;
    InitStrVec(vsettings, settings, ",");
    for(int i=0; i<(int)vsettings.size(); i++) {
        vector<string> pairs;
        InitStrVec(pairs,vsettings[i],"=");
        if(pairs.size()<2) {
            cout<<"Error in setting: can not parse "<<vsettings[i]<<endl;
            exit(-1);
        }
        setmap[setname][pairs[0]]=atoi(pairs[1].c_str());
    }
}

void MyxAODAnalysis :: InitVarVec(vector<string>& map, MapType_Double& list, string cut, string settings) {
    vector<string> vsettings;
    InitStrVec(vsettings, settings, ",");

    for(int i=0; i<(int)vsettings.size(); i++) {
      string new_name = cut+"_"+vsettings[i];
      list[new_name] = atof(vsettings[i].c_str());
    }

    if(map.size()==0) {
      for(int i=0; i<(int)vsettings.size(); i++) {
        string new_name = cut+"_"+vsettings[i];
        map.push_back(new_name);
      }
    }else {
      vector<string> new_settings;
      for(int i=0; i<(int)vsettings.size(); i++) {
        string new_name = cut+"_"+vsettings[i];
        for(int j=0; j<(int)map.size(); j++) {
          string old_name=map[j];
          new_settings.push_back(old_name+";"+new_name);
        }    
      }
      map.clear();
      map=new_settings;
    }
}

