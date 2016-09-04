#include <iostream>
#include <fstream>
#include "TSystem.h"
using namespace std;

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>

#include "xAODMuon/Muon.h"

#include <QuickAna/QuickAna.h>

#include "MyAnalysis/MyxAODAnalysis.h"
#include "MyAnalysis/AnalysisVar.h"

#include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"

#include "AssociationUtils/OverlapRemovalTool.h"
#include "PileupReweighting/PileupReweightingTool.h"

#include "xAODCutFlow/CutBookkeeper.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"

EL::StatusCode MyxAODAnalysis :: initialize ()
{
    time(&start);

    SETNAME.push_back("physics");

    InitSetting(SETTING,"physics","docorr=1,dosys=0,doweight=1,dovbf=1,doopt=0");
    // step cut for obj 
    InitObjSTEP(STEP_obj,"ele","All,ID,LHID,Pt,Eta,ObjQ,Z0,D0,TrkIso,OverLap,Pt_F,Medium");
    InitObjSTEP(STEP_obj,"mu","All,Tool,CB,Eta,Pt,Z0,D0,TrkIso,OverLap,Pt_F,Medium");
    InitObjSTEP(STEP_obj,"jet","All,Pt,Eta,JVT,OverLap,Clean");

    m_event = wk()->xaodEvent();
    Info("initialize()", "Number of event = %lli", m_event->getEntries() );
    m_FileName = wk()->inputFile()->GetName();
    cout << endl << endl << endl
	 << "File name: " << m_FileName << endl;

    const xAOD::EventInfo * eventInfo = 0; 
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() )
    {
	Error("initialize()", "Failed to retrieve EventInfo, Exiting. ");
	return EL::StatusCode::FAILURE;
    }
    isMC = eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION );

    uint64_t nEventsProcessed = 0;
    sumOfWeights = 0.;
    sumOfWeightsSquared = 0.;

    if(isMC)
    {
	TTree * MetaData = dynamic_cast<TTree*>(wk()->inputFile()->Get("MetaData"));
	if(! MetaData)
	{
	    Error("fileExecute()", "MetaData not found! Exiting.");
	    return EL::StatusCode::FAILURE;
	}

	MetaData->LoadTree(0); // ???
	bool m_isDerivation = ! MetaData->GetBranch("StreamAOD"); // ??? !
	if(m_isDerivation)
	{
	    cout << "is derivation" << endl;
	    const xAOD::CutBookkeeperContainer* incompleteCBC = nullptr;
	    if(!m_event->retrieveMetaInput(incompleteCBC, "IncompleteCutBookkeepers").isSuccess())
	    {
		Error("initializeEvent()", "Failed to retrieve IncompleteCutBookkeepers from MetaData! Exiting.");
		return EL::StatusCode::FAILURE;
	    }
	    if ( incompleteCBC->size() != 0 )
	    {
		Error("initializeEvent()","Found incomplete Bookkeepers! Check file for corruption.");
		return EL::StatusCode::FAILURE;
	    }

	    const xAOD::CutBookkeeperContainer* completeCBC = 0;
	    if(!m_event->retrieveMetaInput(completeCBC, "CutBookkeepers").isSuccess())
	    {
		Error("initializeEvent()", "Failed to retrieve CutBookkeepes from MetaData! Exiting.");
		return EL::StatusCode::FAILURE;
	    } // what's CutBookkeeper

	    // find minCycle, what's this ???
	    int minCycle = 10000;
	    for ( auto cbk : *completeCBC )
	    {
		if ( ! cbk->name().empty() && minCycle > cbk->cycle() )
		{
		    minCycle = cbk->cycle();
		}
	    }

	    // find "AllExecutedEvents" from 
	    const xAOD::CutBookkeeper* allEventsCBK = 0;
	    for ( auto cbk : *completeCBC )
	    {
		if ( minCycle == cbk->cycle() && cbk->name() == "AllExecutedEvents" )
		{
		    allEventsCBK = cbk;
		    break;
		}
	    }

	    if(allEventsCBK)  // find 
	    {
		nEventsProcessed = allEventsCBK->nAcceptedEvents();
		sumOfWeights = allEventsCBK->sumOfEventWeights();
		sumOfWeightsSquared = allEventsCBK->sumOfEventWeightsSquared();
	    }
	} // what if it is not derivation
    } // what if it is not MC
    cout << "nEventsProcessed: " << nEventsProcessed << " sumOfWeights: " << sumOfWeights << " sumOfEventWeightsSquared: " << sumOfWeightsSquared << endl;
//    TFile *f = new TFile("SumOfWeights.root","RECREATE");
//    TH1F *h1 = new TH1F("sumOfWeights", "sumOfWeights", 0,10,10);
//    h1->Fill(5, sumOfWeights);
//    h1->Write();
//    f->Close();

    m_grlTool = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
    const char* GRLFilePath = "/afs/cern.ch/user/w/weibin/Test/ExtraTools/grl/data15_13TeV.periodAllYear_DetStatus-v73-pro19-08_DQDefects-00-01-02_PHYS_StandardGRL_All_Good_25ns.xml";  // !!!! important, fill with proper grlfilepath.
    const char* fullGRLFilePath = gSystem->ExpandPathName(GRLFilePath);
    std::vector<std::string> vecStringGRL;
    vecStringGRL.push_back(fullGRLFilePath);
    m_grlTool->setProperty("GoodRunListVec", vecStringGRL);
    m_grlTool->setProperty("PassThrough", false);
    m_grlTool->initialize(); // Is there any selection here ???

    std::vector<std::string> prwFiles;
    string PileConfigFile = "/afs/cern.ch/user/w/weibin/Test/ExtraTools/prw/MC15b.prw.ttbar.root"; // !!! important, fill it with proper value
    prwFiles.push_back(PileConfigFile);
    std::vector<std::string> lumicalcFiles;
    string LumiCalcFile3 = "/afs/cern.ch/user/w/weibin/Test/ExtraTools/prw/ilumicalc_histograms_None_276262-284484.root"; // !!!! important, fill it with proper value 
    lumicalcFiles.push_back(LumiCalcFile3);

    quickAna = new ana::QuickAna("quickana");
    if(!isMC)
    {
	quickAna->isDataFlag = "true";
    }
    quickAna->eventinfoDef = "default";
    if(!isMC)
    {
	quickAna->triggerDef = "HLT_mu20_iloose_L1MU15_OR_HLT_mu50 HLT_e24_lhmedium_L1EM20VH HLT_e60_lhmedium HLT_e120_lhloose";
    }
    else
    {
	quickAna->triggerDef = "HLT_mu20_iloose_L1MU15_OR_HLT_mu50 HLT_e24_lhmedium_L1EM18VH HLT_e60_lhmedium HLT_e120_lhloose";
    }
    quickAna->muonDef="hzhinv_loose hzhinv_medium";
    quickAna->electronDef="hzhinv_loose hzhinv_medium";
    quickAna->jetKine="pt > 20e3 && eta > -4.5 && eta < 4.5";
    quickAna->jetDef="antikt04_HZZ";
    quickAna->tauDef="none";
    quickAna->photonDef="none";
    quickAna->metDef="metZHinv";
    quickAna->orDef="default";
    quickAna->muDataFiles=lumicalcFiles;
    quickAna->muMcFiles=prwFiles;
    quickAna->initialize().ignore();

    m_sysList = CP::make_systematics_vector (quickAna->recommendedSystematics());
    CP::SystematicCode::enableFailure(); // ???? enableFailure()

    SYSNAME.push_back("NOCORR");
    for ( auto sysListItr : m_sysList )
    {
	cout << "sys name is : " << sysListItr.name() << endl;
	if(sysListItr.name() == "" )
	{
	    SYSNAME.push_back("NOMINAL");
        }
	else
	{
	    SYSNAME.push_back(sysListItr.name());
	}
    }
    
    TFile * file1 = wk()->getOutputFile("tree_output");
    if(SETTING["physics"]["dobkgwz"] || SETTING["physics"]["dobkgz"])
    {
	// ??? dobkgwz or dobkgz
	tree_bkg = new TTree("tree_bkg", "tree_bkg");
	tree_bkg->SetDirectory (file1);
	AddVarIntoTree(tree_bkg);
    }
    else
    {
	for(int i=0; i<(int)SYSNAME.size(); i++)
	{
	    string sys=SYSNAME[i];
	    if(sys != "NOCORR" && (!SETTING["physics"]["docorr"] ))
	    {
		continue;
	    }
	    if(sys != "NOMINAL" && SETTING["physics"]["docorr"] && (! SETTING["physics"]["dosys"] ))
	    {
		continue;
	    }
	    if(sys == "NOCORR" && SETTING["physics"]["docorr"] && SETTING["physics"]["dosys"] )
	    {
		continue;
	    }
            // NOMINAL is selected
	    string tree_name = "tree_" +sys;
	    Tree[sys] = new TTree(tree_name.c_str(), "output tree");
	    Tree[sys]->SetDirectory (file1);
	    AddVarIntoTree(Tree[sys]);
	}
    }

    m_eventCounter = 0;
    m_fiducial = 0;

    m_sumOfWeights = 0.;

    CreateCountingMap();

    return EL::StatusCode::SUCCESS;
}

void MyxAODAnalysis :: ClearVariables(MapType2_Int& map)
{
    MapType2_Int::iterator it;
    for(it=map.begin(); it!=map.end(); it++)
    {
	string varname = (*it).first;
	map[varname]["Value"] = -9999;
    }
}

void MyxAODAnalysis :: ClearVariables(MapType2_Float& map)
{
    MapType2_Float::iterator it;
    for(it=map.begin(); it!=map.end(); it++)
    {
	string varname = (*it).first;
	map[varname]["Value"] = -9999.;
    }
}

void MyxAODAnalysis :: ClearVariables(MapType2_VInt& map)
{
    MapType2_VInt::iterator it;
    for(it=map.begin(); it!=map.end(); it++)
    {
	string varname = (*it).first;
	map[varname]["Value"].clear();
    }
}
void MyxAODAnalysis :: ClearVariables(MapType2_VFloat& map)
{
    MapType2_VFloat::iterator it;
    for(it=map.begin(); it!=map.end(); it++)
    {
	string varname = (*it).first;
	map[varname]["Value"].clear();
    }
}
void MyxAODAnalysis :: InitStrVec(vector<string>& out, string in, string de)
{
    int pos=0, pos_pre=0;
    while(true)
    {
	pos = in.find(de,pos_pre);
	if(pos == -1)
	{
	    out.push_back(in.substr(pos_pre, in.size()-pos_pre));
	    break;
	}
	else
	{
	    out.push_back(in.substr(pos_pre, pos-pos_pre));
	    pos_pre = pos + 1;
	}
    }
}

void MyxAODAnalysis :: InitObjSTEP(MapType_VString& STEP, string obj, string steps)
{
    vector<string> str, objstr;
    InitStrVec(str, steps, ",");
    STEP[obj] = str;
}

void MyxAODAnalysis :: InitSetting(MapType2_Int& setmap, string setname, string settings)
{
    vector<string> vsettings;
    InitStrVec(vsettings, settings, ",");
    for(int i=0; i<(int)vsettings.size(); i++)
    {
	vector<string> pairs;
	InitStrVec(pairs, vsettings[i], "=");
	if(pairs.size() < 2)
	{
	    cout << "Error in setting: can not parse " << vsettings[i] << endl;
	    exit(-1);
	}
	setmap[setname][pairs[0]] = atoi(pairs[1].c_str());
    }
}

void MyxAODAnalysis :: InitTreeVar(string varlist, string type)
{
    vector<string> variables;
    InitStrVec(variables, varlist, ",");

    if(type == "I")
    {
	for(int i=0; i<(int)variables.size(); i++)
	{
	    TreeIntVar[variables[i]]["Value"] = -9999;
	}
    }
    if(type == "U")
    {
	for(int i=0; i<(int)variables.size(); i++)
	{
	    TreeUloVar[variables[i]]["Value"] = 0;
	}
    }
    if(type == "F")
    {
	for(int i=0; i<(int)variables.size(); i++)
	{
	    TreeFltVar[variables[i]]["Value"] = -9999.0;
	}
    }
    if(type == "S")
    {
	for(int i=0; i<(int)variables.size(); i++)
	{
	    TreeStrVar[variables[i]]["Value"] = " ";
	}
    }
    if(type == "vector_F")
    {
	for(int i=0; i<(int)variables.size(); i++)
	{
	    TreeFltVVar[variables[i]]["Value"].clear();
	}
    }
    if(type == "vector_I")
    {
	for(int i=0; i<(int)variables.size(); i++)
	{
	    TreeIntVVar[variables[i]]["Value"].clear();
	}
    }
}

void MyxAODAnalysis :: AddVarIntoTree(TTree * tree)
{
    char buf[1000];
    getcwd(buf, sizeof(buf)); // get current working dir
    setenv("LOCAL", buf, 1);  // set env
    FILE* fp;
    char result_buf[1000];
    fp = popen("find $LOCAL -name MiniTree.txt", "r");
    fgets(result_buf, sizeof(result_buf), fp);

    string varfile(result_buf);
    size_t len = varfile.size();
    varfile.erase(len-1);
    ifstream file;
    file.open(varfile.c_str(), ios::out); // ???? out

    if (file.is_open())
    {
	char line[256];
	while ( ! file.eof() )
	{
	    string varname, type;

	    file.getline (line, 100);
	    string sline(line);

	    if(sline.find("Int_t") != string::npos)
	    {
		type = "I";
		varname = sline.substr(9);
		type = varname + "/" + type;
		InitTreeVar(varname, "I");
		tree->Branch(varname.c_str(), &TreeIntVar[varname]["Value"], type.c_str());
	    }
	    if(sline.find("Ulo_t") != string::npos)
	    {
		type = "U";
		varname = sline.substr(9);
		type = varname + "/" +type;
		InitTreeVar(varname, "U");
		tree->Branch(varname.c_str(), &TreeUloVar[varname]["Value"]);
	    }
	    if(sline.find("Float_t") != string::npos)
	    {
		type = "F";
		varname = sline.substr(9);
		type = varname + "/" +type;
		InitTreeVar(varname, "F");
		tree->Branch(varname.c_str(), &TreeFltVar[varname]["Value"], type.c_str());
	    }
	    if(sline.find("Str") != string::npos)
	    {
		type = "S";
		varname = sline.substr(9);
		type = varname + "/" +type;
		InitTreeVar(varname, "S");
		tree->Branch(varname.c_str(), &TreeStrVar[varname]["Value"]);
	    }
	    if(sline.find("Vector_I") != string::npos)
	    {
		varname = sline.substr(9);
		InitTreeVar(varname, "Vector_I");
		tree->Branch(varname.c_str(), &TreeIntVVar[varname]["Value"]);
	    }
	    if(sline.find("Vector_F") != string::npos)
	    {
		varname = sline.substr(9);
		InitTreeVar(varname, "Vector_F");
		tree->Branch(varname.c_str(), &TreeFltVVar[varname]["Value"]);
	    }
	}
    }
}

void MyxAODAnalysis :: CreateCountingMap()
{
    MapType_VString::iterator it;
    for(it=STEP_obj.begin(); it!=STEP_obj.end(); it++)
    {
	string obj = (*it).first;
	for(int i=0; i<(int)(*it).second.size(); i++)
	{
	    string cut = STEP_obj[obj][i];
	    COUNT ini = {0., 0., 0.};

	    for(int j=0; j<(int)SYSNAME.size(); j++)
	    {
		CNT_obj[SYSNAME[j]][obj][cut] = ini;
	    }
	}

    }
}
