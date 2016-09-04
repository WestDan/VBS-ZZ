#ifndef MyxAODAnalysis_H
#define MyxAODAnalysis_H

#include <EventLoop/Algorithm.h>

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

namespace ana
{
    class QuickAna;
}

class MyxAODAnalysis : public EL::Algorithm, public AnalysisVar
{
public:
    time_t start, end; //!

    xAOD::TEvent * m_event; //!

    bool isMC; //! event type
    int m_eventCounter; //! total events and process monitor
    int m_fiducial; //! fiducial events number

    double m_sumOfWeights; //!

    float sumOfWeights; //!
    float sumOfWeightsSquared; //!

    string set; //! treename
    string m_FileName; //!

    ana::QuickAna * quickAna; //! quickAna ???
    std::vector<CP::SystematicSet> m_sysList; //!
    std::vector<CP::SystematicSet> pileup_sysList; //!
    CP::SystematicSet NomSet; //!

#ifndef __CINT__  // ??? Why do we need to judge __CINT__
    AsgElectronLikelihoodTool * m_likelihoodTool; //!
    OverlapRemovalTool * orTool; //!
#endif

    GoodRunsListSelectionTool * m_grlTool; //!

    // variables that don't get filled at submission time should be
    // protected from being send from the submission node to the worker
    // node. ?????
public:
    MyxAODAnalysis ();
    MyxAODAnalysis (string );
//    MyxAODAnalysis (string treename="physics");

    // these are the functions inherited from ???? Algorithm
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
    ClassDef(MyxAODAnalysis, 1);

    void ClearVariables(MapType2_Int& map);
    void ClearVariables(MapType2_Float& map);
    void ClearVariables(MapType2_VInt& map);
    void ClearVariables(MapType2_VFloat& map);

    void InitStrVec(vector<string>& out, string in, string de=",");
    void InitObjSTEP(MapType_VString& STEP, string obj, string steps);
    void InitSetting(MapType2_Int& setmap, string setname, string settings);
    void InitTreeVar(string varlist, string type);
    void AddVarIntoTree(TTree * tree);
    void CreateCountingMap();
};

#endif
