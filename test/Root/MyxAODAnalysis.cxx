#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <EventLoop/OutputStream.h>
#include <EventLoopAlgs/NTupleSvc.h>
#include <EventLoopAlgs/AlgSelect.h>

#include "MyAnalysis/MyxAODAnalysis.h"
#include "MyAnalysis/OBJ_Base.h"

#include "xAODTracking/VertexContainer.h"
#include "xAODJet/JetContainer.h"

#include <xAODMuon/Muon.h>
#include <xAODMuon/MuonContainer.h>
#include <xAODMuon/MuonAuxContainer.h>

// #include <xAODMissingET/MissingET.h>
// #include <xAODMissingET/MissingETAuxContainer.h>
// #include <xAODMissingET/MissingETContainer.h>
// #include <xAODMissingET/MissingETComponentMap.h>

#include "xAODEventInfo/EventInfo.h"

#include <QuickAna/QuickAna.h>

#include "xAODRootAccess/TStore.h"
#include "xAODCore/ShallowCopy.h"

#include "Initialize.h"
#include "Counting.h"

using namespace std;

ClassImp(MyxAODAnalysis)

MyxAODAnalysis :: MyxAODAnalysis()
{
    // variable initialization
}

MyxAODAnalysis :: MyxAODAnalysis ( string treename = "physics" )
{
    set = treename;
}

EL::StatusCode MyxAODAnalysis :: setupJob ( EL::Job& job)
{
    EL::OutputStream output1("tree_output");
    job.outputAdd (output1);

    job.useXAOD ();

    xAOD::Init( "MyxAODAnalysis" ).ignore();
    return EL::StatusCode::SUCCESS;
}

EL::StatusCode MyxAODAnalysis :: histInitialize ()
{
    return EL::StatusCode::SUCCESS;
}

EL::StatusCode MyxAODAnalysis :: fileExecute ()
{
    return EL::StatusCode::SUCCESS;
}

EL::StatusCode MyxAODAnalysis :: changeInput (bool firstFile)
{
    return EL::StatusCode::SUCCESS;
}

// EL::StatusCode MyxAODAnalysis :: initialize ()
// {
    // extract initialization on separate file
    // Initialization, after connecting to file, before event looping
// }

EL::StatusCode MyxAODAnalysis :: execute ()
{
    // execution of every events
    m_eventCounter++;
    if(m_eventCounter == 1 || m_eventCounter%1000 == 0)
    {
	// for monitoring the execution process
	cout << "event counter " << m_eventCounter << endl;
    }

    // get eventinfo
    // what's the difference between the one with the one in Initialized.h file
    const xAOD::EventInfo* eventInfo = 0;
    if( ! m_event->retrieve( eventInfo, "EventInfo").isSuccess() )
    {
	Error("execute()", "Failed to retrieve EventInfo. Exiting.");
	return EL::StatusCode::FAILURE;
    }
    uint32_t run;
    unsigned long long event;
//    float mu = eventInfo->actualInteractionsPerCrossing();
    vector<float> vw;
    float pileWeight = 1.0, mcWeight = 1.0;
    if(isMC) 
    {
	run = eventInfo->mcChannelNumber(); // ??? Channel number
	event = eventInfo->mcEventNumber();
	vw = eventInfo->mcEventWeights();  // vw ???
	if(SETTING[SETNAME[0]]["doweight"] == 1)
	{
	    mcWeight = vw.size() > 0 ? vw[0]:1.0; // vw[0] ???
	}
    }
    else
    {
	if( ! m_grlTool->passRunLB(*eventInfo))
	{
	    return EL::StatusCode::SUCCESS;
	}
	run = eventInfo->runNumber();
	event = eventInfo->eventNumber();
    }
    long double Weight = mcWeight;
    //
    // END OF EVENTINFO
    //

    VOmuon goodm;
    VOelectron goode;
    VOlepton goodl;
    VOjet goodj;

    VOmuon temp_muon;	 
    VOelectron temp_electron;	
    VOjet temp_jet;	 

    // get Vertex
    const xAOD::Vertex * pv(0);
    const xAOD::VertexContainer* vertices = 0;
    if( ! m_event->retrieve( vertices, "PrimaryVertices" ).isSuccess())
    {
	Error("execute()", "Failed to retrieve vertices. Exiting.");
	return EL::StatusCode::FAILURE; // Why failure, it is for event, if this event don't contain primaryvertex, we should continue to next event 
    }
    for(const auto* const vtx_itr : *vertices ) // why more than one primaryvertices ? What is primaryvertex ???
    {
	if(vtx_itr->vertexType() != xAOD::VxType::VertexType::PriVtx)
	    continue;
	else
	{
	    pv = vtx_itr;
	    break;
	}
    }
    if(!pv)
    {
	return EL::StatusCode::SUCCESS;
    }
    double pz0 = pv->z();

    // systematics
    for ( auto sysListItr : m_sysList )
    {
	string sysname = sysListItr.name();

	if( ! isMC && SETTING[SETNAME[0]]["dosys"] == 1 && sysname != "")
	    continue;
	if(SETTING[SETNAME[0]]["dosys"] == 0 && sysname != "")
	    continue;
	if(sysname == "")
	{
	    sysname = "NOMINAL";
	}
	if(SETTING[SETNAME[0]]["docorr"] == 0)
	{
	    sysname = "NOCORR";
	}
	// NOMINAL

	// quickAna process
	if( quickAna->applySystematicVariation (sysListItr) == CP::SystematicCode::Ok)
	{
//	    ClearWeight(Evt_Weight);
	    ClearVariables(TreeIntVar);
	    ClearVariables(TreeIntVVar);
	    ClearVariables(TreeFltVar);
	    ClearVariables(TreeFltVVar);
	    goodm.clear();  temp_muon.clear();
	    goode.clear();  temp_electron.clear();
	    goodj.clear();  temp_jet.clear();

	    quickAna->process(*m_event).ignore(); // any magic here ?

	    auto evtInfo = quickAna->eventinfo();
//	    bool passTrig = evtInfo->auxdata<bool>("passAllTrig"); // auxdata ???
	    if(SETTING[SETNAME[0]]["doweight"] == 1)
	    {
		pileWeight = evtInfo->auxdata<float>("PileupWeight");
		// pileup Weight ??? 
	    }
	    Weight = Weight*pileWeight;
	    
	    int index_muon = 0;
	    for (auto muon : *quickAna->muons())
	    {
//		typedef ElementLink<xAOD::IParticleContainer> LinkType;
//		static const char* linkName = "originalObjectLink";
//		LinkType& auxLink = muon->auxdata<LinkType> (linkName);
//		const xAOD::Muon* origMuon = dynamic_cast<const xAOD::Muon*>(*auxLink.cptr());
		OBJ_MUON muonInfo;

		muonInfo.author = muon->author();
		muonInfo.charge = muon->charge();
		muonInfo.type   = muon->muonType();
//		muonInfo.quality= muon->quality(); // non existent aux dat item

		muonInfo.L.SetPtEtaPhiM(muon->pt(), muon->eta(), muon->phi(), m_mass);

		const xAOD::TrackParticle * trkPart = muon->primaryTrackParticle(); // primary track particle ???
		muonInfo.d0 = trkPart->d0();
		muonInfo.z0 = trkPart->z0() + trkPart->vz()-pz0; // trkPart->vz() ???
		muonInfo.d0err = 0.;
		muonInfo.z0err = 0.;

		if(muonInfo.type != xAOD::Muon::MuonStandAlone)
		{
		    const xAOD::ParametersCovMatrix_t TrkCovMatrix = trkPart->definingParametersCovMatrix(); // ??? Cov Matrix
		    muonInfo.d0err = sqrt(TrkCovMatrix(0,0));
		    muonInfo.z0err = sqrt(TrkCovMatrix(1,1));
		    muonInfo.L_id.SetPtEtaPhiM(trkPart->pt(), trkPart->eta(), trkPart->phi(), m_mass); // L_id ???
		}
		else 
		{
		    muonInfo.L_id = muonInfo.L;
		}

		muonInfo.d0sig = muonInfo.d0err != 0 ? fabs(muonInfo.d0/muonInfo.d0err) : 9999.;
		muonInfo.z0sig = muonInfo.z0err != 0 ? fabs(muonInfo.z0/muonInfo.z0err) : 9999.;

		// Isolation
		float etcone20 = 0, ptcone20 = 0;
		if(muon->isolation(etcone20, xAOD::Iso::etcone20))
		{
		    muonInfo.topoetcone20 = etcone20;
		}
		if(muon->isolation(ptcone20, xAOD::Iso::ptcone20))
		{
		    muonInfo.ptcone20 = ptcone20;
		}

		muon->auxdata < char >("All") = true; // can we change the value of muon object ??? 

		if(muon->auxdata<char> ("ana_select_hzhinv_loose_ID")) muon->auxdata< char >("Tool") = true;
		if(muon->auxdata<char> ("ana_select_hzhinv_loose_CB")) muon->auxdata< char >("CB") = true;
		if(muon->auxdata<char> ("ana_select_hzhinv_loose_Eta")) muon->auxdata< char >("Eta") = true;
		if(muon->auxdata<char> ("ana_select_hzhinv_loose_Pt")) muon->auxdata< char >("Pt") = true;
		if(muon->auxdata<char> ("ana_select_hzhinv_loose_D0")) muon->auxdata< char >("D0") = true;
		if(muon->auxdata<char> ("ana_select_hzhinv_loose_Z0")) muon->auxdata< char >("Z0") = true;
		if(muon->auxdata<char> ("ana_select_hzhinv_loose_Iso")) muon->auxdata< char >("TrkIso") = true;
		if(muon->auxdata<char> ("ana_select_hzhinv_loose")) muon->auxdata< char >("OverLap") = true;

		bool  passMuon = true;
		CountMuObj(muon, STEP_obj, CNT_obj, sysname, passMuon);
		if(passMuon)
		{
		    muonInfo.trigM = muon->auxdata<bool>("HLT_mu20_iloose_L1MU15_OR_HLT_mu50_trigMatch");
		    muonInfo.ismedium = (bool) muon->auxdata<char> ("ana_select_hzhinv_medium_ID");
		    muonInfo.sf = muon->auxdata<float>("ana_weight_hzhinv_medium");
		    muonInfo.index = index_muon;
		    muon->auxdata<int>("index") = index_muon;
		    temp_muon.push_back( muonInfo );
		    index_muon++;
		}
	    }

	    for (int i=0; i<(int) temp_muon.size(); i++ ) 
	    {
		if(temp_muon[i].L.Pt()<4.e3) continue;
		DoCounting(sysname, CNT_obj, "mu", "Pt_F");

		if(fabs( temp_muon[i].L.Eta() ) > 2.8) continue;
//		DoCounting(sysname, CNT_obj, "mu", "Eta");

		if(!temp_muon[i].ismedium) continue;
		DoCounting(sysname, CNT_obj, "mu", "Medium");

		goodm.push_back(temp_muon[i]);
	    }

	    int index_ele = 0;
	    for(auto electron : *quickAna->electrons())
	    {
//		typedef ElementLink<xAOD::IParticleContainer> LinkType;
//		static const char* linkName = "originalObjectLink";
//		LinkType& auxLink = electron->auxdata<LinkType> (linkName);
//		const xAOD::Electron* origEle = dynamic_cast<const xAOD::Electron*>(*auxLink.cptr());

		OBJ_ELECTRON eleInfo;
		eleInfo.author = electron->author();
		eleInfo.charge = electron->charge();

		eleInfo.L.SetPtEtaPhiM( electron->pt(), electron->eta(), electron->phi(), e_mass );
		const xAOD::TrackParticle * trkPart = electron->trackParticle();
		eleInfo.d0 = trkPart->d0();
		eleInfo.z0 = trkPart->z0() + trkPart->vz() - pz0;
		eleInfo.d0err = 0.;
		eleInfo.z0err = 0.;
		const xAOD::ParametersCovMatrix_t TrkCovMatrix = trkPart->definingParametersCovMatrix();
		eleInfo.d0err = sqrt(TrkCovMatrix(0,0));
		eleInfo.z0err = sqrt(TrkCovMatrix(1,1));
		eleInfo.d0sig = eleInfo.d0err!=0 ? fabs(eleInfo.d0/eleInfo.d0err) : 9999.;
		eleInfo.z0sig = eleInfo.z0err!=0 ? fabs(eleInfo.z0/eleInfo.z0err) : 9999.;

		eleInfo.trketa = trkPart->eta();
		eleInfo.trkphi = trkPart->phi();
		eleInfo.trkpt = trkPart->pt();
		eleInfo.qoverp = trkPart->qOverP();
		
		const xAOD::CaloCluster * cluster = electron->caloCluster();
		eleInfo.clE = cluster->e();
		eleInfo.clpt = cluster->pt();
		eleInfo.cleta = cluster->eta();
		eleInfo.clphi = cluster->phi();
		eleInfo.L_trk.SetPtEtaPhiE( trkPart->pt(), trkPart->eta(), trkPart->phi(), cluster->e() );

		float etcone20=0, ptcone20=0;
		if(electron->isolationValue(etcone20, xAOD::Iso::etcone20))
		{
		    eleInfo.etcone20 = etcone20;
		}
		if(electron->isolationValue(ptcone20, xAOD::Iso::ptcone20))
		{
		    eleInfo.ptcone20 = ptcone20;
		}

		eleInfo.passOQ = electron->isGoodOQ(xAOD::EgammaParameters::BADCLUSELECTRON);
		electron->auxdata< char >( "All" ) = true;
		electron->auxdata< char >( "ID" ) = true;
		if(electron->auxdata< char > ("ana_select_hzhinv_loose_Pt"))   electron->auxdata< char >("Pt") = true;
		if(electron->auxdata< char > ("ana_select_hzhinv_loose_Eta"))   electron->auxdata< char >("Eta") = true;
		if(electron->auxdata< char > ("ana_select_hzhinv_loose_OQ"))   electron->auxdata< char >("ObjQ") = true;
		if(electron->auxdata< char > ("ana_select_hzhinv_loose_selectionTool"))   electron->auxdata< char >("LHID") = true;
		if(electron->auxdata< char > ("ana_select_hzhinv_loose_Z0"))   electron->auxdata< char >("Z0") = true;
		if(electron->auxdata< char > ("ana_select_hzhinv_loose_D0"))   electron->auxdata< char >("D0") = true;
		if(electron->auxdata< char > ("ana_select_hzhinv_loose_Iso"))   electron->auxdata< char >("TrkIso") = true;
		if(electron->auxdata< char > ("ana_select_hzhinv_loose"))   electron->auxdata< char >("OverLap") = true;

		bool passEle = true;
		CountEleObj(electron, STEP_obj, CNT_obj, sysname, passEle);

		if(passEle)
		{
		    if(!isMC) eleInfo.trigM = electron->auxdata<bool>("HLT_e24_lhmedium_L1EM20VH_trigMatch") || electron->auxdata<bool>("HLT_e60_lhmedium_trigMatch") || electron->auxdata<bool>("HLT_e120_lhloose_trigMatch");
		    else eleInfo.trigM = electron->auxdata<bool>("HLT_e24_lhmedium_L1EM18VH_trigMatch") || electron->auxdata<bool>("HLT_e60_lhmedium_trigMatch") || electron->auxdata<bool>("HLT_e120_lhloose_trigMatch"); 
		    eleInfo.sf = electron->auxdata<float>("ana_weight_hzhinv_medium");

		    if(eleInfo.L.Pt() > 25.e3)
		    {
			eleInfo.trigSF = electron->auxdata<double>("HLT_e60_lhmedium_TrigSF");
			eleInfo.trigEFF = electron->auxdata<double>("HLT_e60_lhmedium_TrigEff");
		    } // trigSF ??? trigEFF ???
		    else
		    {
			eleInfo.trigSF = 0.0;
			eleInfo.trigEFF = 0.0;
		    }

		    eleInfo.ismedium = (bool)electron->auxdata< char > ("ana_select_hzhinv_medium_selectionTool");
		    eleInfo.index = index_ele;
		    electron->auxdata<int>("index") = index_ele;
		    temp_electron.push_back( eleInfo );
		    index_ele++;
		}
	    }

	    for(int i=0; i<(int)temp_electron.size(); i++)
	    {
		if(temp_electron[i].L.Pt() < 4.e3) continue;
		DoCounting(sysname, CNT_obj, "ele", "Pt_F");

		if(fabs ( temp_electron[i].L.Eta() ) > 2.8 ) continue;
//		DoCounting(sysname, CNT_obj, "ele", "Eta");

		if(!temp_electron[i].ismedium) continue;
		DoCounting(sysname, CNT_obj, "ele", "Medium");

		goode.push_back(temp_electron[i]);
	    }

	    int index_jet = 0;
	    bool passJetCleaning = true;
	    bool bjettag = false;
	    float total_jetsf = 1.0;
	    int nbjetOR = 0;
	    for ( auto jet : *quickAna->jets())
	    {
//		typedef ElementLink<xAOD::IParticleContainer> LinkType;
//		static const char* linkName = "originalObjectLink";
//		LinkType& auxLink = jet->auxdata<LinkType> (linkName);
//		const xAOD::Jet* origJet = dynamic_cast<const xAOD::Jet*> (*auxLink.cptr());

		OBJ_JET jetInfo;
		jetInfo.E = jet->e();
		jetInfo.pt = jet->pt();
		jetInfo.pz = jet->pz();
		jetInfo.eta = jet->eta();
		jetInfo.phi = jet->phi();
		jetInfo.L = jet->p4();

		jetInfo.numTrk0 = jet->getAttribute< std::vector<int> >(xAOD::JetAttribute::NumTrkPt500)[0];

		jet->auxdata< char >("All") = true;

		if( jet->pt() > 15.e3 ) jet->auxdata< char >("Pt") = true;
		if( fabs(jet->eta()) < 9 ) jet->auxdata< char >("Eta") = true;
		if(jet->auxdata<char>("ana_select_jvt")) jet->auxdata< char >("JVT") = true;
		if(jet->auxdata<char>("ana_select")) jet->auxdata< char >("OverLap") = true;

		bool passjet = true;
		CountJetObj(jet, STEP_obj, CNT_obj, sysname, passjet);
		if(passjet)
		{
		    if(jet->auxdecor<char>("bjet")) nbjetOR++;
		    bjettag = bjettag || (bool)jet->auxdecor<char>("bjet");

		    if(jet->auxdata<char>("ana_select_cleaning_tool"))
		    {
			DoCounting(sysname, CNT_obj, "jet", "Clean");
			jetInfo.index = index_jet;
			jetInfo.sf = jet->auxdata<float>("ana_weight");
			total_jetsf = total_jetsf*jetInfo.sf;
			jet->auxdata<int>("index") = index_jet;
			goodj.push_back( jetInfo );
			index_jet++;
		    }
		    else passJetCleaning = false;
		}
	    }

	    /*
	    // missing Et
	    xAOD::MissingET * met_rebuild = quickAna->met();
	    double rmet = met_rebuild->met();
	    double rmpx = met_rebuild->mpx();
	    double rmpy = met_rebuild->mpy();
	    double rmphi= met_rebuild->phi();

	    OBJ_MET met_obj;
	    met_obj.L_RefFinal.SetPxPyPzE(rmpx, rmpy, 0, rmet);
	    met_obj.L = met_obj.L_RefFinal;
	    met_obj.met_rebuild = rmet;
	    goodmet.push_back(met_obj);
	    */

	    int nmuon = goodm.size();
	    int nele  = goode.size();
	    int njet  = goodj.size();
	    if ( (nmuon + nele) < 4 ) return EL::StatusCode::SUCCESS;
	    if ( njet < 2 ) return EL::StatusCode::SUCCESS;
	    m_fiducial++;

	    // extract leptons
	    OBJ_LEPTON lepton;
	    for(int i=0; i<(int)goodm.size(); i++)
	    {
		lepton.pt = goodm[i].L.Pt();
		lepton.eta = goodm[i].L.Eta();
		lepton.phi = goodm[i].L.Phi();
		lepton.E = goodm[i].L.E();
//		lepton.type = goodm[i].type();
//		lepton.quality = goodm[i].quality();
		if(goodm[i].charge == 1)
		{
		    lepton.ID = 13;
		}
		else
		{
		    lepton.ID = -13;
		}
		goodl.push_back(lepton);
	    }
	    for(int i=0; i<(int)goode.size(); i++)
	    {
		lepton.pt = goode[i].L.Pt();
		lepton.eta = goode[i].L.Eta();
		lepton.phi = goode[i].L.Phi();
		lepton.E = goode[i].L.E();
//		lepton.type = goode[i].type();
//		lepton.quality = goode[i].quality();
		if(goode[i].charge == 1)
		{
		    lepton.ID = 11;
		}
		else
		{
		    lepton.ID = -11;
		}
		goodl.push_back(lepton);
	    }

	    /*
	    // find the first 4 largest pt leptons
	    int lep[4] = {0, 0, 0, 0};
	    for(int i=0; i<(int)goodl.size(); i++)
	    {
		double pt = goodl[i].pt;
		if (pt > goodl[lep[0]].pt)
		{
		    lep[3] = lep[2];
		    lep[2] = lep[1];
		    lep[1] = lep[0];
		    lep[0] = i;
		}
		else if ( pt > goodl[lep[1]].pt)
		{
		    lep[3] = lep[2];
		    lep[2] = lep[1];
		    lep[1] = i;
		}
		else if ( pt > goodl[lep[2]].pt)
		{
		    lep[3] = lep[2];
		    lep[2] =i;
		}
		else if ( pt > goodl[lep[3]].pt)
		{
		    lep[3] = i;
		}
	    }
	     */

	    // assignment of leptons properties 
	    vector<float> l_pt, l_eta, l_phi, l_e;
	    vector<int> l_id;
	    for (int i=0; i<(int)goodl.size(); i++)
	    {
		l_pt.push_back( goodl[i].pt );
		l_eta.push_back( goodl[i].eta );
		l_phi.push_back( goodl[i].phi );
		l_e.push_back( goodl[i].E );
//		type_l[i] = goodl[i].type;
//		quality_l[i] = goodl[i].quality;
		l_id.push_back( goodl[i].ID );
	    }

	    /*
	    // find the first 2 largest pt jets
	    int jet[2] = {0 ,0};
	    for (int i=0; i<(int)goodj.size(); i++)
	    {
		double pt = goodj[i].pt;
		if( pt > goodj[jet[0]].pt )
		{
		    jet[1] = jet[0];
		    jet[0] = i;
		}
		else if ( pt > goodj[jet[1]].pt )
		{
		    jet[1] = i;
		}
	    }
	    */
	    // assignment of jet properties
	    vector<float> j_pt, j_eta, j_phi, j_e;
	    for (int i=0; i<(int)goodj.size(); i++)
	    {
		j_pt.push_back(  goodj[i].pt );
		j_eta.push_back( goodj[i].eta );
		j_phi.push_back( goodj[i].phi );
		j_e.push_back(   goodj[i].E );
		// any other jet properties that need to be recorded ?
	    }

	    TreeFltVVar["l_pt"]["Value"] = l_pt;
	    TreeFltVVar["l_eta"]["Value"] = l_eta;
	    TreeFltVVar["l_phi"]["Value"] = l_phi;
	    TreeFltVVar["l_e"]["Value"] = l_e;
	    TreeIntVVar["l_id"]["Value"] = l_id;
	    TreeFltVVar["j_pt"]["Value"] = j_pt;
	    TreeFltVVar["j_eta"]["Value"] = j_eta;
	    TreeFltVVar["j_phi"]["Value"] = j_phi;
	    TreeFltVVar["j_e"]["Value"] = j_e;

	    TreeStrVar["filename"]["Value"] = m_FileName;
	    TreeIntVar["signal_bkg"]["Value"] = SB::_signal;
	    TreeIntVar["run"]["Value"] = run;
	    TreeUloVar["event"]["Value"] = event;
	    TreeIntVar["njet"]["Value"] = njet;
	    TreeIntVar["nmuon"]["Value"] = nmuon;
	    TreeIntVar["nele"]["Value"] = nele;
	    TreeFltVar["nbjetOR"]["Value"] = nbjetOR;
//	    TreeFltVar["met"]["Value"] = rmet*0.001;
//	    TreeFltVar["met_x"]["Value"] = rmpx*0.001;
//	    TreeFltVar["met_y"]["Value"] = rmpy*0.001;
	    TreeFltVar["pileupweight"]["Value"] = pileWeight;
	    TreeFltVar["mcweight"]["Value"] = mcWeight;
//	    TreeFltVar["lep_sf"]["Value"] = lep_sf;
//	    TreeFltVar["trigsf"]["Value"] = trigger_sf;
	    Tree[sysname]->Fill();

	}
    }
    return EL::StatusCode::SUCCESS;
}

EL::StatusCode MyxAODAnalysis :: postExecute ()
{
    return EL::StatusCode::SUCCESS;
}

EL::StatusCode MyxAODAnalysis :: finalize ()
{
    return EL::StatusCode::SUCCESS;
}

EL::StatusCode MyxAODAnalysis :: histFinalize ()
{
    return EL::StatusCode::SUCCESS;
}
