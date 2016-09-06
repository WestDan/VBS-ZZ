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
if(m_eventCounter == 1796 || m_eventCounter ==5613)
{
    cout << m_eventCounter << "beginning of execution" << endl;
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

if(m_eventCounter == 1796 || m_eventCounter ==5613)
{
    cout << m_eventCounter << "before quickAna loop" << endl;
}
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
		{    muonInfo.topoetcone20 = etcone20; 	}
		if(muon->isolation(ptcone20, xAOD::Iso::ptcone20))
		{    muonInfo.ptcone20 = ptcone20; }

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
// Is this needed ??? Is there selection in auxdata<char>("Pt")
		DoCounting(sysname, CNT_obj, "mu", "Pt_F");

		if(fabs( temp_muon[i].L.Eta() ) > 2.8) continue; 
// Is this needed ??? Is there selection in auxdata<char>("Eta")
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

		    if(eleInfo.L.Pt() > 25.e3) // Why only apply for electrons with pt larger than 25GeV
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
	    if( ! ( (nmuon == 4 && nele == 0) || ( nmuon == 2 && nele == 2) || ( nmuon == 0 && nele == 4) ) ) 
		return EL::StatusCode::SUCCESS;
//	    if ( (nmuon + nele) != 4 ) return EL::StatusCode::SUCCESS; // at least 4 leptons needed.
	    m_fiducial++; // used for selection efficiency

	    // extract leptons
	    vector<Pair> ll_temp, ll_good;
	    ll_temp.clear();
	    ll_good.clear();
if(m_eventCounter == 1796 || m_eventCounter ==5613)
{
    cout << "Muon: " << endl;
    for(int i=0; i<nmuon; i++)
    {
	cout << i << "\t" << goodm[i].L.Pt() << "\t" << goodm[i].L.Eta() << "\t" << goodm[i].charge << endl;
    }
    cout << "electron: " << endl;
    for(int i=0; i<nele; i++)
    {
	cout << i << "\t" << goode[i].L.Pt() << "\t" << goode[i].L.Eta() << "\t" << goode[i].charge << endl;
    }
}
	    for(int i=0; i<nmuon; i++)
	    {
		for(int j=i+1; j<nmuon; j++)
		{
		    // opposite charged
		    if( (goodm[i].charge + goodm[j].charge) == 0)
		    { 
			Pair temp;
			temp.flavor = 1;
			temp.index = {i, j};
			temp.Z = goodm[i].L + goodm[j].L;
			temp.mass = temp.Z.M();
			temp.sf = goodm[i].sf * goodm[j].sf;
			temp.trigSF = evtInfo->auxdata<double>("HLT_mu20_iloose_L1MU15_OR_HLT_mu50_Mu_TrigSF");
			ll_temp.push_back(temp);
		    }
		}
	    }
            for(int i=0; i<nele; i++)
	    {
		for(int j=i+1; j<nele; j++)
		{
		    if( (goode[i].charge + goode[j].charge) == 0)
		    {
			Pair temp;
			temp.flavor = 0;
			temp.index  = {i, j};
			temp.Z      = goode[i].L + goode[j].L;
			temp.mass   = temp.Z.M();
			temp.sf     = goode[i].sf*goode[j].sf;
			float ineff_data = 1.0;
			float ineff_MC   = 1.0;
			if(goode[i].L.Pt() > 25.e3)
			{
			    ineff_data *= (1-goode[i].trigSF*goode[i].trigEFF);
			    ineff_MC   *= (1-goode[i].trigEFF);
			}
			if(goode[j].L.Pt() > 25.e3)
			{
			    ineff_data *= (1-goode[j].trigSF*goode[j].trigEFF);
			    ineff_MC   *= (1-goode[j].trigEFF);
			}
			if ( ineff_MC == 1 )
			{   temp.trigSF = 0.0;  }
			else
			{   temp.trigSF = (1-ineff_data)/(1-ineff_MC);   }
			ll_temp.push_back(temp);
		    }
		}
	    }
	    int z_index1 = -1, z_index2 = -1;
	    float diff = 9999.e9; 
	    int z_index = -1; // for case of only one suitable pair;
	    float dis = 9999.e9;
	    for(int i=0; i<(int)ll_temp.size(); i++)
	    {
		int pair1_1 = ll_temp[i].index.first;
		int pair1_2 = ll_temp[i].index.second;
		for(int j=i+1; j<(int)ll_temp.size(); j++)
		{
		    if(ll_temp[j].flavor == ll_temp[i].flavor )
		    {
			int pair2_1 = ll_temp[j].index.first;
			int pair2_2 = ll_temp[j].index.second;
			if( pair2_1 == pair1_1 || pair2_1 == pair1_2 || pair2_2 == pair1_1 || pair2_2 == pair1_2) continue;
		    }
		    float diff_temp = fabs(ll_temp[i].mass - ZMass) + fabs(ll_temp[j].mass - ZMass);
		    if( diff_temp < diff)
		    {
			diff = diff_temp;
			z_index1 = i;
			z_index2 = j;
		    }
		}
		float dis_temp = fabs(ll_temp[i].mass - ZMass);
		if(dis_temp < dis)
		{
		    dis = dis_temp;
		    z_index = i;
		}
	    }
	    if(z_index1 != -1 && z_index2 != -1) // at least 2 pairs
	    {
		if( fabs(ll_temp[z_index1].mass - ZMass) < fabs(ll_temp[z_index2].mass - ZMass))
		{
		    ll_good.push_back(ll_temp[z_index1]);
		    ll_good.push_back(ll_temp[z_index2]);
		}
		else
		{
		    ll_good.push_back(ll_temp[z_index2]);
		    ll_good.push_back(ll_temp[z_index1]);
		}
	    }
	    else if(z_index != -1) // only one good pair;
	    {
		ll_good.push_back(ll_temp[z_index]);
	    }

	    float l1_pt = -9999., l1_eta = -9999., l1_phi = -9999., l1_e = -9999.;
	    float l2_pt = -9999., l2_eta = -9999., l2_phi = -9999., l2_e = -9999.;
	    float l3_pt = -9999., l3_eta = -9999., l3_phi = -9999., l3_e = -9999.;
	    float l4_pt = -9999., l4_eta = -9999., l4_phi = -9999., l4_e = -9999.;
	    float z1_pt = -9999., z1_eta = -9999., z1_phi = -9999., z1_e = -9999., z1_m = -9999.;
	    float z2_pt = -9999., z2_eta = -9999., z2_phi = -9999., z2_e = -9999., z2_m = -9999.;
	    float zz_pt = -9999., zz_eta = -9999., zz_phi = -9999., zz_e = -9999., zz_m = -9999.;
	    float zz_rap = -9999.;
	    float delta_R_ll_1 = -9999., delta_R_ll_2 = -9999.;
	    float delta_R_zz = -9999., delta_R_jj = -9999.;
	    float delta_R_zz_jj = -9999.;
	    float lep_sf = -9999.;
	    float trigger_sf = -9999.;
	    string channel = "None";
	    if(ll_good.size() == 1)
	    {
		z1_pt  = ll_good[0].Z.Pt();
		z1_eta = ll_good[0].Z.Eta();
		z1_phi = ll_good[0].Z.Phi();
		z1_e   = ll_good[0].Z.E();
		z1_m   = ll_good[0].mass;
		lep_sf = ll_good[0].sf;
		trigger_sf = ll_good[0].trigSF;
	    }
	    else if(ll_good.size() == 2)
	    {
		{
		    z1_pt  = ll_good[0].Z.Pt();
		    z1_eta = ll_good[0].Z.Eta();
		    z1_phi = ll_good[0].Z.Phi();
		    z1_e   = ll_good[0].Z.E();
		    z1_m   = ll_good[0].mass;
		    z2_pt  = ll_good[1].Z.Pt();
		    z2_eta = ll_good[1].Z.Eta();
		    z2_phi = ll_good[1].Z.Phi();
		    z2_e   = ll_good[1].Z.E();
		    z2_m   = ll_good[1].mass;
		    lep_sf = ll_good[0].sf * ll_good[1].sf;
		}
		TLorentzVector zz = ll_good[0].Z + ll_good[1].Z;
		zz_pt  = zz.Pt();
		zz_eta = zz.Eta();
		zz_phi = zz.Phi();
		zz_e   = zz.E();
		zz_m   = zz.M();
		zz_rap = 0.5*log( (zz.E() + zz.Pz()) / ( zz.E() - zz.Pz()));
		delta_R_zz = sqrt( pow((z1_eta - z2_eta), 2) + pow((z1_phi - z2_phi), 2));
		// channel
		if(ll_good[0].flavor == 0 && ll_good[1].flavor == 0)
		{   
		    channel = "eeee";	
		    float ineff_data = 1.0, ineff_MC = 1.0;
		    for(int i=0; i<2; i++)
		    {
			if(goode[ll_good[i].index.first].L.Pt()>25.e3)
			{
			    ineff_data *= (1-goode[ll_good[i].index.first].trigSF*goode[ll_good[i].index.first].trigEFF);
			    ineff_MC   *= (1-goode[ll_good[i].index.first].trigEFF);
			}
			if(goode[ll_good[i].index.second].L.Pt()>25.e3)
			{
			    ineff_data *= (1-goode[ll_good[i].index.second].trigSF*goode[ll_good[i].index.second].trigEFF);
			    ineff_MC   *= (1-goode[ll_good[i].index.second].trigEFF);
			}
		    }
		    if(ineff_MC == 1) trigger_sf = 0.0;
		    else trigger_sf = (1-ineff_data)/(1-ineff_MC);
		}
		else if(ll_good[0].flavor == 1 && ll_good[1].flavor==1)
		{   
		    channel = "mmmm";	
		    trigger_sf = evtInfo->auxdata<double>("HLT_mu20_iloose_L1MU15_OR_HLT_mu50_Mu_TrigSF");
		}
		else
		{  
		    channel = "eemm";	
		    trigger_sf = ll_good[0].trigSF * ll_good[1].trigSF;
		}

		TLorentzVector leptons[4];
		if(ll_good[0].flavor == 0) // electron
		{
		    leptons[0] = goode[ll_good[0].index.first].L;
		    leptons[1] = goode[ll_good[0].index.second].L;
		}
		else if(ll_good[0].flavor == 1)
		{
		    leptons[0] = goodm[ll_good[0].index.first].L;
		    leptons[1] = goodm[ll_good[0].index.second].L;
		}
		delta_R_ll_1 = sqrt(pow((leptons[0].Eta() - leptons[1].Eta()), 2) + pow((leptons[0].Phi() - leptons[1].Phi()), 2));
		if(ll_good[1].flavor == 0)
		{
		    leptons[2] = goode[ll_good[1].index.first].L;
		    leptons[3] = goode[ll_good[1].index.second].L;
		}
		else if(ll_good[1].flavor == 1)
		{
		    leptons[2] = goodm[ll_good[1].index.first].L;
		    leptons[3] = goodm[ll_good[1].index.second].L;
		}
		delta_R_ll_2 = sqrt(pow((leptons[2].Eta() - leptons[3].Eta()), 2) + pow((leptons[2].Phi() - leptons[3].Phi()), 2));
		for(int i=0; i<4; i++)
		{
		    for(int j=i+1; j<4; j++)
		    {
			if(leptons[j].Pt() > leptons[i].Pt())
			{
			    TLorentzVector temp = leptons[i];
			    leptons[i] = leptons[j];
			    leptons[j] = temp;
			}
		    }
		}
		// extraction of pt, eta, phi, e, value
		{
		    l1_pt  = leptons[0].Pt();
		    l1_eta = leptons[0].Eta();
		    l1_phi = leptons[0].Phi();
		    l1_e   = leptons[0].E();
		    l2_pt  = leptons[1].Pt();
		    l2_eta = leptons[1].Eta();
		    l2_phi = leptons[1].Phi();
		    l2_e   = leptons[1].E();
		    l3_pt  = leptons[2].Pt();
		    l3_eta = leptons[2].Eta();
		    l3_phi = leptons[2].Phi();
		    l3_e   = leptons[2].E();
		    l4_pt  = leptons[3].Pt();
		    l4_eta = leptons[3].Eta();
		    l4_phi = leptons[3].Phi();
		    l4_e   = leptons[3].E();
		}
	    }
	    // assignment of jet properties
	    float leadingj_pt  = -9999.;
	    float leadingj_eta = -9999.;
	    float leadingj_phi = -9999.;
	    float leadingj_e   = -9999.;
	    float leadingj_rap = -9999.;
	    float subleadingj_pt  = -9999.;
	    float subleadingj_eta = -9999.;
	    float subleadingj_phi = -9999.;
	    float subleadingj_e   = -9999.;
	    float subleadingj_rap = -9999.;
	    float jj_pt = -9999., jj_eta = -9999., jj_phi = -9999., jj_e = -9999., jj_m = -9999.;
	    float delta_eta_jj = -9999.;
	    float centrality = -9999.;
	    pair<int, int> jj_index = {-1, -1};
	    for (int i=0; i<(int)goodj.size(); i++)
	    {
		if(goodj[i].L.Pt() > leadingj_pt)
		{
		    subleadingj_pt  = leadingj_pt;
		    subleadingj_eta = leadingj_eta;
		    subleadingj_phi = leadingj_phi;
		    subleadingj_e   = leadingj_e;
		    subleadingj_rap = leadingj_rap;
		    jj_index.second = jj_index.first;
		    leadingj_pt	    = goodj[i].L.Pt();
		    leadingj_eta    = goodj[i].L.Eta();
		    leadingj_phi    = goodj[i].L.Phi();
		    leadingj_e      = goodj[i].L.E();
		    jj_index.first  = i;
		    leadingj_rap = 0.5*log( (goodj[i].L.E() + goodj[i].L.Pz()) / ( goodj[i].L.E() - goodj[i].L.Pz()));
		}
		else if ( goodj[i].L.Pt() > subleadingj_pt)
		{
		    subleadingj_pt  = goodj[i].L.Pt();
		    subleadingj_eta = goodj[i].L.Eta();
		    subleadingj_phi = goodj[i].L.Phi();
		    subleadingj_e   = goodj[i].L.E();
		    jj_index.second = i;
		    subleadingj_rap = 0.5*log( (goodj[i].L.E() + goodj[i].L.Pz()) / ( goodj[i].L.E() - goodj[i].L.Pz()));
		}
	    }
	    if( subleadingj_pt > 0)
	    {
		delta_eta_jj = fabs( leadingj_eta - subleadingj_eta );
		delta_R_jj = sqrt( pow((leadingj_eta -subleadingj_eta),2) + pow((leadingj_phi - subleadingj_phi), 2));
		TLorentzVector jj = goodj[jj_index.first].L + goodj[jj_index.second].L;
		jj_pt  = jj.Pt();
		jj_eta = jj.Eta();
		jj_phi = jj.Phi();
		jj_e   = jj.E();
		jj_m   = jj.M();
	    }

	    if( subleadingj_pt > 0 && ll_good.size() == 2 )
	    {
		delta_R_zz_jj = sqrt( pow((zz_eta - jj_eta), 2) + pow((zz_phi - jj_phi), 2));
		centrality = (zz_rap - (leadingj_rap + subleadingj_rap) / 2) / fabs( leadingj_rap - subleadingj_rap );
	    }

	    TreeStrVar["filename"]["Value"]     = m_FileName;
	    TreeIntVar["signal_bkg"]["Value"]   = SB::_signal;
	    TreeIntVar["run"]["Value"]          = run;
	    TreeUloVar["event"]["Value"]        = event;
	    TreeIntVar["njet"]["Value"]         = njet;
	    TreeIntVar["nmuon"]["Value"]        = nmuon;
	    TreeIntVar["nele"]["Value"]         = nele;
	    TreeFltVar["pileupweight"]["Value"] = pileWeight;
	    TreeFltVar["mcweight"]["Value"]     = mcWeight;
	    TreeFltVar["lep_sf"]["Value"]       = lep_sf;
	    TreeFltVar["trigsf"]["Value"]       = trigger_sf;
	    TreeStrVar["channel"]["Value"]      = channel;

	    TreeFltVar["l1_pt"]["Value"]        = l1_pt;
	    TreeFltVar["l1_eta"]["Value"]       = l1_eta;
	    TreeFltVar["l1_phi"]["Value"]       = l1_phi;
	    TreeFltVar["l1_e"]["Value"]         = l1_e;
	    TreeFltVar["l2_pt"]["Value"]        = l2_pt;
	    TreeFltVar["l2_eta"]["Value"]       = l2_eta;
	    TreeFltVar["l2_phi"]["Value"]       = l2_phi;
	    TreeFltVar["l2_e"]["Value"]         = l2_e;
	    TreeFltVar["l3_pt"]["Value"]        = l3_pt;
	    TreeFltVar["l3_eta"]["Value"]       = l3_eta;
	    TreeFltVar["l3_phi"]["Value"]       = l3_phi;
	    TreeFltVar["l3_e"]["Value"]         = l3_e;
	    TreeFltVar["l4_pt"]["Value"]        = l4_pt;
	    TreeFltVar["l4_eta"]["Value"]       = l4_eta;
	    TreeFltVar["l4_phi"]["Value"]       = l4_phi;
	    TreeFltVar["l4_e"]["Value"]         = l4_e;
	    TreeFltVar["leadingj_pt"]["Value"]  = leadingj_pt;
	    TreeFltVar["leadingj_eta"]["Value"] = leadingj_eta;
	    TreeFltVar["leadingj_phi"]["Value"] = leadingj_phi;
	    TreeFltVar["leadingj_e"]["Value"]   = leadingj_e;
	    TreeFltVar["subleadingj_pt"]["Value"]  = subleadingj_pt;
	    TreeFltVar["subleadingj_eta"]["Value"] = subleadingj_eta;
	    TreeFltVar["subleadingj_phi"]["Value"] = subleadingj_phi;
	    TreeFltVar["subleadingj_e"]["Value"]   = subleadingj_e;
	    TreeFltVar["z1_pt"]["Value"]        = z1_pt;
	    TreeFltVar["z1_eta"]["Value"]       = z1_eta;
	    TreeFltVar["z1_phi"]["Value"]       = z1_phi;
	    TreeFltVar["z1_e"]["Value"]         = z1_e;
	    TreeFltVar["z1_m"]["Value"]         = z1_m;
	    TreeFltVar["z2_pt"]["Value"]        = z2_pt;
	    TreeFltVar["z2_eta"]["Value"]       = z2_eta;
	    TreeFltVar["z2_phi"]["Value"]       = z2_phi;
	    TreeFltVar["z2_e"]["Value"]         = z2_e;
	    TreeFltVar["z2_m"]["Value"]         = z2_m;
	    TreeFltVar["zz_pt"]["Value"]        = zz_pt;
	    TreeFltVar["zz_eta"]["Value"]       = zz_eta;
	    TreeFltVar["zz_phi"]["Value"]       = zz_phi;
	    TreeFltVar["zz_e"]["Value"]         = zz_e;
	    TreeFltVar["zz_m"]["Value"]         = zz_m;
	    TreeFltVar["jj_pt"]["Value"]        = jj_pt;
	    TreeFltVar["jj_eta"]["Value"]       = jj_eta;
	    TreeFltVar["jj_phi"]["Value"]       = jj_phi;
	    TreeFltVar["jj_e"]["Value"]         = jj_e;
	    TreeFltVar["jj_m"]["Value"]         = jj_m;
	    TreeFltVar["delta_eta_jj"]["Value"] = delta_eta_jj;
	    TreeFltVar["delta_R_ll_1"]["Value"] = delta_R_ll_1;
	    TreeFltVar["delta_R_ll_2"]["Value"] = delta_R_ll_2;
	    TreeFltVar["delta_R_zz"]["Value"]   = delta_R_zz;
	    TreeFltVar["delta_R_jj"]["Value"]   = delta_R_jj;
	    TreeFltVar["delta_R_zz_jj"]["Value"]= delta_R_zz_jj;
	    TreeFltVar["centrality"]["Value"]   = centrality;
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
    cout << "fiducial events: " << m_fiducial << endl;
    cout << "selection efficiency: " << (m_fiducial/(m_eventCounter + 0.0) ) << endl;
    return EL::StatusCode::SUCCESS;
}

EL::StatusCode MyxAODAnalysis :: histFinalize ()
{
    return EL::StatusCode::SUCCESS;
}
