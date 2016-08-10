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

#include <xAODMissingET/MissingET.h>
#include <xAODMissingET/MissingETAuxContainer.h>
#include <xAODMissingET/MissingETContainer.h>
#include <xAODMissingET/MissingETComponentMap.h>

#include "xAODEventInfo/EventInfo.h"

#include "xAODTruth/TruthEvent.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthVertex.h"
#include "xAODTruth/TruthVertexContainer.h"

#include <QuickAna/QuickAna.h>

#include "xAODRootAccess/TStore.h"
#include "xAODCore/ShallowCopy.h"

#include "Initialize.h"
#include "Counting.h"


using namespace std;



// this is needed to distribute the algorithm to the workers
ClassImp(MyxAODAnalysis)	



MyxAODAnalysis :: MyxAODAnalysis ()
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
  //m_setSysList=new std::vector<std::string>();
}

MyxAODAnalysis :: MyxAODAnalysis (string treename="physics")
{

  set = treename;
  //m_setSysList=new std::vector<std::string>();
  //

}



EL::StatusCode MyxAODAnalysis :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.
  EL::OutputStream output1("tree_output");
  job.outputAdd (output1);

  EL::OutputStream output2("hist_output");
  job.outputAdd (output2);
  
  EL::OutputStream output3("truth_hist_output");
  job.outputAdd (output3); // add by weibin for truth output

  EL::OutputStream output4("cutflow");
  job.outputAdd (output4);  // add output files  --Zhang

  
  job.useXAOD ();   // filetype is xAOD, tell EventLoop that we actually want to use the xAODRootAccess   --Zhang

  // let's initialize the algorithm to use the xAODRootAccess package
  xAOD::Init( "MyxAODAnalysis" ).ignore(); // call before opening first file   
  // xAOD::Init(const char * appname = "xAOD::Init") return a TReturnCode type var, which call the ignore method to ignore the return code, mark it as checked, Why so??? --Zhang

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: histInitialize ()
{
    // Here you do everything that needs to be done at the very
    // beginning on each worker node, e.g. create histograms and 
    // output trees. This method gets called before any input files
    // are connected.   --Zhang
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: fileExecute ()
{
    // Here you do everything you needs to be done exactly once for 
    // every single file, e.g. collect a list of all lumi-blocks
    // processed     -- Zhnag
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: changeInput (bool firstFile)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}

// < move the "initialize ()" to a separate file "initialize.h"

EL::StatusCode MyxAODAnalysis :: execute ()
{
  m_eventCounter++;
//  if(m_eventCounter == 20)
//  {
//      return EL::StatusCode::FAILURE;
//  }
  if(m_eventCounter==1 || m_eventCounter%1000==0)  // for monitering the dealing precess
    cout << "event counter " << m_eventCounter << endl;

  // get eventinfo  
  const xAOD::EventInfo* eventInfo = 0;
  if( !m_event->retrieve( eventInfo, "EventInfo").isSuccess() )
  {
    Error("execute ()", "Failed to retrieve EventInfo. Exiting." );
    return EL::StatusCode::FAILURE;
  }

//  float aveIntPC = eventInfo->actualInteractionsPerCrossing() + eventInfo->averageInteractionsPerCrossing();
  uint32_t run;
  unsigned long long event;
  float mu = eventInfo->actualInteractionsPerCrossing(); 
  vector<float> vw;
  float pileWeight=1.0, mcWeight = 1.0;
  if(isMC) {
    run = eventInfo->mcChannelNumber(); // Channel Number
    event = eventInfo->mcEventNumber();
    vw = eventInfo->mcEventWeights(); //evnetweights, why is a vector ??? more than one weight ??? --Zhang 
    // SETNAME[0] = "physics"
    if(SETTING[SETNAME[0]]["doweight"]==1) mcWeight = vw.size()>0?vw[0]:1.0;
    if(mcWeight>5) mcWeight = 1;
  }else {   
    if(!m_grl->passRunLB(*eventInfo)){
       return EL::StatusCode::SUCCESS; // go to next event
    }
    run = eventInfo->runNumber();
    event = eventInfo->eventNumber();  // for data, weight is always 1
  }
  long double Weight = mcWeight;
  //cout << "weight " <<  mcweight << endl;
  //
  // END OF EVENTINFO
  // 

  //
  // muon, electron, jet, met
  //
  VOmuon goodm; 
  VOelectron goode;
  VOmet goodmet;
  VOjet goodj;
//  VOmet goodmet;  

  // get Vertex
  const xAOD::Vertex *pv(0);
  const xAOD::VertexContainer* vertices = 0;
  if( !m_event->retrieve( vertices, "PrimaryVertices" ).isSuccess() )
  {
    Error("execute ()", "Failed to retrieve verteices. Exiting." );
    return EL::StatusCode::FAILURE;
  } //  primaryVertices
  for ( const auto* const vtx_itr : *vertices )  // const auto* const ??? --Zhang
  {
    if (vtx_itr->vertexType() != xAOD::VxType::VertexType::PriVtx) continue;
    else { pv = vtx_itr; break;}
  }   // until we find the first priamryVertex, why only first is enough??? --Zhang
  if(!pv) return EL::StatusCode::SUCCESS;
  double pz0 = pv->z();         // assignment of pv->z();


  // get TruthEvent
  const xAOD::TruthEventContainer* TruthEvtContainer = 0;
  if( !m_event->retrieve( TruthEvtContainer, "TruthEvents").isSuccess() ){
    Error("excute()", "Failed to retrieve Truth info. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  TLorentzVector truth_Z;
  TLorentzVector truth_Higgs;
  int Z_num = 0;
  int truth_type = -1;
  float Truth_Higgs_Mass = 0;
  xAOD::TruthEventContainer::const_iterator trEvt_itr = TruthEvtContainer->begin();
  xAOD::TruthEventContainer::const_iterator trEvt_end = TruthEvtContainer->end();
  // what's the meaning of the following code looping through TruthEvtContainer;
  for( ; trEvt_itr != trEvt_end; ++trEvt_itr ) { 
    int nTruPart = (*trEvt_itr)->nTruthParticles();
    for(int j=0; j<nTruPart; j++) { 
      const xAOD::TruthParticle* trPart = (*trEvt_itr)->truthParticle(j);
      if(!trPart) continue;
      int status = trPart->status();
      int PDG = trPart->pdgId();  
      bool ProdVx = trPart->hasProdVtx();  
      bool DecayVx = trPart->hasDecayVtx(); 

      if (PDG==23){  // Z0 particle
        if(!DecayVx) continue;
        int nChi = trPart->nChildren();
        for (int i = 0; i < nChi; i++){
          const xAOD::TruthParticle* chiPart = trPart->child(i);
          if(!chiPart) continue;
          int PDGchi = chiPart->pdgId();
          int statuschi = chiPart->status(); 
          if(fabs(PDGchi)==11) // e-
            truth_type = 0;
          if(fabs(PDGchi)==13)    // u-
            truth_type = 1; // truth_Type ???
          if(fabs(PDGchi)==11 || fabs(PDGchi)==13){
            truth_Z = trPart->p4(); 
          }
        } 
      } // end of Z0 decay mode
      if (PDG==25){  // H0(Higgs)
        truth_Higgs = trPart->p4();
        Truth_Higgs_Mass = truth_Higgs.M()*0.001;
        //cout<<"Truth Higgs Mass: "<<truth_Higgs.M()*0.001<<endl;
      } 
    } 
  } // END OF EVENTCONTAINER loop


  //
  // Truth Particles Dealing 
  //

  // retrieve truth leptons
  const xAOD::TruthParticleContainer* mc_particle = 0;
  if (isMC) {
    if ( ! m_event->retrieve(mc_particle,"TruthParticles").isSuccess() ) {
      Error("execute()", "Failed to retrieve TruthParticleContainer container. Exiting.");
      return EL::StatusCode::FAILURE;
    }
  }

  ClearFlags(FLAG_Truth_cut_temp);  
  ClearFlags(FLAG_Truth_cut);       
  ClearVariables(TruthHistVar);    

  VOtruth_muon truth_muons;
  truth_muons.clear();
//  double index = 0;
  if(isMC) {
    for (auto truth = mc_particle->begin(); truth != mc_particle->end(); ++truth ) {

      bool findZ = false;
      int pdg = (*truth)->pdgId();
      int status = (*truth)->status();
      int barcode = (*truth)->barcode();

      // choose prompt(barcode) stable (status == 1) muon (pdg == 13)
      if(abs(pdg)!=13 || status!=1 || barcode > 1.e5) continue;

      /*
      // must have a parent and eventually pointing back to the Z boson
      int nparent = (*truth)->nParents(); // could any particle has more than one parent ???
      if(nparent!=1) continue;

      const xAOD::TruthParticle* theparent = (*truth)->parent(0);

      while(nparent==1) {
	int pdgparent = theparent->pdgId();
	int status = theparent->status();
	cout << "Parent pid: " << pdgparent << "\t" << "status: " << status << endl;
	nparent = theparent->nParents();
	theparent = theparent->parent(0);
      } 

       */
      // for Sherpa Generator, there may be no Z particles in the process
      int nparent = (*truth)->nParents();
      const xAOD::TruthParticle* theparent;

      if(nparent > 0)
      {
        for(int i=0; i<nparent; i++)
        {
          theparent = (*truth)->parent(i);
          int parent_pdg = theparent->pdgId();
          int parent_status = theparent->status();
          int parent_barcode = theparent->barcode();

          if(parent_pdg == pdg)
          {
            if(parent_status == 3 && barcode < 1.e5)
            {
              findZ = true;
              break;
            }
            else 
            {
              nparent = theparent->nParents();
	      break;
            }
          }
	  if(i == (nparent - 1))
	  {
	    nparent = 0;
	  }
        }
      }

      if(!findZ)
      {
        while(nparent > 0)
        {
          for(int i=0; i<nparent; i++)
          {
            const xAOD::TruthParticle* temp = theparent->parent(i);
            int parent_pdg = temp->pdgId();
            int parent_status = temp->status();
            int parent_barcode = temp->barcode();

            if(parent_pdg == pdg)
            {
              if(parent_status == 3 && barcode < 1.e5)
              {
                findZ = true;
                break;
              }
              else 
              {
                nparent = temp->nParents();
                break;
              }
            }
	    if( i == (nparent - 1))
	    {
	      nparent = 0;
	    }
          }
	  if(findZ)
	  {
	    break;
	  }
        }
      }

      if(!findZ) continue;

      // find Z
      int charge = (*truth)->charge();
      double pt = (*truth)->pt();
      double eta = (*truth)->eta();
      double phi = (*truth)->phi();
//      if( pt < 5.e3 || abs(eta) > 2.5) continue;

      TLorentzVector temp_muon;
      temp_muon.SetPtEtaPhiM(pt, eta, phi, 105.65837);

      // dressing
      for(auto truth2 = mc_particle->begin(); truth2 != mc_particle->end(); ++truth2){
	int pdg2 = (*truth2)->pdgId();
	int status2 = (*truth2)->status();
	int barcode2 = (*truth2)->barcode();
	if(pdg2!=22 || status2!=1 || barcode2>1.e5) continue;

	TLorentzVector temp_photon;
	double pt = (*truth2)->pt() ;
	double eta = (*truth2)->eta();
	double phi = (*truth2)->phi();
	temp_photon.SetPtEtaPhiM(pt, eta, phi, 0);
	double dR = temp_photon.DeltaR(temp_muon);
	if(dR<0.1) temp_muon += temp_photon;
      }

      OBJ_TRUTH_MUON muonInfo;
      muonInfo.charge = charge;
      muonInfo.pt = temp_muon.Pt();
      muonInfo.eta = temp_muon.Eta();
      muonInfo.phi = temp_muon.Phi();
      muonInfo.m = temp_muon.M();
      truth_muons.push_back(muonInfo);

      /*
      cout << "Index: " << index << "\t" << pdg << "\t" << status << "\t" << barcode << "\t" << pt << "\t" << eta << "\t" << phi <<  endl;
      nparent = (*truth)->nParents();
      for (int i=0; i<nparent; i++)
      {
	const xAOD::TruthParticle* theparent = (*truth)->parent(i);
	int parent_pdg = theparent->pdgId();
	int parent_status = theparent->status();
	int parent_barcode = theparent->barcode();
//	int parent_index = theparent - mc_particle->begin();
        double parent_pt = theparent->pt() ;
        double parent_eta = theparent->eta();
        double parent_phi = theparent->phi();
	cout << "\t\t" << i << " : " << parent_pdg << "\t" << parent_status << "\t" << parent_barcode << "\t" << parent_pt << "\t" << parent_eta << "\t" << parent_phi << endl;
	cout << endl;
      }
       */
    }
  }

//  cout << "muons: " << truth_muons.size() << endl;
  VOtruth_electron truth_electrons;
  truth_electrons.clear();
//  index = 0;
//  cout << "Event: " << m_eventCounter << endl;
//  const xAOD::TruthParticle* begin = mc_particle->begin();
  if ( isMC ) {
    for( auto truth = mc_particle->begin(); truth != mc_particle->end(); ++truth) {
      
      int pdg = (*truth)->pdgId();
      int status = (*truth)->status();
      int barcode = (*truth)->barcode();
      
      int charge = (*truth)->charge();
      double pt = (*truth)->pt() ;
      double eta = (*truth)->eta();
      double phi = (*truth)->phi();

      int nparent = (*truth)->nParents();
      /*
      index++;
      cout << "Index: " << index << "\t" << pdg << "\t" << status << "\t" << barcode << "\t" << pt << "\t" << eta << "\t" << phi <<  endl;
      for (int i=0; i<nparent; i++)
      {
	const xAOD::TruthParticle* theparent = (*truth)->parent(i);
	int parent_pdg = theparent->pdgId();
	int parent_status = theparent->status();
	int parent_barcode = theparent->barcode();
//	int parent_index = theparent - mc_particle->begin();
        double parent_pt = theparent->pt() ;
        double parent_eta = theparent->eta();
        double parent_phi = theparent->phi();
	cout << "\t\t" << i << " : " << parent_pdg << "\t" << parent_status << "\t" << parent_barcode << "\t" << parent_pt << "\t" << parent_eta << "\t" << parent_phi << endl;
      }
       */
    
      // choose stable prompt electron( pdg == 11)
//      if(abs(pdg)!=11 || status!=1 || barcode > 1.e5) continue;
      bool findZ = false;
      if(abs(pdg)!=11 || status != 1 || barcode > 1.e5) continue;

//      int nparent = (*truth)->nParents();
      const xAOD::TruthParticle* theparent;

      if(nparent > 0)
      {
        for(int i=0; i<nparent; i++)
        {
          theparent = (*truth)->parent(i);
          int parent_pdg = theparent->pdgId();
          int parent_status = theparent->status();
          int parent_barcode = theparent->barcode();

          if(parent_pdg == pdg)
          {
            if(parent_status == 3 && barcode < 1.e5)
            {
              findZ = true;
              break;
            }
            else 
            {
              nparent = theparent->nParents();
	      break;
            }
          }
	  if(i == (nparent - 1))
	  {
	    nparent = 0;
	  }
        }
      }

      if(!findZ)
      {
        while(nparent > 0)
        {
          for(int i=0; i<nparent; i++)
          {
            const xAOD::TruthParticle* temp = theparent->parent(i);
            int parent_status = temp->status();
            int parent_pdg = temp->pdgId();
            int parent_barcode = temp->barcode();

            if(parent_pdg == pdg)
            {
              if(parent_status == 3 && barcode < 1.e5)
              {
                findZ = true;
                break;
              }
              else 
              {
                nparent = temp->nParents();
                break;
              }
            }
	  if(i == (nparent - 1) )
	  {
	    nparent = 0;
	  }
          }
	  if(findZ)
	  {
	    break;
	  }
        }
      }

      if(!findZ) continue;

      /*
      // Z parent
      int nparent = (*truth)->nParents();
      if(nparent!=1) continue;

      bool findZ = false;
      const xAOD::TruthParticle* theparent = (*truth)->parent(0);

      while(nparent == 1) {
        int pdgparent = theparent->pdgId();
	if(pdgparent == 23) findZ = true;
	else {
	  nparent = theparent->nParents();
	  theparent = theparent->parent(0);
	}
	if(findZ) break;
      }

      if (!findZ) continue;
       */

//      if(pt<5.e3 || abs(eta)>2.47) continue;

      OBJ_TRUTH_ELECTRON elecInfo;

      TLorentzVector temp_electron;
      temp_electron.SetPtEtaPhiM(pt, eta, phi, 0.510999);

      // dressing 
      for( auto truth2 = mc_particle->begin(); truth2 != mc_particle->end(); ++truth2 ) {
        int pdg2 = (*truth2)->pdgId();
	int status2 = (*truth2)->status();
	int barcode2 = (*truth2)->barcode();
	if(pdg2!=22 || status2!=1 || barcode2>1.e5) continue;

	TLorentzVector temp_photon;
	double pt = (*truth2)->pt();
	double eta = (*truth2)->eta();
	double phi = (*truth2)->phi();
	temp_photon.SetPtEtaPhiM(pt, eta, phi, 0);
	double dR = temp_photon.DeltaR(temp_electron);
	if(dR<0.1) temp_electron+=temp_photon;
      }

      elecInfo.charge = charge;
      elecInfo.pt = temp_electron.Pt();
      elecInfo.eta = temp_electron.Eta();
      elecInfo.phi = temp_electron.Phi();
      elecInfo.m = temp_electron.M();
      truth_electrons.push_back(elecInfo);

      /*
      cout << "Index: " << index << "\t" << pdg << "\t" << status << "\t" << barcode << "\t" << pt << "\t" << eta << "\t" << phi <<  endl;
      nparent = (*truth)->nParents();
      for (int i=0; i<nparent; i++)
      {
	const xAOD::TruthParticle* theparent = (*truth)->parent(i);
	int parent_pdg = theparent->pdgId();
	int parent_status = theparent->status();
	int parent_barcode = theparent->barcode();
//	int parent_index = theparent - mc_particle->begin();
        double parent_pt = theparent->pt() ;
        double parent_eta = theparent->eta();
        double parent_phi = theparent->phi();
	cout << "\t\t" << i << " : " << parent_pdg << "\t" << parent_status << "\t" << parent_barcode << "\t" << parent_pt << "\t" << parent_eta << "\t" << parent_phi << endl;
	cout << endl;
      }
       */
    }
  }

//  int lep_num = truth_muons.size() + truth_electrons.size();
//  cout << "lep_num: " << lep_num << endl;
//  cout << endl;

  // retrieve truth jets
  const xAOD::JetContainer* TruthJets = 0;
  if(! m_event->retrieve(TruthJets, "AntiKt4TruthJets").isSuccess())
  {
    Error("excute()", "Failed to retrieve Truth Jets info. Exiting.");
    return EL::StatusCode::FAILURE;
  }

  vector<OBJ_TRUTH_JET> good_truth_jets;
  good_truth_jets.clear();
  xAOD::JetContainer::const_iterator trJets_itr = TruthJets->begin();
  xAOD::JetContainer::const_iterator trJets_end = TruthJets->end();
  for(; trJets_itr != trJets_end; ++trJets_itr)
  {
//    if(((*trJets_itr)->pt() > 25.e3 && fabs((*trJets_itr)->eta()) < 2.4) || ((*trJets_itr)->pt() > 30.e3 && fabs((*trJets_itr)->eta()) < 4.5))
    if((*trJets_itr)->pt() > 25.e3)
    {
      OBJ_TRUTH_JET jetInfo;

      jetInfo.px = (*trJets_itr)->px();
      jetInfo.py = (*trJets_itr)->py();
      jetInfo.pz = (*trJets_itr)->pz();
      jetInfo.E =  (*trJets_itr)->e();

      good_truth_jets.push_back(jetInfo);
    }
  }


  // a simple way dealing with truth particle
    string TRUTH_CHANNEL = "All";
    int truth_ele_num = -1, truth_muon_num = -1, truth_lep_num = -1;
    double truth_Pt1 = -9999.0, truth_Pt2 = -9999.0, truth_Pt3 = -9999.0, truth_Pt4 = -9999.0;
    double truth_Eta1 = -9999.0, truth_Eta2 = -9999.0, truth_Eta3 = -9999.0, truth_Eta4 = -9999.0;
    double truth_Phi1 = -9999.0, truth_Phi2 = -9999.0, truth_Phi3 = -9999.0, truth_Phi4 = -9999.0;
    double truth_Z1 = -9999.0, truth_Z2 = -9999.0;
    double truth_ZZ = -9999.0;
    int truth_jet_num = -1;
    double truth_leadingJet_Pt = -9999.0, truth_subleadingJet_Pt = -9999.0;
    double truth_leadingJet_Eta = -9999.0, truth_subleadingJet_Eta = -9999.0;
    double truth_leadingJet_Phi = -9999.0, truth_subleadingJet_Phi = -9999.0;
    double truth_Mjj = -9999.0;
    double truth_Delta_Jet_Eta = -9999.0;
    double truth_leadingJet_Rap = -9999.0, truth_subleadingJet_Rap = -9999.0;
    double truth_4leptons_Rap = -9999.0;
    double truth_Centrality = -9999.0;
    bool Pass_Truth_4Leptons_mmmm = false;
    bool Pass_Truth_4Leptons_eeee = false;
    bool Pass_Truth_4Leptons_eemm = false;
    bool Pass_Truth_4Leptons = false;
    bool Pass_Truth_LepPt = false;
    bool Pass_Truth_LepEta = false;
    bool Pass_Truth_LeadingLepPt = false;
    bool Pass_Truth_2Pairs = false;
    bool Pass_Truth_Mll1 = false;
    bool Pass_Truth_Mll2 = false;
    bool Pass_Truth_NumJet = false;
    bool Pass_Truth_Mjj = false;
    bool Pass_Truth_Delta_Jet_Eta = false;
    bool Pass_Truth_Centrality = false;

    truth_ele_num = truth_electrons.size();
    truth_muon_num = truth_muons.size();
    truth_lep_num = truth_ele_num + truth_muon_num;

    if(truth_muons.size() == 4 && truth_electrons.size() == 0) {
      TRUTH_CHANNEL = "mmmm";
      Pass_Truth_4Leptons_mmmm = true;

      truth_Pt1 = truth_muons[0].pt;
      truth_Pt2 = truth_muons[1].pt;
      truth_Pt3 = truth_muons[2].pt;
      truth_Pt4 = truth_muons[3].pt;

      truth_Eta1 = truth_muons[0].eta;
      truth_Eta2 = truth_muons[1].eta;
      truth_Eta3 = truth_muons[2].eta;
      truth_Eta4 = truth_muons[3].eta;

      truth_Phi1 = truth_muons[0].phi;
      truth_Phi2 = truth_muons[1].phi;
      truth_Phi3 = truth_muons[2].phi;
      truth_Phi4 = truth_muons[3].phi;
    }
    else if (truth_electrons.size() == 4 && truth_muons.size() == 0) {
      TRUTH_CHANNEL = "eeee";
      Pass_Truth_4Leptons_eeee = true;

      truth_Pt1 = truth_electrons[0].pt;
      truth_Pt2 = truth_electrons[1].pt;
      truth_Pt3 = truth_electrons[2].pt;
      truth_Pt4 = truth_electrons[3].pt;

      truth_Eta1 = truth_electrons[0].eta;
      truth_Eta2 = truth_electrons[1].eta;
      truth_Eta3 = truth_electrons[2].eta;
      truth_Eta4 = truth_electrons[3].eta;

      truth_Phi1 = truth_electrons[0].phi;
      truth_Phi2 = truth_electrons[1].phi;
      truth_Phi3 = truth_electrons[2].phi;
      truth_Phi4 = truth_electrons[3].phi;
    }
    else if (truth_electrons.size() == 2 && truth_muons.size() == 2) {
      TRUTH_CHANNEL = "eemm";
      Pass_Truth_4Leptons_eemm = true;

      truth_Pt1 = truth_electrons[0].pt;
      truth_Pt2 = truth_electrons[1].pt;
      truth_Pt3 = truth_muons[0].pt;
      truth_Pt4 = truth_muons[1].pt;
    
      truth_Eta1 = truth_electrons[0].eta;
      truth_Eta2 = truth_electrons[1].eta;
      truth_Eta3 = truth_muons[0].eta;
      truth_Eta4 = truth_muons[1].eta;
    
      truth_Phi1 = truth_electrons[0].phi;
      truth_Phi2 = truth_electrons[1].phi;
      truth_Phi3 = truth_muons[0].phi;
      truth_Phi4 = truth_muons[1].phi;
    }

    {
      int index[4] = {0,1,2,3};
      double Pt[4] = {truth_Pt1, truth_Pt2, truth_Pt3, truth_Pt4};
      double Eta[4] = {truth_Eta1, truth_Eta2, truth_Eta3, truth_Eta4};
      double Phi[4] = {truth_Phi1, truth_Phi2, truth_Phi3, truth_Phi4};
      for ( int i=0; i<4; i++) {
        for ( int j=i+1; j<4; j++) {
	  if(Pt[index[j]] > Pt[index[i]]) {
	    int temp =  index[i];
	    index[i] = index[j];
	    index[j] = temp;
	  }
	}
      }
      truth_Pt1 = Pt[index[0]];
      truth_Pt2 = Pt[index[1]];
      truth_Pt3 = Pt[index[2]];
      truth_Pt4 = Pt[index[3]];
      truth_Eta1 = Eta[index[0]];
      truth_Eta2 = Eta[index[1]];
      truth_Eta3 = Eta[index[2]];
      truth_Eta4 = Eta[index[3]];
      truth_Phi1 = Phi[index[0]];
      truth_Phi2 = Phi[index[1]];
      truth_Phi3 = Phi[index[2]];
      truth_Phi4 = Phi[index[3]];
    }

    if(Pass_Truth_4Leptons_mmmm || Pass_Truth_4Leptons_eeee || Pass_Truth_4Leptons_eemm) 
    {
      Pass_Truth_4Leptons = true;
      m_truth_4leptons++;
      w_truth_4leptons += Weight;
    }

    if(Pass_Truth_4Leptons && truth_Eta1 < 2.5 && truth_Eta2 < 2.5 && truth_Eta3 < 2.5 && truth_Eta4 < 2.5)
    {
      Pass_Truth_LepEta = true;
      m_truth_lepeta++;
      w_truth_lepeta += Weight;
    }

    if(Pass_Truth_LepEta && truth_Pt4 > 7.e3)
    {
      Pass_Truth_LepPt = true;
      m_truth_leppt++;
      w_truth_leppt += Weight;
    }


    if(Pass_Truth_LepPt && truth_Pt1 > 25.e3)
    {
      Pass_Truth_LeadingLepPt = true;
      m_truth_leadingleppt++;
      w_truth_leadingleppt += Weight;
    }

    if(Pass_Truth_LeadingLepPt)
    {
      TLorentzVector total,lep;
      for(int i=0; i<truth_muons.size(); i++)
      {
	lep.SetPtEtaPhiM(truth_muons[i].pt, truth_muons[i].eta, truth_muons[i].phi, truth_muons[i].m);
	total += lep;
      }
      for(int i=0; i<truth_electrons.size(); i++)
      {
	lep.SetPtEtaPhiM(truth_electrons[i].pt, truth_electrons[i].eta, truth_electrons[i].phi, truth_electrons[i].m);
	total += lep;
      }
      truth_ZZ = total.M();
      truth_4leptons_Rap = total.Rapidity();
    }

    vector<Truth_Pair> truth_pairs;
    truth_pairs.clear();
    for (int i=0; i<truth_muons.size(); i++)
    {
      for(int j=i+1; j<truth_muons.size(); j++)
      {
	if(truth_muons[i].charge * truth_muons[j].charge == 1) continue;
	Truth_Pair temp_pair;
	temp_pair.flavour = 1; // flavour: muon
	temp_pair.opposign = true;
	temp_pair.index.push_back(i);
	temp_pair.index.push_back(j);
	TLorentzVector muon1, muon2, muons;
	muon1.SetPtEtaPhiM(truth_muons[i].pt, truth_muons[i].eta, truth_muons[i].phi, truth_muons[i].m);
	muon2.SetPtEtaPhiM(truth_muons[j].pt, truth_muons[j].eta, truth_muons[j].phi, truth_muons[j].m);
	muons = muon1 + muon2;
	temp_pair.mass = muons.M();
	truth_pairs.push_back(temp_pair);
      }
    }

    for (int i=0; i<truth_electrons.size(); i++)
    {
      for(int j=i+1; j<truth_electrons.size(); j++)
      {
	if(truth_electrons[i].charge*truth_electrons[j].charge == 1) continue;
	Truth_Pair temp_pair;
	temp_pair.flavour = 2; // flavour: electron
	temp_pair.opposign = true;
	temp_pair.index.push_back(i);
	temp_pair.index.push_back(j);
	TLorentzVector electron1, electron2, electrons;
	electron1.SetPtEtaPhiM(truth_electrons[i].pt, truth_electrons[i].eta, truth_electrons[i].phi, truth_electrons[i].m);
	electron2.SetPtEtaPhiM(truth_electrons[j].pt, truth_electrons[j].eta, truth_electrons[j].phi, truth_electrons[j].m);
	electrons = electron1 + electron2;
	temp_pair.mass = electrons.M();
	truth_pairs.push_back(temp_pair);
      }
    }

    /*
    for(int i=0; i<truth_pairs.size(); i++)
    {
      for(int j=0; j>truth_pairs.size(); j++)
      {
	if(truth_pairs[i].flavour != truth_pairs[j].flavour) 
	{
	  Pass2Pairs = true;
	  break;
	}
	else
	{
	  int indexi_0 = -1, indexi_1 = -1;
	  int indexj_0 = -1, indexj_1 = -1;
	  indexi_0 = truth_pairs[i].index[0];
	  indexi_1 = truth_pairs[i].index[1];
	  indexj_0 = truth_pairs[j].index[0];
	  indexj_1 = truth_pairs[j].index[1];
	  if(indexi_0 != indexj_0 && indexi_0 != indexj_1 && indexi_1 != indexj_0 && indexi_0 != indexj_1)
	  {
	    Pass2Pairs = true;
	    break;
	  }
	}
      }
      if(Pass2Pairs) break;
    } // 2 same-flavour opposign charge pairs
    */

    vector<Truth_Pair> good_truth_pairs;
    good_truth_pairs.clear();
    int tag = -1, sectag = -1;
    double tag_dm = 99999999999.e3, sectag_dm = 999999999999.e3;
    for(int i=0; i<truth_pairs.size(); i++)
    {
      int indexi_0 = truth_pairs[i].index[0];
      int indexi_1 = truth_pairs[i].index[1];
      double temp_tag_dm = abs(truth_pairs[i].mass - ZMass);
      double temp_sectag_dm = 999999999999.e3;
      if ((temp_sectag_dm + temp_tag_dm) < (tag_dm + sectag_dm) )
      {
        tag_dm = temp_tag_dm;
        tag = i;
      }
      for (int j=i+1; j<truth_pairs.size(); j++)
      {
	int indexj_0 = truth_pairs[j].index[0];
	int indexj_1 = truth_pairs[j].index[1];
	if(truth_pairs[j].flavour == truth_pairs[i].flavour)
	{
	  if(indexj_0 == indexi_0 || indexj_0 == indexi_1 || indexj_1 == indexi_0 || indexj_1 == indexi_0)
	    continue;
	}
	temp_sectag_dm = abs(truth_pairs[j].mass - ZMass);
	if ((temp_sectag_dm + temp_tag_dm) < (tag_dm + sectag_dm) )
	{
          tag_dm = temp_tag_dm;
	  sectag_dm = temp_sectag_dm;
          tag = i;
	  sectag = j;
	}
      }
    }
    
    if(tag != -1)
    {
      truth_Z1 = truth_pairs[tag].mass;
      good_truth_pairs.push_back(truth_pairs[tag]);
      if(sectag != -1)
      {
        if(sectag_dm < tag_dm )
        {
          double temp_dm = sectag_dm;
          tag_dm = sectag_dm;
          sectag_dm = temp_dm;
          tag += sectag;
          sectag = tag - sectag;
          tag = tag - sectag;
        } // assign the index with small deltaM of pairs to tag
        truth_Z1 = truth_pairs[tag].mass;
	truth_Z2 = truth_pairs[sectag].mass;
        good_truth_pairs.push_back(truth_pairs[tag]);
        good_truth_pairs.push_back(truth_pairs[sectag]);
	if(Pass_Truth_LeadingLepPt)
	{
          Pass_Truth_2Pairs = true;
	  m_truth_2pairs++;
	  w_truth_2pairs += Weight;
	}
      }
    }

    // truth jets
    if(good_truth_jets.size() > 1) // have at least 2 good jets
    {
      truth_jet_num = good_truth_jets.size();
      Pass_Truth_NumJet = true;
      TLorentzVector leadingJet, subleadingJet;
      int leading_tag=0, subleading_tag=0;
      double leadingPt=-9999.0, subleadingPt=-9999.0;
      for(int i=0; i<good_truth_jets.size(); i++)
      {
	double pt = sqrt(pow(good_truth_jets[i].px,2) + pow(good_truth_jets[i].py, 2));
	if(pt>leadingPt)
	{
	  subleadingPt = leadingPt;
	  subleading_tag = leading_tag;
	  leadingPt = pt;
	  leading_tag = i;
	  continue;
	}
	if(pt>subleadingPt)
	{
	  subleadingPt = pt;
	  subleading_tag = i;
	  continue;
	}
      }
      leadingJet.SetPxPyPzE(good_truth_jets[leading_tag].px, good_truth_jets[leading_tag].py, good_truth_jets[leading_tag].pz, good_truth_jets[leading_tag].E);
      subleadingJet.SetPxPyPzE(good_truth_jets[subleading_tag].px, good_truth_jets[subleading_tag].py, good_truth_jets[subleading_tag].pz, good_truth_jets[subleading_tag].E);

      truth_leadingJet_Pt = leadingJet.Pt();
      truth_leadingJet_Eta = leadingJet.Eta();
      truth_leadingJet_Phi = leadingJet.Phi();
      truth_subleadingJet_Pt = subleadingJet.Pt();
      truth_subleadingJet_Eta = subleadingJet.Eta();
      truth_subleadingJet_Phi = subleadingJet.Phi();
      truth_Delta_Jet_Eta = fabs(leadingJet.Eta() - subleadingJet.Eta());
      TLorentzVector total = leadingJet + subleadingJet;
      truth_Mjj = total.M();
      truth_leadingJet_Rap = leadingJet.Rapidity();
      truth_subleadingJet_Rap = subleadingJet.Rapidity();

      truth_Centrality = ( truth_4leptons_Rap - (truth_leadingJet_Rap + truth_subleadingJet_Rap)/2)/fabs(truth_leadingJet_Rap - truth_subleadingJet_Rap);
    }

    if(Pass_Truth_2Pairs && tag_dm < 25.e3 )
    {
      Pass_Truth_Mll1 = true;
      m_truth_mll1++;
      w_truth_mll1 += Weight;
      if(sectag_dm < 25.e3)
	Pass_Truth_Mll2 = true;
        m_truth_mll2++;
        w_truth_mll2 += Weight;
    }

    if(Pass_Truth_Mll2 && good_truth_jets.size() > 1)
    {
      Pass_Truth_NumJet = true;
      m_truth_numjets++;
      w_truth_numjets += Weight;
    }

//    if(Pass_Truth_NumJet && truth_Mjj > 500.e3)
//      Pass_Truth_Mjj = true;

    SetFlag(FLAG_Truth_cut,"xAOD","All", 1);
    m_truth_xAOD++;
    w_truth_xAOD += Weight;
    SetFlag(FLAG_Truth_cut,"4Leptons", TRUTH_CHANNEL, Pass_Truth_4Leptons);
    SetFlag(FLAG_Truth_cut,"LepPt", TRUTH_CHANNEL, Pass_Truth_LepPt);
    SetFlag(FLAG_Truth_cut,"LepEta", TRUTH_CHANNEL, Pass_Truth_LepEta);
    SetFlag(FLAG_Truth_cut,"LeadingLepPt", TRUTH_CHANNEL, Pass_Truth_LeadingLepPt);
    SetFlag(FLAG_Truth_cut,"2Pairs", TRUTH_CHANNEL, Pass_Truth_2Pairs);
    SetFlag(FLAG_Truth_cut,"Mll1", TRUTH_CHANNEL, Pass_Truth_Mll1);
    SetFlag(FLAG_Truth_cut,"Mll2", TRUTH_CHANNEL, Pass_Truth_Mll2);
    SetFlag(FLAG_Truth_cut, "NumJets", TRUTH_CHANNEL, Pass_Truth_NumJet);

    TruthHistVar["truth_Pt1"]["Value"] = truth_Pt1;
    TruthHistVar["truth_Pt2"]["Value"] = truth_Pt2;
    TruthHistVar["truth_Pt3"]["Value"] = truth_Pt3;
    TruthHistVar["truth_Pt4"]["Value"] = truth_Pt4;
    TruthHistVar["truth_Eta1"]["Value"] = truth_Eta1;
    TruthHistVar["truth_Eta2"]["Value"] = truth_Eta2;
    TruthHistVar["truth_Eta3"]["Value"] = truth_Eta3;
    TruthHistVar["truth_Eta4"]["Value"] = truth_Eta4;
    TruthHistVar["truth_Phi1"]["Value"] = truth_Phi1;
    TruthHistVar["truth_Phi2"]["Value"] = truth_Phi2;
    TruthHistVar["truth_Phi3"]["Value"] = truth_Phi3;
    TruthHistVar["truth_Phi4"]["Value"] = truth_Phi4;
    TruthHistVar["truth_lep_num"]["Value"] = truth_lep_num;
    TruthHistVar["truth_ele_num"]["Value"] = truth_ele_num;
    TruthHistVar["truth_muon_num"]["Value"] = truth_muon_num;
    TruthHistVar["truth_Z1"]["Value"] = truth_Z1;
    TruthHistVar["truth_Z2"]["Value"] = truth_Z2;
    TruthHistVar["truth_ZZ"]["Value"] = truth_ZZ;
    TruthHistVar["truth_jet_num"]["Value"] = truth_jet_num;
    TruthHistVar["truth_leadingJet_Pt"]["Value"] = truth_leadingJet_Pt;
    TruthHistVar["truth_leadingJet_Eta"]["Value"] = truth_leadingJet_Eta;
    TruthHistVar["truth_leadingJet_Phi"]["Value"] = truth_leadingJet_Phi;
    TruthHistVar["truth_subleadingJet_Pt"]["Value"] = truth_subleadingJet_Pt;
    TruthHistVar["truth_subleadingJet_Eta"]["Value"] = truth_subleadingJet_Eta;
    TruthHistVar["truth_subleadingJet_Phi"]["Value"] = truth_subleadingJet_Phi;
    TruthHistVar["truth_Mjj"]["Value"] = truth_Mjj;
    TruthHistVar["truth_Delta_Jet_Eta"]["Value"] = truth_Delta_Jet_Eta;
    TruthHistVar["truth_Centrality"]["Value"] = truth_Centrality;


    string sysname = "NOMINAL";
    FillTruthHistograms(sysname);

//    cout << "Truth Pt: " << truth_Pt1 << "\t" << truth_Pt2 << "\t"<< truth_Pt3 << "\t"<< truth_Pt4 << endl;
//    cout << "Truth Eta: " << truth_Eta1 << "\t"<< truth_Eta2 << "\t"<< truth_Eta3 << "\t"<< truth_Eta4 << endl;
//    cout << "Truth Phi: " << truth_Phi1 << "\t"<< truth_Phi2 << "\t"<< truth_Phi3 << "\t"<< truth_Phi4 << endl;
  // END

  // systematics 
  for (auto sysListItr : m_sysList){
    string sysname = sysListItr.name();

    if(!isMC && SETTING[SETNAME[0]]["dosys"]==1 && sysname != "") continue;
    if(SETTING[SETNAME[0]]["dosys"]==0 && sysname != "") continue;

    if(sysname == "") sysname = "NOMINAL";

    if(SETTING[SETNAME[0]]["docorr"]==0) sysname="NOCORR";
    // so, actually, only "NOMINAL" systematic are processed.

    // quickAna process
    if( quickAna->applySystematicVariation (sysListItr) == CP::SystematicCode::Ok) {

      ClearFlags(FLAG_cut_temp);  // what's the difference between
      ClearFlags(FLAG_cut);       // these ClearFlags();
      ClearWeight(Evt_Weight);    // ClearWeight() and ClearVariables
      ClearVariables(HistVar);    // functions? --Zhang
      ClearVariables(VVar);
      ClearVariables(V2DVar);
      ClearVariables(TreeIntVar);
      ClearVariables(TreeFltVar);
      ClearVariables(TreeFltVVar);
      goodm.clear();
      goode.clear();
      goodj.clear();
      goodmet.clear();

      VOmuon temp_muon; temp_muon.clear(); 
      VOelectron temp_electron; temp_electron.clear(); 
      VOjet temp_jet; temp_jet.clear(); 

      quickAna->process(*m_event).ignore();   

      auto evtInfo = quickAna->eventinfo();  
      bool passTrig = evtInfo->auxdata<bool>("passAllTrig"); // auxdata ??? --Zhang
      if(SETTING[SETNAME[0]]["doweight"]==1) pileWeight = evtInfo->auxdata<float>("PileupWeight"); 
      //cout << "pile weight is " << pileWeight << endl;
      Weight = Weight*pileWeight; 

      SetFlag(FLAG_cut_temp,"xAOD","All", 1);
      if(SETTING[SETNAME[0]]["doweight"]==1)
        SetWeight(Evt_Weight, "xAOD", "All", Weight);


      SetFlag(FLAG_cut_temp,"Trigger","All", passTrig);
      if(SETTING[SETNAME[0]]["doweight"]==1)
        SetWeight(Evt_Weight, "Trigger", "All", Weight);


      int index_muon = 0;
      for(auto muon : *quickAna->muons()) {

        // without correction
        typedef ElementLink<xAOD::IParticleContainer> LinkType;
        static const char* linkName = "originalObjectLink";
        LinkType& auxLink = muon->auxdata<LinkType> (linkName);
        const xAOD::Muon* origMuon = dynamic_cast<const xAOD::Muon*>(*auxLink.cptr());  
        OBJ_MUON muonInfo;

        muonInfo.author = muon->author();
        muonInfo.charge = muon->charge();
        muonInfo.type =   muon->muonType();

        //muonInfo.quality = muon->quality();// why comment it --Zhang
     
        muonInfo.L.SetPtEtaPhiM( muon->pt(), muon->eta(), muon->phi(), m_mass ); 
        const xAOD::TrackParticle * trkPart = muon->primaryTrackParticle();
        muonInfo.d0 = trkPart->d0();
        muonInfo.z0 = trkPart->z0()+trkPart->vz()-pz0;
        muonInfo.d0err = 0;
        muonInfo.z0err = 0;
        if(muonInfo.type!=xAOD::Muon::MuonStandAlone) {
          const xAOD::ParametersCovMatrix_t TrkCovMatrix = trkPart->definingParametersCovMatrix();
          muonInfo.d0err = sqrt(TrkCovMatrix(0,0));
          muonInfo.z0err = sqrt(TrkCovMatrix(1,1));

          muonInfo.L_id.SetPtEtaPhiM(trkPart->pt(), trkPart->eta(), trkPart->phi(), m_mass);
          muonInfo.qoverp= trkPart->qOverP();
        }else muonInfo.L_id = muonInfo.L;

        muonInfo.d0sig = muonInfo.d0err!=0 ? fabs(muonInfo.d0/muonInfo.d0err) : 9999.;
        muonInfo.z0sig = muonInfo.z0err!=0 ? fabs(muonInfo.z0/muonInfo.z0err) : 9999.;
        

// isolation reconstruction is not good, kept for future

        float etcone20=0, ptcone20=0;
        if(muon->isolation(etcone20, xAOD::Iso::etcone20))
          muonInfo.topoetcone20 = etcone20; 
        if(muon->isolation(ptcone20, xAOD::Iso::ptcone20))
          muonInfo.ptcone20 = ptcone20; 

        muon->auxdata< char >( "All" ) = true;  // strange syntax ??? a function should be the left side of an assignment ??? --Zhang

        if(muon->auxdata<char> ("ana_select_hzhinv_loose_ID"))  muon->auxdata< char >( "Tool" ) = true;
        if(muon->auxdata<char> ("ana_select_hzhinv_loose_CB"))  muon->auxdata< char >( "CB" ) = true;
        if(muon->auxdata<char> ("ana_select_hzhinv_loose_Eta"))  muon->auxdata< char >( "Eta" ) = true;
        if(muon->auxdata<char> ("ana_select_hzhinv_loose_Pt"))  muon->auxdata< char >( "Pt" ) = true;
        if(muon->auxdata<char> ("ana_select_hzhinv_loose_D0"))  muon->auxdata< char >( "D0" ) = true;
        if(muon->auxdata<char> ("ana_select_hzhinv_loose_Z0"))  muon->auxdata< char >( "Z0" ) = true;
        if(muon->auxdata<char> ("ana_select_hzhinv_loose_Iso"))  muon->auxdata< char >( "TrkIso" ) = true;
        if(muon->auxdata<char> ("ana_select_hzhinv_loose")) muon->auxdata< char >( "OverLap" ) = true;

        bool passMuon=true;
	// in Initialize.h 55, 
	// STEP_obj["mu"] = {"All","Tool","CB","Eta","Pt","Z0","D0","TrkIso","OverLap","Pt15","Medium"}
        CountMuObj(muon, STEP_obj, CNT_obj, sysname, passMuon);
        if(passMuon) {
          muonInfo.trigM = muon->auxdata<bool>("HLT_mu20_iloose_L1MU15_OR_HLT_mu50_trigMatch");
          muonInfo.ismedium = (bool)muon->auxdata<char> ("ana_select_hzhinv_medium_ID");
          muonInfo.sf=muon->auxdata<float>("ana_weight_hzhinv_medium");
          //cout << "muon trig match is " << muonInfo.trigM << endl;
          muonInfo.index = index_muon;
          muon->auxdata<int>("index") = index_muon;
          temp_muon.push_back( muonInfo );
          index_muon++;
        }  // we should do passMuon judgement firstly, if denied, there is no need to record the corresponding muoninfo at all

      } // END OF TEMP_MUON SELECTION    

      for(int i=0; i<(int)temp_muon.size(); i++) {


        if(temp_muon[i].L.Pt()<7.e3) continue;
        DoCounting(sysname, CNT_obj, "mu", "Pt15");

        if(!temp_muon[i].ismedium) continue;
        DoCounting(sysname, CNT_obj, "mu", "Medium");

        goodm.push_back(temp_muon[i]);
      }

      int index_ele = 0;
      for(auto electron : *quickAna->electrons()) {
//        cout << "index ele " << index_ele <<endl;
        typedef ElementLink<xAOD::IParticleContainer> LinkType;
        static const char* linkName = "originalObjectLink";
        LinkType& auxLink = electron->auxdata<LinkType> (linkName);
        const xAOD::Electron* origEle = dynamic_cast<const xAOD::Electron*>(*auxLink.cptr());

        OBJ_ELECTRON eleInfo;

        eleInfo.author = electron->author();
        eleInfo.charge = electron->charge();

        eleInfo.L.SetPtEtaPhiM( electron->pt(), electron->eta(), electron->phi(), e_mass ); 

        const xAOD::TrackParticle * trkPart = electron->trackParticle();
        eleInfo.d0 = trkPart->d0();
        eleInfo.z0 = trkPart->z0()+trkPart->vz()-pz0;
        eleInfo.d0err = 0;
        eleInfo.z0err = 0;
        const xAOD::ParametersCovMatrix_t TrkCovMatrix = trkPart->definingParametersCovMatrix();
        eleInfo.d0err = sqrt(TrkCovMatrix(0,0));
        eleInfo.z0err = sqrt(TrkCovMatrix(1,1));
        eleInfo.d0sig = eleInfo.d0err!=0 ? fabs(eleInfo.d0/eleInfo.d0err) : 9999.;
        eleInfo.z0sig = eleInfo.z0err!=0 ? fabs(eleInfo.z0/eleInfo.z0err) : 9999.;

        eleInfo.trketa = trkPart->eta();
        eleInfo.trkphi = trkPart->phi();
        eleInfo.trkpt = trkPart->pt();
        eleInfo.qoverp= trkPart->qOverP();



        const xAOD::CaloCluster * cluster = electron->caloCluster();
        eleInfo.clE = cluster->e();
        eleInfo.clpt = cluster->pt();
        eleInfo.cleta = cluster->eta();
        eleInfo.clphi = cluster->phi();
        eleInfo.L_trk.SetPtEtaPhiE( trkPart->pt(), trkPart->eta(), trkPart->phi(), cluster->e() );

        float etcone20=0, ptcone20=0;
        if(electron->isolationValue(etcone20, xAOD::Iso::etcone20))
          eleInfo.etcone20 = etcone20;
        if(electron->isolationValue(ptcone20, xAOD::Iso::ptcone20))
          eleInfo.ptcone20 = ptcone20;

    
        //eleInfo.islikelihood = static_cast<bool> (myLikelihood->accept(origEle));
        eleInfo.passOQ = electron->isGoodOQ(xAOD::EgammaParameters::BADCLUSELECTRON);

        electron->auxdata< char >( "All" ) = true;

        electron->auxdata< char >( "ID" ) = true;
        if(electron->auxdata<char> ("ana_select_hzhinv_loose_Pt")) electron->auxdata< char >( "Pt" ) = true;
        if(electron->auxdata<char> ("ana_select_hzhinv_loose_Eta")) electron->auxdata< char >( "Eta" ) = true;
        if(electron->auxdata<char> ("ana_select_hzhinv_loose_OQ")) electron->auxdata< char >( "ObjQ" ) = true;
        if(electron->auxdata<char> ("ana_select_hzhinv_loose_selectionTool")) electron->auxdata< char >( "LHID" ) = true;
        if(electron->auxdata<char> ("ana_select_hzhinv_loose_Z0")) electron->auxdata< char >( "Z0" ) = true;
        if(electron->auxdata<char> ("ana_select_hzhinv_loose_D0")) electron->auxdata< char >( "D0" ) = true;
        if(electron->auxdata<char> ("ana_select_hzhinv_loose_Iso")) electron->auxdata< char >( "TrkIso" ) = true;
        if(electron->auxdata<char> ("ana_select_hzhinv_loose")) electron->auxdata< char >( "OverLap" ) = true;


        bool passEle=true;
	// STEP_obj["ele"]="All,ID,LHID,Pt,Eta,ObjQ,Z0,D0,TrkIso,OverLap,Pt15,Medium"
        CountEleObj(electron, STEP_obj, CNT_obj, sysname, passEle);

        if(passEle) {
          if(!isMC) eleInfo.trigM = electron->auxdata<bool>("HLT_e24_lhmedium_L1EM20VH_trigMatch") || electron->auxdata<bool>("HLT_e60_lhmedium_trigMatch") || electron->auxdata<bool>("HLT_e120_lhloose_trigMatch");
          else eleInfo.trigM = electron->auxdata<bool>("HLT_e24_lhmedium_L1EM18VH_trigMatch") || electron->auxdata<bool>("HLT_e60_lhmedium_trigMatch") || electron->auxdata<bool>("HLT_e120_lhloose_trigMatch");

          eleInfo.sf=electron->auxdata<float>("ana_weight_hzhinv_medium");

          if(eleInfo.L.Pt()>25.e3) // && eleInfo.trigM) 
          {
            eleInfo.trigSF=electron->auxdata<double>("HLT_e60_lhmedium_TrigSF");
            eleInfo.trigEFF=electron->auxdata<double>("HLT_e60_lhmedium_TrigEff");
          }
          else {
            eleInfo.trigSF = 0.0;
            eleInfo.trigEFF = 0.0;
          }    
 
          eleInfo.ismedium = (bool)electron->auxdata<char> ("ana_select_hzhinv_medium_selectionTool");
          eleInfo.index = index_ele;
          electron->auxdata<int>("index") = index_ele;
          temp_electron.push_back( eleInfo );
          index_ele++;
          //cout << "electron trig match is " << eleInfo.trigM << "  trg sf is " << electron->auxdata<double>("HLT_e60_lhmedium_TrigSF") <<" sec sf " << electron->auxdata<double>("HLT_e120_lhloose_TrigSF") << endl;
        }
      }

      for(int i=0; i<(int)temp_electron.size(); i++) {

        if(temp_electron[i].L.Pt()<7.e3) continue;
        DoCounting(sysname, CNT_obj, "ele", "Pt15");

        if(!temp_electron[i].ismedium) continue;
        DoCounting(sysname, CNT_obj, "ele", "Medium");

        goode.push_back(temp_electron[i]);
      }

      int index_jet = 0;
      bool passJetCleaning = true;
      bool bjettag = false;
      float total_jetsf = 1.0;
      int  nbjetOR = 0;
      for(auto jet : *quickAna->jets()) {

        typedef ElementLink<xAOD::IParticleContainer> LinkType;
        static const char* linkName = "originalObjectLink";
        LinkType& auxLink = jet->auxdata<LinkType> (linkName);
        const xAOD::Jet* origJet = dynamic_cast<const xAOD::Jet*>(*auxLink.cptr());

        OBJ_JET jetInfo;
	jetInfo.E = jet->e(); // in class xAOD::Jet_v1
        jetInfo.pt = jet->pt();
	jetInfo.pz = jet->pz(); // unsure --Zhang
        jetInfo.eta= jet->eta();
	jetInfo.phi = jet->phi(); 
        jetInfo.L  = jet->p4();

        jetInfo.numTrk0 = jet->getAttribute< std::vector<int> >(xAOD::JetAttribute::NumTrkPt500)[0];
        //cout << "#### bjet or not ####" << (bool)jet->auxdecor<char>("bjet") << endl;


        jet->auxdata< char >( "All" ) = true;

	// requires jet pt > 25GeV
        if( jet->pt()>25.e3 ) jet->auxdata< char >( "Pt" ) = true;
        if( fabs(jet->eta())<4.5 ) jet->auxdata< char >( "Eta" ) = true;
        if(jet->auxdata<char>("ana_select_jvt")) jet->auxdata< char >( "JVT" ) = true;
        if(jet->auxdata<char> ("ana_select")) jet->auxdata< char >( "OverLap") = true;

        bool passJet=true;
	// STEP_obj["jet"] = "All,Pt,Eta,JVT,OverLap,Clean"
        CountJetObj(jet, STEP_obj, CNT_obj, sysname, passJet);
        if(passJet) {
          if(jet->auxdecor<char>("bjet")) nbjetOR++;
          bjettag = bjettag || (bool)jet->auxdecor<char>("bjet");

          if(jet->auxdata<char>("ana_select_cleaning_tool")){
            DoCounting(sysname, CNT_obj, "jet", "Clean");
            jetInfo.index = index_jet;
            jetInfo.sf=jet->auxdata<float>("ana_weight");
            total_jetsf = total_jetsf*jetInfo.sf;
            jet->auxdata<int>("index") = index_jet;
            goodj.push_back( jetInfo );
            index_jet++;
          }
          else passJetCleaning=false;
        }
      }
//      if(goodj.size() != 2) return EL::StatusCode::SUCCESS;

      // only one missing Et
      xAOD::MissingET* met_rebuild = quickAna->met();
      double rmet = met_rebuild->met();
      double rmpx = met_rebuild->mpx();
      double rmpy = met_rebuild->mpy();
      //double rmphi= met_rebuild->phi();

      OBJ_MET met_obj;
      met_obj.L_RefFinal.SetPxPyPzE(rmpx, rmpy, 0, rmet);
      met_obj.L=met_obj.L_RefFinal;
      met_obj.met_rebuild = rmet;
      //met_obj.met = vmet;
      goodmet.push_back(met_obj);


      SetFlag(FLAG_cut_temp,"JetClean","All", passJetCleaning);
      if(SETTING[SETNAME[0]]["doweight"]==1)
        SetWeight(Evt_Weight, "JetClean", "All", Weight);

      int nMuons = goodm.size(); 
      int nElectrons = goode.size(); 
      int nJetOR = goodj.size(); 
      int nJet = temp_jet.size(); 

      vector<Pair> pairs;
      for(int i=0; i<(int)goodm.size(); i++) {
        for(int j=i+1; j<(int)goodm.size(); j++) {
          int index1, index2;
          if(goodm[i].L.Pt()>=goodm[j].L.Pt()) {index1=i; index2=j;}
          else {index1=j; index2=i;}  // index1.Pt > index2.Pt
          Pair temp;
          temp.flavor=1; 
          if(goodm[i].charge*goodm[j].charge == -1) temp.opposign=true;
          else temp.opposign=false;
          temp.index.push_back(index1); temp.index.push_back(index2);
          temp.lepton.push_back(goodm[index1].L); temp.lepton.push_back(goodm[index2].L);
          temp.lepton_trk.push_back(goodm[index1].L); temp.lepton_trk.push_back(goodm[index2].L);
          temp.charge.push_back((int)goodm[index1].charge); temp.charge.push_back((int)goodm[index2].charge);
          temp.Z = temp.lepton[0] + temp.lepton[1];
          temp.sf = goodm[index1].sf*goodm[index2].sf;
          temp.trigSF = evtInfo->auxdata<double>("HLT_mu20_iloose_L1MU15_OR_HLT_mu50_Mu_TrigSF");

          pairs.push_back(temp);
        }
      } 

      for(int i=0; i<(int)goode.size(); i++) {
        for(int j=i+1; j<(int)goode.size(); j++) {
          int index1, index2;
          if(goode[i].L.Pt()>=goode[j].L.Pt()) {index1=i; index2=j;}
          else {index1=j; index2=i;}
          Pair temp;
          temp.flavor=0;
          if(goode[i].charge*goode[j].charge == -1) temp.opposign=true;
          else temp.opposign=false;
          temp.index.push_back(index1); temp.index.push_back(index2);
          temp.lepton.push_back(goode[index1].L); temp.lepton.push_back(goode[index2].L);
          temp.lepton_trk.push_back(goode[index1].L_trk); temp.lepton_trk.push_back(goode[index2].L_trk);
          temp.charge.push_back((int)goode[index1].charge); temp.charge.push_back((int)goode[index2].charge);
          temp.Z = temp.lepton[0] + temp.lepton[1];
          temp.sf = goode[index1].sf*goode[index2].sf;
          float ineff_data = 1.0;
          float ineff_MC = 1.0;
	  // why only consider electron without muon 
          if (goode[index1].L.Pt()>25.e3){ // ????
            ineff_data *= (1-goode[index1].trigSF*goode[index1].trigEFF);
            ineff_MC   *= (1-goode[index1].trigEFF);
          }
          if (goode[index2].L.Pt()>25.e3){
            ineff_data *= (1-goode[index2].trigSF*goode[index2].trigEFF);
            ineff_MC   *= (1-goode[index2].trigEFF);
          }
	  // for float comparasion, it's better to use a range , rather than equal sign.
          if (ineff_MC==1)
            temp.trigSF = 0.0;
          else
            temp.trigSF = (1-ineff_data)/(1-ineff_MC);
          pairs.push_back(temp);
        }
      }

      for(vector<Pair>::const_iterator it = pairs.begin(); it != pairs.end(); ){
        if (!(*it).opposign) pairs.erase(it);
	else ++it;
      } // same flavour OS pairs

      vector<Pair> good_pairs;
      if(goode.size() == 2 && goodm.size() ==2){  // 2e2m
	for(int i=0; i<pairs.size(); i++)
	  good_pairs.push_back(pairs[i]);
      }
      else if(goode.size() == 4 || goodm.size() ==4){
        double tag_dmZ = 9999999999999999.e3;
	double dmZ = 999999999999999999.e3;
	int tag = -1;
	int index1 = -1, index2 = -1;
	for ( int i=0; i<pairs.size(); i++){
	  dmZ = fabs(pairs[i].Z.M() - ZMass);
	  if(dmZ < tag_dmZ) {
	    tag_dmZ = dmZ;
	    tag = i;
	  }
	}
	if(tag_dmZ < 25.e3) {
	  good_pairs.push_back(pairs[tag]);
	  index1 = pairs[tag].index[0];
	  index2 = pairs[tag].index[1];
	  vector<Pair>::const_iterator it = pairs.begin();
	  for(; it != pairs.end(); ) {
	    int  index = (*it).index[0];
	    if(index == index1 || index == index2)
	      pairs.erase(it);
	    else{
	      index = (*it).index[1];
	      if(index == index1 || index == index2)
		pairs.erase(it);
	      else ++it;
	    }
	  }
	  if(pairs.size() == 1) good_pairs.push_back(pairs[0]);
	}
      } // good_pairs

      double Pt1=-9999.0, Pt2=-9999.0, Pt3=-9999.0,Pt4=-9999.0;
      double Eta1=-9999.0, Eta2=-9999.0, Eta3=-9999.0, Eta4=-9999.0;
      double Phi1=-9999.0, Phi2=-9999.0, Phi3=-9999.0, Phi4=-9999.0;
      double DiLepton_Mass1=-9999.0, DiLepton_Mass2=-9999.0;
      double DiLepton_Pt1=-9999.0, DiLepton_Pt2=-9999.0;
      double DiLepton_Rap1=-9999.0, DiLepton_Rap2=-9999.0; // DiLep rapidity
      double DiJet_Mass=-9999.0;
      int LeadingJetTag=-1, SubLeadingJetTag=-1;
      double LeadingJetPt=-9999.0, SubLeadingJetPt=-9999.0;
      double LeadingJetEta=-9999.0, SubLeadingJetEta=-9999.0;
      double LeadingJetPhi=-9999.0, SubLeadingJetPhi=-9999.0;
      double LeadingJetRap=-9999.0, SubLeadingJetRap=-9999.0;
      double DeltaJetEta=-9999.0; // △(η) = |η1 - η2|
      double Centrality=-9999.0; // C=[y(4l) - (y(j1) + y(j2))/2]/|y(j1)-y(j2)|
      double _4Lepton_Mass=-9999.0;
      double _4Lepton_Pt=-9999.0;
      double _4Lepton_Rap=-9999.0;
      int Number_Jets=-1;
      bool Pass4Leptons_mmmm = false;
      bool Pass4Leptons_eeee = false;
      bool Pass4Leptons_eemm = false;
      bool Pass4Leptons = false;
      bool Pass2Pairs = false;
      bool PassMll1 = false;
      bool PassMll2 = false;
      bool PassJet = false;
      double weight = 1.0;
      double lep_sf = 1.0, trigger_sf = 1.0;

      if(goodm.size() == 4 && goode.size() == 0){  // 4 muons channel
	Pass4Leptons_mmmm = true;
	lep_sf = goodm[0].sf*goodm[1].sf*goodm[2].sf*goodm[3].sf;
	trigger_sf = evtInfo->auxdata<double>("HLT_mu20_iloose_L1MU15_OR_HLT_mu50_Mu_TrigSF");

	Pt1 = goodm[0].L.Pt();
	Pt2 = goodm[1].L.Pt();
	Pt3 = goodm[2].L.Pt();
	Pt4 = goodm[3].L.Pt();
	Eta1 = goodm[0].L.Eta();
	Eta2 = goodm[1].L.Eta();
	Eta3 = goodm[2].L.Eta();
	Eta4 = goodm[3].L.Eta();
	Phi1 = goodm[0].L.Phi();
	Phi2 = goodm[1].L.Phi();
	Phi3 = goodm[2].L.Phi();
	Phi4 = goodm[3].L.Phi();
      }else if (goode.size() == 4 && goodm.size() == 0){  // 4 electrons channel
	Pass4Leptons_eeee = true;
	lep_sf = goode[0].sf*goode[1].sf*goode[2].sf*goode[3].sf;
	float ineff_data=1.0, ineff_MC=1.0;
	for(int i=0; i<4; i++){
	  if (goode[i].L.Pt() > 25.e3){
	    ineff_data *= (1-goode[i].trigSF*goode[i].trigEFF);
	    ineff_MC *= (1-goode[i].trigEFF);
	  }
	}
	if (ineff_MC==1) trigger_sf = 0.0;
	else trigger_sf = (1-ineff_data)/(1-ineff_MC);
	

	Pt1 = goode[0].L.Pt();
	Pt2 = goode[1].L.Pt();
	Pt3 = goode[2].L.Pt();
	Pt4 = goode[3].L.Pt();
	Eta1 = goode[0].L.Eta();
	Eta2 = goode[1].L.Eta();
	Eta3 = goode[2].L.Eta();
	Eta4 = goode[3].L.Eta();
	Phi1 = goode[0].L.Phi();
	Phi2 = goode[1].L.Phi();
	Phi3 = goode[2].L.Phi();
	Phi4 = goode[3].L.Phi();
      }else if( goodm.size() == 2 && goode.size() == 2){
	Pass4Leptons_eemm = true;
	lep_sf = goodm[0].sf*goodm[1].sf*goode[0].sf*goode[1].sf;
	float ineff_data=1.0, ineff_MC=1.0;
	trigger_sf = evtInfo->auxdata<double>("HLT_mu20_iloose_L1MU15_OR_HLT_mu50_Mu_TrigSF");
	
        if (goode[0].L.Pt()>25.e3){
          ineff_data *= (1-goode[0].trigSF*goode[0].trigEFF);
          ineff_MC   *= (1-goode[0].trigEFF);
        }
        if (goode[1].L.Pt()>25.e3){
          ineff_data *= (1-goode[1].trigSF*goode[1].trigEFF);
          ineff_MC   *= (1-goode[1].trigEFF);
        }
        if (ineff_MC==1)
          trigger_sf = 0.0;
        else
          trigger_sf *= (1-ineff_data)/(1-ineff_MC);
    

	Pt1 = goodm[0].L.Pt();
	Pt2 = goodm[1].L.Pt();
	Pt3 = goode[0].L.Pt();
	Pt4 = goode[1].L.Pt();
	Eta1 = goodm[0].L.Eta();
	Eta2 = goodm[1].L.Eta();
	Eta3 = goode[0].L.Eta();
	Eta4 = goode[1].L.Eta();
	Phi1 = goodm[0].L.Phi();
	Phi2 = goodm[1].L.Phi();
	Phi3 = goode[0].L.Phi();
	Phi4 = goode[1].L.Phi();
      } // extract Pt, Eta and Phi value from 4 leptons

      // sort Pt in descending order
      { int index[4]={0,1,2,3};
	double Pt[4]={Pt1,Pt2,Pt3,Pt4};
	double Eta[4]={Eta1,Eta2,Eta3,Eta4};
	double Phi[4]={Phi1,Phi2,Phi3,Phi4};
	for( int i=0; i<4; i++ ){
	  for( int j=i+1; j<4; j++) {
	    if (Pt[index[j]] > Pt[index[i]]) {
	      int temp = index[i];
	      index[i] = index[j];
	      index[j] = temp;
	    }
	  }
	}
	Pt1 = Pt[index[0]];
	Pt2 = Pt[index[1]];
	Pt3 = Pt[index[2]];
	Pt4 = Pt[index[3]];
	Eta1 = Eta[index[0]];
	Eta2 = Eta[index[1]];
	Eta3 = Eta[index[2]];
	Eta4 = Eta[index[3]];
	Phi1 = Phi[index[0]];
	Phi2 = Phi[index[1]];
	Phi3 = Phi[index[2]];
	Phi4 = Phi[index[3]];
      }

      if(good_pairs.size() == 2){
	Pass2Pairs = true;
        DiLepton_Mass1 = good_pairs[0].Z.M();
        DiLepton_Pt1 = good_pairs[0].Z.Pt();
        DiLepton_Rap1 = good_pairs[0].Z.Rapidity();
        DiLepton_Mass2 = good_pairs[1].Z.M();
        DiLepton_Pt2 = good_pairs[1].Z.Pt();
        DiLepton_Rap2 = good_pairs[1].Z.Rapidity();

	TLorentzVector Total = good_pairs[0].Z + good_pairs[1].Z;
	_4Lepton_Mass = Total.M();
	_4Lepton_Pt = Total.Pt();
	_4Lepton_Rap = Total.Rapidity();
      }

      Number_Jets = goodj.size();
      if(Number_Jets > 1){
	if(goodj[0].pt > goodj[1].pt) {
	  LeadingJetTag = 0;
	  SubLeadingJetTag = 1;
	}
	else { LeadingJetTag = 1; SubLeadingJetTag = 0; }
        for(int i=2; i<Number_Jets; i++){
	  if(goodj[i].pt > goodj[LeadingJetTag].pt){
	    SubLeadingJetTag = LeadingJetTag;
	    LeadingJetTag = i;
	  }    
	  else if(goodj[i].pt > goodj[SubLeadingJetTag].pt) SubLeadingJetTag = i;
	}
	LeadingJetPt = goodj[LeadingJetTag].pt;
	LeadingJetEta = goodj[LeadingJetTag].eta;
	LeadingJetPhi = goodj[LeadingJetTag].phi;
	LeadingJetRap = goodj[LeadingJetTag].L.Rapidity();
	SubLeadingJetPt = goodj[SubLeadingJetTag].pt;
	SubLeadingJetEta = goodj[SubLeadingJetTag].eta;
	SubLeadingJetPhi = goodj[SubLeadingJetTag].phi;
	SubLeadingJetRap = goodj[SubLeadingJetTag].L.Rapidity();
	DeltaJetEta = LeadingJetEta - SubLeadingJetEta;
	DiJet_Mass = (goodj[LeadingJetTag].L + goodj[SubLeadingJetTag].L).M();

	Centrality = ( _4Lepton_Rap - (LeadingJetRap + SubLeadingJetRap)/2)/fabs(LeadingJetRap - SubLeadingJetRap);
      }

      HistVar["Pt1"]["Value"] = Pt1;
      HistVar["Pt2"]["Value"] = Pt2;
      HistVar["Pt3"]["Value"] = Pt3;
      HistVar["Pt4"]["Value"] = Pt4;
      HistVar["Eta1"]["Value"] = Eta1;
      HistVar["Eta2"]["Value"] = Eta2;
      HistVar["Eta3"]["Value"] = Eta3;
      HistVar["Eta4"]["Value"] = Eta4;
      HistVar["Phi1"]["Value"] = Phi1;
      HistVar["Phi2"]["Value"] = Phi2;
      HistVar["Phi3"]["Value"] = Phi3;
      HistVar["Phi4"]["Value"] = Phi4;
      HistVar["NumJets"]["Value"] = Number_Jets;
      HistVar["LeadingJetPt"]["Value"] = LeadingJetPt;
      HistVar["LeadingJetEta"]["Value"] = LeadingJetEta;
      HistVar["LeadingJetPhi"]["Value"] = LeadingJetPhi;
      HistVar["SubLeadingJetPt"]["Value"] = SubLeadingJetPt;
      HistVar["SubLeadingJetEta"]["Value"] = SubLeadingJetEta;
      HistVar["SubLeadingJetPhi"]["Value"] = SubLeadingJetPhi;
      HistVar["DiJetMass"]["Value"] = DiJet_Mass;
      HistVar["DeltaJetEta"]["Value"] = DeltaJetEta;
      HistVar["DiLepMass1"]["Value"] = DiLepton_Mass1;
      HistVar["DiLepMass2"]["Value"] = DiLepton_Mass2;
      HistVar["DiLepPt1"]["Value"] = DiLepton_Pt1;
      HistVar["DiLepPt2"]["Value"] = DiLepton_Pt2;
      HistVar["DiLepRap1"]["Value"] = DiLepton_Rap1;
      HistVar["DiLepRap2"]["Value"] = DiLepton_Rap2;
      HistVar["4LepMass"]["Value"] = _4Lepton_Mass;
      HistVar["4LepPt"]["Value"] = _4Lepton_Pt;
      HistVar["Centrality"]["Value"] = Centrality;

//    HistVar["truth_lep_num"]["Value"] = truth_lep_num;
//    HistVar["truth_ele_num"]["Value"] = truth_ele_num;
//    HistVar["truth_muon_num"]["Value"] = truth_muon_num;
//    HistVar["truth_Pt1"]["Value"] = truth_Pt1;
//    HistVar["truth_Pt2"]["Value"] = truth_Pt2;
//    HistVar["truth_Pt3"]["Value"] = truth_Pt3;
//    HistVar["truth_Pt4"]["Value"] = truth_Pt4;
//    HistVar["truth_Eta1"]["Value"] = truth_Eta1;
//    HistVar["truth_Eta2"]["Value"] = truth_Eta2;
//    HistVar["truth_Eta3"]["Value"] = truth_Eta3;
//    HistVar["truth_Eta4"]["Value"] = truth_Eta4;
//    HistVar["truth_Phi1"]["Value"] = truth_Phi1;
//    HistVar["truth_Phi2"]["Value"] = truth_Phi2;
//    HistVar["truth_Phi3"]["Value"] = truth_Phi3;
//    HistVar["truth_Phi4"]["Value"] = truth_Phi4;
      // Setflags
      
      // set the channel
      string CHANNEL="All";
      if(Pass4Leptons_mmmm) CHANNEL = "mmmm";
      else if(Pass4Leptons_eemm) CHANNEL = "eemm";
	   else if(Pass4Leptons_eeee) CHANNEL = "eeee";
      if (Pass4Leptons_eeee || Pass4Leptons_eemm || Pass4Leptons_mmmm ){
	Pass4Leptons =  true;
	SetFlag(FLAG_cut_temp, "4Leptons", CHANNEL, Pass4Leptons);
        if(SETTING[SETNAME[0]]["doweight"] == 1)
	    SetWeight(Evt_Weight, "4Leptons", "All", Weight);
      }

      if (good_pairs.size() == 2){
	Pass2Pairs = true;
	Weight = Weight*lep_sf*trigger_sf;
	SetFlag(FLAG_cut_temp, "2Pairs", CHANNEL, Pass2Pairs);
        if(SETTING[SETNAME[0]]["doweight"] == 1)
	    SetWeight(Evt_Weight, "2Pairs", CHANNEL, Weight);
      }

      if (fabs(DiLepton_Mass1 - ZMass) < 25.e3){
	PassMll1 = true;
	SetFlag(FLAG_cut_temp, "Mll1", CHANNEL, 1); // one on-shell Z
        if(SETTING[SETNAME[0]]["doweight"] == 1)
	    SetWeight(Evt_Weight, "Mll1", CHANNEL, Weight);
      }
      if(fabs(DiLepton_Mass2 -  ZMass) < 25.e3){
	PassMll2 = true;
	SetFlag(FLAG_cut_temp, "Mll2", CHANNEL, 1); // two on-shell Z
        if(SETTING[SETNAME[0]]["doweight"] == 1)
	    SetWeight(Evt_Weight, "Mll2", CHANNEL, Weight);
      }

      if( PassMll2 && Number_Jets >= 2 ){
	PassJet = true;
	Weight *= total_jetsf;
	SetFlag(FLAG_cut_temp,"NumJets", CHANNEL, 1);
        if(SETTING[SETNAME[0]]["doweight"] == 1)
            SetWeight(Evt_Weight, "NumJets", CHANNEL, Weight);
      }
      
      
      if(SETTING[SETNAME[0]]["doopt"]==1)
        CountEvt(sysname, CHN, STEP_cut, FLAG_cut_temp, FLAG_cut, CNT_cut, Evt_Weight, Var_Map);
      else
        CountEvt(sysname, CHN, STEP_cut, FLAG_cut_temp, FLAG_cut, CNT_cut, Evt_Weight);      
 
      FillHistograms(sysname);

        TreeStrVar["filename"]["Value"] = m_FileName;
        TreeIntVar["signal_bkg"]["Value"] = SB::_signal;
        TreeIntVar["run"]["Value"] = run;
        TreeUloVar["event"]["Value"] = event;
        TreeIntVar["nbjetOR"]["Value"] = nbjetOR;
        TreeFltVar["met"]["Value"] = rmet*0.001;
        TreeFltVar["met_x"]["Value"] = rmpx*0.001;
        TreeFltVar["met_y"]["Value"] = rmpy*0.001;
        TreeFltVar["pileupweight"]["Value"] = pileWeight;
        TreeFltVar["mcweight"]["Value"] = mcWeight;
	TreeFltVar["lep_sf"]["Value"] = lep_sf;
	TreeFltVar["trigsf"]["Value"] = trigger_sf;
        Tree[sysname]->Fill();

      if(Pass4Leptons){
	m_4leptons++;
	if(Pass2Pairs){
	  m_2pairs++;
	  if(PassMll1){
	    m_mll1++;
	    if(PassMll2){
	      m_mll2++;
	      if(PassJet){
		m_numjets++;
	      }
	    }
	  }
	}
      }
    }
    m_sumOfWeights += Weight;
  }
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: finalize ()
{

  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.

  delete quickAna; 

  system("rm -rf cutflow"); 
  system("mkdir -p cutflow"); 

  //<Cut Flow Output>
  for(int i=0; i<(int)CHN.size(); i++) {
    string chn=CHN[i];
    string filename = "log_eff_" + chn + "_physics.txt";
    ofstream file(filename.c_str());

    if(file.is_open()) {
	// according to the output of this function, we will find that
	// this for loop only want to output one particular "sys",
	// then we can output it explicitly without the engagement of 
	// for loop
      for(int n=0; n<(int)SYSNAME.size(); n++) {
        string sys = SYSNAME[n];

        if(sys!="NOCORR" && (!SETTING["physics"]["docorr"])) continue;

        if(sys!="NOMINAL" && SETTING["physics"]["docorr"] && (!SETTING["physics"]["dosys"])) continue;

        if(sys=="NOCORR" && SETTING["physics"]["docorr"] && SETTING["physics"]["dosys"]) continue;


        file << "### <Channel : " << chn << "> ; Systematic type is : " 
             << sys << " ###" << endl;
        file <<endl;
        file << "WeightSum = " << sumOfWeights << endl;

	cout << "### <Channel : " << chn << "> ; Systematic type is :  " << sys << " ###" << endl;
	cout << endl;
	cout << "WeightSum = " << sumOfWeights << endl;

        if(SETTING[SETNAME[0]]["doopt"]==1) {
          for(int n=0; n<(int)Var_Map.size(); n++) {
            string fullvar=Var_Map[n];

            file << endl;
            file << "~~~~~ variation : " << fullvar << "  ~~~~~" << endl;
            file << endl;

            cout << endl;
            cout << "~~~~~ variation : " << fullvar << "  ~~~~~" << endl;
            cout << endl;

            for(int j=0; j<(int)STEP_cut.size(); j++) {
              string cut = STEP_cut[j];
              file<<cut<<" = "<<CNT_cut[sys][chn][cut][fullvar].num
                  << " ; " << CNT_cut[sys][chn][cut][fullvar].wnum <<
                  " +/- "<<CNT_cut[sys][chn][cut][fullvar].err<<endl;

              cout <<cut<<" = "<<CNT_cut[sys][chn][cut][fullvar].num
                  << " ; " << CNT_cut[sys][chn][cut][fullvar].wnum <<
                  " +/- "<<CNT_cut[sys][chn][cut][fullvar].err<<endl;
            }
            file << endl;
            file << " ======== border ======== " << endl;
            file << endl;

            cout << endl;
            cout << " ======== border ======== " << endl;
            cout << endl;
          }
        }else {
          for(int j=0; j<(int)STEP_cut.size(); j++) {
            string cut = STEP_cut[j];
            file<<cut<<" = "<<CNT_cut[sys][chn][cut]["default"].num
                << " ; " << CNT_cut[sys][chn][cut]["default"].wnum <<
                " +/- "<<CNT_cut[sys][chn][cut]["default"].err<<endl;

            cout <<cut<<" = "<<CNT_cut[sys][chn][cut]["default"].num
                << " ; " << CNT_cut[sys][chn][cut]["default"].wnum <<
                " +/- "<<CNT_cut[sys][chn][cut]["default"].err<<endl;
          }
	  // this for loop is provided in both consequential sentenses
	  // so it can be extracted from the if conditions
          file << endl;
          file << " ======== border ======== " << endl;
          file << endl;

          cout << endl;
          cout << " ======== border ======== " << endl;
          cout << endl;
        }
      }
    }
    else {
      cout<<"Can not open file "<<filename<<endl;
      exit(-1);
    }
    file.close();

    TMacro* m = new TMacro(filename.c_str());
    m->Write();
    delete m;

    string command1 = "mv " + filename + " cutflow";
    system(command1.c_str());
    // isn't there any way to write file directly in a dir.
  }

  string filename_obj = "log_eff_obj_physics.txt";
  ofstream file_obj(filename_obj.c_str());
  for(int n=0; n<(int)SYSNAME.size(); n++) {
    string sys = SYSNAME[n];

    if(sys!="NOCORR" && (!SETTING["physics"]["docorr"])) continue;

    if(sys!="NOMINAL" && SETTING["physics"]["docorr"] && (!SETTING["physics"]["dosys"])) continue;

    if(sys=="NOCORR" && SETTING["physics"]["docorr"] && SETTING["physics"]["dosys"]) continue;


    file_obj << "### <Object Selection> Systematic type is : "
	     << sys << " ###" << endl;
    file_obj <<endl;

    cout << "### <Object Selection> Systematic type is : "
	     << sys << " ###" << endl;
    cout << endl;

    MapType_VString::iterator it;
    for(it=STEP_obj.begin(); it!=STEP_obj.end(); it++) {
      string obj=(*it).first;
      for(int i=0; i<(int)STEP_obj[obj].size(); i++) {
        string cut=STEP_obj[obj][i];
        string cutname=obj+"_"+cut;
        file_obj<<cutname<<" = "<<CNT_obj[sys][obj][cut].num
                <<" +/- "<<CNT_obj[sys][obj][cut].err<<endl;

        cout <<cutname<<" = "<<CNT_obj[sys][obj][cut].num
                <<" +/- "<<CNT_obj[sys][obj][cut].err<<endl;
      }
      file_obj<< endl;
      cout << endl;
    }
    file_obj << endl;
    file_obj << " ======== border ======== " << endl;
    file_obj << endl;

    cout << endl;
    cout << " ======== border ======== " << endl;
    cout << endl;
  }  
  file_obj.close();


  TMacro* m = new TMacro(filename_obj.c_str());
  m->Write();
  delete m;

  string command2= "mv " + filename_obj + " cutflow";
  system(command2.c_str());

  string truth_cutflow = "truth_cutflow";
  ofstream file(truth_cutflow.c_str());
  file << "all: " << m_truth_xAOD << endl;
  file << "4Leptons: " << m_truth_4leptons << endl;
  file << "LepPt: " << m_truth_leppt << endl;
  file << "LepEta: " << m_truth_lepeta << endl;
  file << "LeadingLepPt: " << m_truth_leadingleppt << endl;
  file << "2Pairs: " << m_truth_2pairs << endl;
  file << "Mll1: " << m_truth_mll1 << endl;
  file << "Mll2: " << m_truth_mll2 << endl;
  file << "NumJets: " << m_truth_numjets << endl;
  file << "===== border =====" << endl;
  file << "all: " << w_truth_xAOD << endl;
  file << "4Leptons: " << w_truth_4leptons << endl;
  file << "LepPt: " << w_truth_leppt << endl;
  file << "LepEta: " << w_truth_lepeta << endl;
  file << "LeadingLepPt: " << w_truth_leadingleppt << endl;
  file << "2Pairs: " << w_truth_2pairs << endl;
  file << "Mll1: " << w_truth_mll1 << endl;
  file << "Mll2: " << w_truth_mll2 << endl;
  file << "NumJets: " << w_truth_numjets << endl;
  file << "===== END =====" << endl;
  file.close();

  cout << endl;
  cout << "####" << endl;
  printf("Finalize : %i events have been processed !\n", m_eventCounter);
  printf("Finalize : %f weights have been processed !\n", m_filter);
  printf("Finalize MyxAODAnalysis !\n");
  cout << "####" << endl;
  cout << endl;

  time(&end);  // time(&start) in initial() file
  double timedif = difftime(end,start);
  if(timedif>3600) { printf("Finalize : Time Cost: %f hours\n", timedif/3600.); }
  else if(timedif>60) { printf("Finalize : Time Cost: %f minutes\n", timedif/60.); }
  else { printf("Finalize : Time Cost: %f second\n", timedif); }


  cout << "Events pass 4 leptons check: " << m_4leptons << endl;
  cout << "Events pass 2 pairs check: " << m_2pairs << endl;
  cout << "Events pass 1 M.Z check: " << m_mll1 << endl;
  cout << "Events pass 2 M.Z check: " << m_mll2 << endl;
  cout << "Events pass 2 jets check: " << m_numjets << endl;
  cout << "SumofWeight: " << m_sumOfWeights << endl;

  cout << "all event: " << m_truth_xAOD << endl;
  cout << "Truth Events pass 4 leptons check: " << m_truth_4leptons << endl;
  cout << "Truth Events pass LepPt cut : " << m_truth_leppt << endl;
  cout << "Truth Events pass LepEta cut : " << m_truth_lepeta << endl;
  cout << "Truth Events pass leadingLepPt cut : " << m_truth_leadingleppt << endl;
  cout << "Truth Events pass 2 pairs check: " << m_truth_2pairs << endl;
  cout << "Truth Events pass 1 M.Z check: " << m_truth_mll1 << endl;
  cout << "Truth Events pass 2 M.Z check: " << m_truth_mll2 << endl;
  cout << "Truth Events pass 2 jets check: " << m_truth_numjets << endl;

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODAnalysis :: histFinalize ()
{


  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.
  return EL::StatusCode::SUCCESS;
}
