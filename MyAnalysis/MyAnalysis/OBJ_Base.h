#ifndef MyAnalysis_OBJ_H
#define MyAnalysis_OBJ_H

//<Common Header>
#include <TLorentzVector.h>
#include <string>
#include <vector>
using namespace std;

//<Predefined Constants>
const double PI=3.1415926;
const double ZMass=91.1876e3; //MeV
//const double m_mass=105.658367; //MeV
//const double e_mass = 0.510998910; //MeV
const double m_mass=105.658; //MeV
const double e_mass = 0.511; //MeV
const double Unit_GeV = 1000.; //MeV

// Truth Object
// Truth_Electron, Truth_Muon
class TRUTH_OBJ
{
    public:
	double E, Et;
	double p, pt, px, py, pz;
	double eta, phi, m;
	int charge;
};

class OBJ_TRUTH_MUON : public TRUTH_OBJ 
{
    public:
	double id;
};

class OBJ_TRUTH_ELECTRON : public TRUTH_OBJ
{
    public:
	double id;
};

class OBJ_TRUTH_JET : public TRUTH_OBJ
{
    public:
	double id;
};

class Truth_Pair
{
    public:
	int flavour; // muon: 1; electron: 2;
	bool opposign;
	double mass;
	vector<int> index;
};


//<Object Class>
//<General, Muon, Electron, Jet, MET>
class OBJ{
    public:
        double E,Et;
        double p,pt,px,py,pz;
        double eta,phi,m;
        double d0,z0,d0err,z0err,d0sig,z0sig;
        double d0_old,z0_old,d0err_old,z0err_old,d0sig_old,z0sig_old;
        float sf,trigSF,trigEFF;
        int index;
        bool trigM;
        TLorentzVector L,Lori,Lcorr,LF;       
};    

class OBJ_MUON : public OBJ{
    public:
        int author, charge, type, quality; // type is added by Cong
        double iscombined, isloose, ismedium, istight, istag, issa; // classified into type and quality
        double id; /* 1: staco, 2: ST, 3: SA, 4: calo */
        double caloId, caloLikelihood;
        int passhits_staco, passhits_calo, passhits_sa, passhits_3rdChain;
        int passhits_blayer, passhits_pix, passhits_sct, passhits_hole, passhits_trt;
        int isgood, iscr; //too define control region. isgood==1 for std good muon
        double ptme, ptid, ptms, pt_corr, ptme_corr, ptid_corr,qoverp;
        double etcone20, etcone30, etcone40, ptcone20, ptcone30, ptcone40, etcone20_corr, ptcone20_corr;
        double topoetcone20;
        double mu_calo_eta, mu_calo_phi, mu_staco_eta, mu_staco_phi;
        double met_corrx, met_corry;
        double sf_loose, sf_calo;
        double npv;
        vector<int> fsr_index;
        vector<double> fsr_dR;
        vector<string> fsr_type;
        TLorentzVector L_me, L_ms, L_id, Lcorr_me, Lcorr_ms, Lcorr_id, LF_me, LF_id;
};

class OBJ_ELECTRON: public OBJ{
    public:
        int source; //1: from Z, 2: from b/c, 3: from hadron, 4: from photon conversion, 5: unknown, 6: from W
        int iscr, isgood; //too define control region. isgood==1 for std good electron, iscr==1 for cr
        double charge;
        uint16_t author;
        double isloose, ismedium, istight, isloosepp, ismediumpp, istightpp, islooseppzz, ismultilep, islikelihood;
        bool passOQ;
        double clE, clpt, cleta, clphi, trkpt, trketa, trkphi, qoverp;
        double trkpt2, trketa2, trkphi2;
        double clE_corr, pt2_corr, clE_noEP;
        double etcone20, etcone20_corr, etcone30, etcone40, ptcone20, ptcone30, ptcone40, ptcone20_corr;
        double el_etas2, el_phis2, el_rawcl_E;
        double met_corrx, met_corry;
        double sf_loosepp, sf_mediumpp, sf_tightpp, sf_id, sf_reco;
        int npv;
        TLorentzVector Lcl,L_trk,LF_trk,LF_tri; 
};       

class OBJ_JET : public OBJ{
    public:
        double pt_em, eta_em, phi_em, m_em;
        double isbad, jes, isugly, islooserbad;
        double vtxf, fmax, smax;
        int numTrk0;
        TLorentzVector LEM;
};

class OBJ_MET : public OBJ{
    public:
        double pt_corr,px_corr,py_corr;
        double met, met_rebuild;
        TLorentzVector L_RefFinal,L_Topo;
        TLorentzVector L_RefFinal_Corr,L_Topo_Corr;
};

//class for lepton pairs
class Pair {
    public:
        int flavor; //0-electron, 1-muon
        vector<int> index;
        vector<int> charge;
        vector<TLorentzVector> lepton;
        vector<TLorentzVector> lepton_trk;
        TLorentzVector Z;
        double mass;
        float sf,trigSF,trigEFF;
        bool opposign; //charge*charge of two leptons, if -1 then true, else false

};

//class for quads, i.e. Z pairs
class Quad {
    public:
        vector<int> index;
        vector<Pair> pair;
        vector<TLorentzVector> alter_pairs;
        vector<TLorentzVector> lepton;
        vector<TLorentzVector> lepton_trk;
        vector<int> index_lep;
        vector<int> flavor_lep;
        vector<int> charge_lep;
        vector<double> ptcone20, etcone20, ptcone20_corr, etcone20_corr, d0sig, z0sig, d0, z0;
        TLorentzVector ZZ;
};


typedef vector<OBJ_MUON> VOmuon;
typedef vector<OBJ_ELECTRON> VOelectron;
typedef vector<OBJ_JET> VOjet;
typedef vector<OBJ_MET> VOmet;
typedef vector<OBJ_TRUTH_ELECTRON> VOtruth_electron;
typedef vector<OBJ_TRUTH_MUON> VOtruth_muon;



#endif
