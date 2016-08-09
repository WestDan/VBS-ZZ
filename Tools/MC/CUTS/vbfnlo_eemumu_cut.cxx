#include <iostream>
#include <string>
#include <time.h>
using namespace std;

// Is there any way to directly operate trees, rather than copy
// part of one tree to another.

// this program only suitable for those 2e2mu events.
//
//	 what this program do 
// tree =====================> cutted tree
// cut events, using macro for all cut conditions.
typedef map<string, bool> MapType1_Bool;
typedef map<string, map<string, double> > MapType2_Double;

time_t start, finish;

MapType1_Bool FLAG_cut;

int convert(Double_t px, Double_t py, Double_t pz, Double_t & p, Double_t & pt, Double_t & eta, Double_t & phi)
{
    // θ ∈ [0, π], ψ ∈ [-π, π], 
    // px = p*sin(θ)*cos(ψ), py = p*sin(θ)*sin(ψ), pz = p*cos(θ)
    // η = -ln( tan(θ/2) ) = -log( sin(θ)/( 1 + cos(θ) ) ) 
    // log() base on e, log10() base on 10
    // sin(θ) = pt/p
    p = sqrt(px*px + py*py + pz*pz);
    pt = sqrt(px*px + py*py);
    double sin_theta = pt/p;
    double cos_theta = pz/p;
    eta = -log(sin_theta/( 1+cos_theta) );
    // the return value range of atan() is [-π, π]
    phi = atan(py/px);
    if(abs(pt*cos(phi) - px) < 1 ) ;
    else 
      if( phi > 0) phi -= TMath::Pi();
      else phi += TMath::Pi();
    return 0;
}

void vbfnlo_eemumu_cut() 
{

    // input root file
    string input_root = "output_tree_1092.root";
    TFile * file = new TFile(input_root.c_str(), "read");
    if( ! file ) 
    {
	cout << "Error, file " << input_root << "doesn't exit" << endl;
	return ;
    }

    // input tree
    string input_tree_name = "vbfnlo";
    TTree * input_tree = (TTree *) file->Get(input_tree_name.c_str());
    if ( ! input_tree ) 
    {
	cout << "input_tree " << input_tree_name << " doesn't exist in file " << input_root << endl;
	return ;
    }

    // output root file
    string output_root_name = "output.root";
    TFile * output_root = new TFile(output_root_name.c_str(), "recreate");
    if( ! output_root )
    {
	cout << "Error! unable to create file: output.root" << endl;
	return;
    }

    Double_t px_e[2], px_m[2];
    Double_t py_e[2], py_m[2];
    Double_t pz_e[2], pz_m[2];
    Double_t p_e[2], pt_e[2], p_m[2], pt_m[2];
    Double_t eta_e[2], eta_m[2];
    Double_t phi_e[2], phi_m[2];
    Double_t E_e[2], E_m[2];
    Double_t M_e[2], M_m[2];
    Double_t p_ll[2], pt_ll[2], eta_ll[2], phi_ll[2];
    Double_t E_ll[2], M_ll[2];
    Double_t px_llll, py_llll, pz_llll;
    Double_t p_llll, pt_llll;
    Double_t E_llll, M_llll;
    Double_t eta_llll, phi_llll;
    Double_t px_j[2];   // Jets
    Double_t py_j[2];
    Double_t pz_j[2];
    Double_t p_j[2], pt_j[2];
    Double_t E_j[2];
    Double_t M_j[2];
    Double_t eta_j[2], phi_j[2];
    Double_t px_jj, py_jj, pz_jj;
    Double_t p_jj, pt_jj;
    Double_t eta_jj, phi_jj;
    Double_t E_jj, M_jj;
    Double_t deltaEta_jj;
    Double_t centrality;

    Double_t px[10];
    Double_t py[10];
    Double_t pz[10];
    Double_t E[10];
    Double_t M[10];
    Int_t pdgID[10];
    Int_t status[10];
    input_tree->SetBranchAddress("mc_px", px);
    input_tree->SetBranchAddress("mc_py", py);
    input_tree->SetBranchAddress("mc_pz", pz);
    input_tree->SetBranchAddress("mc_E", E);
    input_tree->SetBranchAddress("mc_m", M);
    input_tree->SetBranchAddress("mc_pdgID", pdgID);
    input_tree->SetBranchAddress("mc_status", status);

    // output tree
    string output_tree_name = input_tree_name;
    TTree * output_tree = input_tree->CloneTree(0);

    for(int i=0; i<(int)input_tree->GetEntries(); i++)
    {
	input_tree->GetEntry(i);

	int index_e=0, index_m=0, index_j=0;
	for (int j=0; j<10; j++)
	{
	    if(status[j] != 1) continue;
	    else if (abs(pdgID[j]) == 1 || abs(pdgID[j]) == 2)  // jet
	    {
		px_j[index_j] = px[j];
		py_j[index_j] = py[j];
		pz_j[index_j] = pz[j];
		E_j[index_j] = E[j];
		M_j[index_j] = M[j];
		index_j++;
	    }
	    else if ( abs(pdgID[j]) == 11) // electron
	    {
		px_e[index_e] = px[j];
		py_e[index_e] = py[j];
		pz_e[index_e] = pz[j];
		E_e[index_e] = E[j];
		M_e[index_e] = M[j];
		index_e++;
	    }
	    else if ( abs(pdgID[j]) == 13) // muon
	    {
		px_m[index_m] = px[j];
		py_m[index_m] = py[j];
		pz_m[index_m] = pz[j];
		E_m[index_m] = E[j];
		M_m[index_m] = M[j];
		index_m++;
	    }
	}

	for (int i=0; i<2; i++)
	{
	    convert(px_e[i], py_e[i], pz_e[i], p_e[i], pt_e[i], eta_e[i], phi_e[i]);
	}
	for (int i=0; i<2; i++)
	{
	    convert(px_m[i], py_m[i], pz_m[i], p_m[i], pt_m[i], eta_m[i], phi_m[i]);
	}
	for (int i=0; i<2; i++)
	{
	    convert(px_j[i], py_j[i], pz_j[i], p_j[i], pt_j[i], eta_j[i], phi_j[i]);
	}
	// sort the four leptons pt and two jet pt
	if ( pt_e[0] < pt_e[1])
	{
	    double temp_px, temp_py, temp_pz, temp_p, temp_pt, temp_eta, temp_phi, temp_E, temp_M;
	    temp_px = px_e[0];
	    px_e[0] = px_e[1];
	    px_e[1] = temp_px;

	    temp_py = py_e[0];
	    py_e[0] = py_e[1];
	    py_e[1] = temp_py;

	    temp_pz = pz_e[0];
	    pz_e[0] = pz_e[1];
	    pz_e[1] = temp_pz;

	    temp_p = p_e[0];
	    p_e[0] = p_e[1];
	    p_e[1] = temp_p;

	    temp_pt = pt_e[0];
	    pt_e[0] = pt_e[1];
	    pt_e[1] = temp_pt;

	    temp_eta = eta_e[0];
	    eta_e[0] = eta_e[1];
	    eta_e[1] = temp_eta;

	    temp_phi = phi_e[0];
	    phi_e[0] = phi_e[1];
	    phi_e[1] = temp_phi;

	    temp_E = E_e[0];
	    E_e[0] = E_e[1];
	    E_e[1] = temp_E;
	    
	    temp_M = M_e[0];
	    M_e[0] = M_e[1];
	    M_e[1] = temp_M;
	}
	
	if ( pt_m[0] < pt_m[1])
	{
	    double temp_px, temp_py, temp_pz, temp_p, temp_pt, temp_eta, temp_phi, temp_E, temp_M;
	    temp_px = px_m[0];
	    px_m[0] = px_m[1];
	    px_m[1] = temp_px;

	    temp_py = py_m[0];
	    py_m[0] = py_m[1];
	    py_m[1] = temp_py;

	    temp_pz = pz_m[0];
	    pz_m[0] = pz_m[1];
	    pz_m[1] = temp_pz;

	    temp_p = p_m[0];
	    p_m[0] = p_m[1];
	    p_m[1] = temp_p;

	    temp_pt = pt_m[0];
	    pt_m[0] = pt_m[1];
	    pt_m[1] = temp_pt;

	    temp_eta = eta_m[0];
	    eta_m[0] = eta_m[1];
	    eta_m[1] = temp_eta;

	    temp_phi = phi_m[0];
	    phi_m[0] = phi_m[1];
	    phi_m[1] = temp_phi;

	    temp_E = E_m[0];
	    E_m[0] = E_m[1];
	    E_m[1] = temp_E;
	    
	    temp_M = M_m[0];
	    M_m[0] = M_m[1];
	    M_m[1] = temp_M;
	}

	if(pt_j[0] < pt_j[1])
	{

	    double temp_px, temp_py, temp_pz, temp_p, temp_pt, temp_eta, temp_phi, temp_E, temp_M;
	    temp_px = px_j[0];
	    px_j[0] = px_j[1];
	    px_j[1] = temp_px;

	    temp_py = py_j[0];
	    py_j[0] = py_j[1];
	    py_j[1] = temp_py;

	    temp_pz = pz_j[0];
	    pz_j[0] = pz_j[1];
	    pz_j[1] = temp_pz;

	    temp_p = p_j[0];
	    p_j[0] = p_j[1];
	    p_j[1] = temp_p;

	    temp_pt = pt_j[0];
	    pt_j[0] = pt_j[1];
	    pt_j[1] = temp_pt;

	    temp_eta = eta_j[0];
	    eta_j[0] = eta_j[1];
	    eta_j[1] = temp_eta;

	    temp_phi = phi_j[0];
	    phi_j[0] = phi_j[1];
	    phi_j[1] = temp_phi;

	    temp_E = E_j[0];
	    E_j[0] = E_j[1];
	    E_j[1] = temp_E;

	    temp_M = M_j[0];
	    M_j[0] = M_j[1];
	    M_j[1] = temp_M;
	}

	double px_ll[2], py_ll[2], pz_ll[2];
	double rap_j[2];
	double rap_llll;
	px_ll[0] = px_e[0] + px_e[1];
	py_ll[0] = py_e[0] + py_e[1];
	pz_ll[0] = pz_e[0] + pz_e[1];
	E_ll[0] = E_e[0] + E_e[1];
	convert(px_ll[0], py_ll[0], pz_ll[0], p_ll[0], pt_ll[0], eta_ll[0], phi_ll[0]);
	M_ll[0] = sqrt(pow(E_ll[0],2) - pow(p_ll[0], 2));
	
	px_ll[1] = px_m[0] + px_m[1];
	py_ll[1] = py_m[0] + py_m[1];
	pz_ll[1] = pz_m[0] + pz_m[1];
	E_ll[1] = E_m[0] + E_m[1];
	convert(px_ll[1], py_ll[1], pz_ll[1], p_ll[1], pt_ll[1], eta_ll[1], phi_ll[1]);
	M_ll[1] = sqrt(pow(E_ll[1],2) - pow(p_ll[1], 2));

	px_llll = px_ll[0] + px_ll[1];
	py_llll = py_ll[0] + py_ll[1];
	pz_llll = pz_ll[0] + pz_ll[1];
	E_llll = E_ll[0] + E_ll[1];
	convert(px_llll, py_llll, pz_llll, p_llll, pt_llll, eta_llll, phi_llll);
	M_llll = sqrt(pow(E_llll,2) - pow(p_llll, 2));
	rap_llll = 0.5*log((E_llll + pz_llll) / (E_llll - pz_llll));

	for (int i=0; i<2; i++)
	{
	    rap_j[i] = 0.5*log((E_j[i] + pz_j[i]) / (E_j[i] - pz_j[i]));
	}
	px_jj = px_j[0] + px_j[1];
	py_jj = py_j[0] + py_j[1];
	pz_jj = pz_j[0] + pz_j[1];
	E_jj = E_j[0] + E_j[1];
	convert(px_jj, py_jj, pz_jj, p_jj, pt_jj, eta_jj, phi_jj);
	M_jj = sqrt(pow(E_jj, 2) - pow(p_jj, 2));

	deltaEta_jj = abs(eta_j[0] - eta_j[1]);
	centrality = (rap_llll - (rap_j[0] + rap_j[1])/2 ) / abs(rap_j[0] - rap_j[1]);

	// lepton Eta < 2.5
	if( abs(eta_e[0])<2.5 && abs(eta_e[1])<2.5 && abs(eta_m[0])<2.5 && abs(eta_m[1]) < 2.5 )
	{
	    FLAG_cut["lepEta"] = true;
	}

	// lepton Pt > 7Gev
	if( pt_e[0]>7 && pt_e[1]>7 && pt_m[0]>7 && pt_m[1]>7 )
	{
	    FLAG_cut["lepPt"] = true;
	}

	// number of Jet (Jet Pt > 25Gev)
	if( pt_j[0]>25 && pt_j[1]>25 )
	{
	    FLAG_cut["numJet"] = true;
	}

	// Mjj(>500GeV)
	if( M_jj>500 )
	{
	    FLAG_cut["mjj"] = true;
	}

	// DeltaJetEta( > 3)
	if( deltaEta_jj>3 )
	{
	    FLAG_cut["deltaJetEta"] = true;
	}

	// Centrality(< 3) 
	if( abs(centrality)<3 )
	{
	    FLAG_cut["centrality"] = true;
	}

	// choose wanted cuts 
	if (FLAG_cut["lepEta"] && FLAG_cut["lepPt"] && FLAG_cut["numJet"] && FLAG_cut["mjj"] && FLAG_cut["deltaJetEta"] && FLAG_cut["centrality"] )
	{
	    output_tree->Fill();
	}

    }
    output_tree->AutoSave();

    time(&finish);
    double timedif = difftime(finish, start);

    cout << "Finalize: Time Cost: " << timedif << "seconds." << endl;
    cout << "Input events: " << input_tree->GetEntries() << endl;
    cout << "Output events: " << output_tree->GetEntries() << endl;
}





