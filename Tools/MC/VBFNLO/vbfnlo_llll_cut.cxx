#include <iostream>
#include <string>
#include <time.h>
using namespace std;

// Is there any way to directly operate trees, rather than copy
// part of one tree to another.

//	 what this program do 
// tree =====================> cutted tree
// cut events, using macro for all cut conditions.
typedef map<string, bool> MapType1_Bool;
typedef map<string, map<string, double> > MapType2_Double;

static double Z_MASS = 91.1876;  // unit in GeV
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

void swap(double & x, double & y)
{
    double temp = x;
    x = y;
    y = temp;
}
void vbfnlo_llll_cut() 
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

    Double_t px_l[4], py_l[4], pz_l[4];
    Double_t p_l[4], pt_l[4];
    Double_t eta_l[4], phi_l[4];
    Double_t E_l[4], M_l[4];
    Int_t    parent1_index_l[4], parent2_index_l[4];
    Double_t p_Z[2], pt_Z[2], eta_Z[2], phi_Z[2];
    Double_t E_Z[2], M_Z[2];
    Double_t p_ZZ, pt_ZZ;
    Double_t E_ZZ, M_ZZ;
    Double_t eta_ZZ, phi_ZZ;
    Double_t px_j[2], py_j[2], pz_j[2];   // Jets
    Double_t p_j[2], pt_j[2];
    Double_t E_j[2], M_j[2];
    Double_t eta_j[2], phi_j[2];
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
    Int_t parent1_index[10];
    Int_t parent2_index[10];
    input_tree->SetBranchAddress("mc_px", px);
    input_tree->SetBranchAddress("mc_py", py);
    input_tree->SetBranchAddress("mc_pz", pz);
    input_tree->SetBranchAddress("mc_E", E);
    input_tree->SetBranchAddress("mc_m", M);
    input_tree->SetBranchAddress("mc_pdgID", pdgID);
    input_tree->SetBranchAddress("mc_status", status);
    input_tree->SetBranchAddress("mc_parent1_index", parent1_index);
    input_tree->SetBranchAddress("mc_parent2_index", parent2_index);

    // output tree
    string output_tree_name = input_tree_name;
    // output trees with branch in standard name.
//    TTree * output_tree = input_tree->CloneTree(0);
    TTree * output_tree = new TTree("vbfnlo_ZZ", " vbfnlo 4leptons events");
    // set output tree branches according to standard interface.
    output_tree->Branch("pt_l1", pt_l,   "pt_l1/D");
    output_tree->Branch("pt_l2", pt_l+1, "pt_l2/D");
    output_tree->Branch("pt_l3", pt_l+2, "pt_l3/D");
    output_tree->Branch("pt_l4", pt_l+3, "pt_l4/D");
    output_tree->Branch("eta_l1", eta_l,   "eta_l1/D");
    output_tree->Branch("eta_l2", eta_l+1, "eta_l2/D");
    output_tree->Branch("eta_l3", eta_l+2, "eta_l3/D");
    output_tree->Branch("eta_l4", eta_l+3, "eta_l4/D");
    output_tree->Branch("phi_l1", phi_l,   "phi_l1/D");
    output_tree->Branch("phi_l2", phi_l+1, "phi_l2/D");
    output_tree->Branch("phi_l3", phi_l+2, "phi_l3/D");
    output_tree->Branch("phi_l4", phi_l+3, "phi_l4/D");
    output_tree->Branch("E_l1", E_l,   "E_l1/D");
    output_tree->Branch("E_l2", E_l+1, "E_l2/D");
    output_tree->Branch("E_l3", E_l+2, "E_l3/D");
    output_tree->Branch("E_l4", E_l+3, "E_l4/D");

    output_tree->Branch("pt_Z1", pt_Z,   "pt_Z1/D");
    output_tree->Branch("pt_Z2", pt_Z+1, "pt_Z2/D");
    output_tree->Branch("eta_Z1", eta_Z,   "eta_Z1/D");
    output_tree->Branch("eta_Z2", eta_Z+1, "eta_Z2/D");
    output_tree->Branch("phi_Z1", phi_Z,   "phi_Z1/D");
    output_tree->Branch("phi_Z2", phi_Z+1, "phi_Z2/D");
    output_tree->Branch("E_Z1", E_Z,   "E_Z1/D");
    output_tree->Branch("E_Z2", E_Z+1, "E_Z2/D");
    output_tree->Branch("M_Z1", M_Z,   "M_Z1/D");
    output_tree->Branch("M_Z2", M_Z+1, "M_Z2/D");

    output_tree->Branch("pt_ZZ", &pt_ZZ, "pt_ZZ/D");
    output_tree->Branch("eta_ZZ", &eta_ZZ, "eta_ZZ/D");
    output_tree->Branch("phi_ZZ", &phi_ZZ, "phi_ZZ/D");
    output_tree->Branch("E_ZZ", &E_ZZ, "E_ZZ/D");
    output_tree->Branch("M_ZZ", &M_ZZ, "M_ZZ/D");

    output_tree->Branch("pt_j1", pt_j,   "pt_j1/D");
    output_tree->Branch("pt_j2", pt_j+1, "pt_j2/D");
    output_tree->Branch("eta_j1", eta_j,   "eta_j1/D");
    output_tree->Branch("eta_j2", eta_j+1, "eta_j2/D");
    output_tree->Branch("phi_j1", phi_j,   "phi_j1/D");
    output_tree->Branch("phi_j2", phi_j+1, "phi_j2/D");
    output_tree->Branch("E_j1", E_j,   "E_j1/D");
    output_tree->Branch("E_j2", E_j+1, "E_j2/D");
    output_tree->Branch("M_j1", M_j,   "M_j1/D");
    output_tree->Branch("M_j2", M_j+1, "M_j2/D");

    output_tree->Branch("pt_jj", &pt_jj, "pt_jj/D");
    output_tree->Branch("eta_jj", &eta_jj, "eta_jj/D");
    output_tree->Branch("phi_jj", &phi_jj, "phi_jj/D");
    output_tree->Branch("E_jj", &E_jj, "E_jj/D");
    output_tree->Branch("M_jj", &M_jj, "M_jj/D");

    output_tree->Branch("deltaEta_jj", &deltaEta_jj, "deltaEta_jj/D");
    output_tree->Branch("centrality", &centrality, "centrality/D");

    for(int i=0; i<(int)input_tree->GetEntries(); i++)
    {
	// Initialize variables
	for (int i=0; i<4; i++)
	{
	    px_l[i] = -9999E10;
	    py_l[i] = -9999E10;
	    pz_l[i] = -9999E10;
	    p_l[i] = -9999E10;
	    pt_l[i] = -9999E10;
	    eta_l[i] = -9999E10;
	    phi_l[i] = -9999E10;
	    E_l[i] = -9999E10;
	    M_l[i] = -9999E10;
	    parent1_index_l[i] = -1;
	    parent2_index_l[i] = -1;
	}

	for (int i=0; i<2; i++)
	{
	    p_Z[i] = -9999E10;
	    pt_Z[i] = -9999E10;
	    eta_Z[i] = -9999E10;
	    phi_Z[i] = -9999E10;
	    E_Z[i] = -9999E10;
	    M_Z[i] = -9999E10;

	    px_j[i] = -9999E10;
	    py_j[i] = -9999E10;
	    pz_j[i] = -9999E10;
	    p_j[i] = -9999E10;
	    pt_j[i] = -9999E10;
	    eta_j[i] = -9999E10;
	    phi_j[i] = -9999E10;
	    E_j[i] = -9999E10;
	    M_j[i] = -9999E10;
	}
	
	p_ZZ = -9999E10;
	pt_ZZ = -9999E10;
	eta_ZZ = -9999E10;
	phi_ZZ = -9999E10;
	E_ZZ = -9999E10;
	M_ZZ = -9999E10;

	p_jj = -9999E10;
	pt_jj = -9999E10;
	eta_jj = -9999E10;
	phi_jj = -9999E10;
	E_jj = -9999E10;
	M_jj = -9999E10;

	deltaEta_jj = -9999E10;
	centrality = -9999E10;

	double rap_j[2]={-9999E10, -9999E10}, rap_ZZ=-9999E10;

	FLAG_cut["lepEta"] = false;
	FLAG_cut["lepPt"] = false;
	FLAG_cut["mZ1"] = false;
	FLAG_cut["mZ2"] = false;
	FLAG_cut["numJet"] = false;
	FLAG_cut["mjj"] = false;
	FLAG_cut["deltaJetEta"] = false;
	FLAG_cut["centrality"] = false;

	input_tree->GetEntry(i);

	int index_l=0, index_e=0, index_m=0, index_j=0;
	for (int j=0; j<10; j++)
	{
	    if(status[j] != 1) continue;
	    else if (abs(pdgID[j]) == 1 || abs(pdgID[j]) == 2 || abs(pdgID[j]) == 3 || abs(pdgID[j]) == 4)  // jet
	    {
		px_j[index_j] = px[j];
		py_j[index_j] = py[j];
		pz_j[index_j] = pz[j];
		E_j[index_j] = E[j];
		M_j[index_j] = M[j];
		index_j++;
	    }
	    else if ( abs(pdgID[j]) == 11 ) // electron
	    {
		px_l[index_l] = px[j];
		py_l[index_l] = py[j];
		pz_l[index_l] = pz[j];
		E_l[index_l] = E[j];
		M_l[index_l] = M[j];
                parent1_index_l[index_l] = parent1_index[j];
                parent2_index_l[index_l] = parent2_index[j];
		index_e++;
		index_l++;
	    }
	    else if ( abs(pdgID[j]) == 13) // muon
	    {
		px_l[index_l] = px[j];
		py_l[index_l] = py[j];
		pz_l[index_l] = pz[j];
		E_l[index_l] = E[j];
		M_l[index_l] = M[j];
                parent1_index_l[index_l] = parent1_index[j];
                parent2_index_l[index_l] = parent2_index[j];
		index_m++;
		index_l++;
	    }
	}


	if( (index_e == 2 && index_m ==2) || (index_e == 4 && index_m == 0) || (index_e == 0 && index_m == 4) )
	{
	    for (int i=0; i<4; i++)
	    {
		convert(px_l[i], py_l[i], pz_l[i], p_l[i], pt_l[i], eta_l[i], phi_l[i]);
	    }

	// pairs leptons, look firstly at their pID, then their parent 
	    double px_Z[2], py_Z[2], pz_Z[2];
	    int index_Z = 0;
	    for ( int i=0; i<4; i++)
	    {
		for (int j=i+1; j<4; j++)
		{
		    if(parent1_index_l[j] == parent1_index_l[i]) // do we need to compare their pID ???
		    {
			px_Z[index_Z] = px_l[i] + px_l[j];
			py_Z[index_Z] = py_l[i] + py_l[j];
			pz_Z[index_Z] = pz_l[i] + pz_l[j];
			E_Z[index_Z] = E_l[i] + E_l[j];
			convert(px_Z[index_Z], py_Z[index_Z], pz_Z[index_Z], p_Z[index_Z], pt_Z[index_Z], eta_Z[index_Z], phi_Z[index_Z]);
			M_Z[index_Z] = sqrt(pow(E_Z[index_Z],2) - pow(p_Z[index_Z], 2));
			index_Z++;
		    }
		}
	    }

	    if(index_Z == 2)
	    {
		double px_ZZ = px_Z[0] + px_Z[1];
		double py_ZZ = py_Z[0] + py_Z[1];
		double pz_ZZ = pz_Z[0] + pz_Z[1];
		E_ZZ = E_Z[0] + E_Z[1];
		convert(px_ZZ, py_ZZ, pz_ZZ, p_ZZ, pt_ZZ, eta_ZZ, phi_ZZ);
		M_ZZ = sqrt(pow(E_ZZ,2) - pow(p_ZZ, 2));
		rap_ZZ = 0.5*log((E_ZZ + pz_ZZ) / (E_ZZ - pz_ZZ));

		// sort the two Z particles
		if ( abs(M_Z[1] - Z_MASS) < abs(M_Z[0] - Z_MASS))
		{
		    swap(pt_Z[0], pt_Z[1]);
		    swap(eta_Z[0], eta_Z[1]);
		    swap(phi_Z[0], phi_Z[1]);
		    swap(E_Z[0], E_Z[1]);
		    swap(M_Z[0], M_Z[1]);
		}
	    }

	    // sort leptons
	    for( int i=0; i<4; i++ )
	    {
		for( int j=i+1; j<4; j++ )
		{
		    if ( pt_l[j] > pt_l[i] )
		    {
			// only care about those important parameters, px, py , pz and parent_index are discarded.
			swap(pt_l[i], pt_l[j]);
			swap(eta_l[i], eta_l[j]);
			swap(phi_l[i], phi_l[j]);
			swap(E_l[i], E_l[j]);
			swap(M_l[i], M_l[j]);
		    }
		}
	    }
	}
	else 
	{
	    cout << " event " << i << " doesn't have 4 leptons";
	    continue;
	}

    
	if( index_j == 2)
	{
	    for (int i=0; i<2; i++)
	    {
		convert(px_j[i], py_j[i], pz_j[i], p_j[i], pt_j[i], eta_j[i], phi_j[i]);
	    }

	    double px_jj = px_j[0] + px_j[1];
	    double py_jj = py_j[0] + py_j[1];
	    double pz_jj = pz_j[0] + pz_j[1];
	    E_jj = E_j[0] + E_j[1];
	    convert(px_jj, py_jj, pz_jj, p_jj, pt_jj, eta_jj, phi_jj);
	    M_jj = sqrt(pow(E_jj, 2) - pow(p_jj, 2));

	    if(pt_j[0] < pt_j[1])
	    {
		swap(pt_j[0], pt_j[1]);
		swap(eta_j[0], eta_j[1]);
		swap(phi_j[0], phi_j[1]);
		swap(E_j[0], E_j[1]);
		swap(M_j[0], M_j[1]);
	    }

	    for (int i=0; i<2; i++)
	    {
		rap_j[i] = 0.5*log((E_j[i] + pz_j[i]) / (E_j[i] - pz_j[i]));
	    }

	    deltaEta_jj = abs(eta_j[0] - eta_j[1]);
	}
	else 
	{
	    cout << "jet number less than 2";
	    continue;
	}

	centrality = (rap_ZZ - (rap_j[0] + rap_j[1])/2 ) / abs(rap_j[0] - rap_j[1]);

	// lepton Eta < 2.5
	if( abs(eta_l[0])<2.5 && abs(eta_l[1])<2.5 && abs(eta_l[2])<2.5 && abs(eta_l[3]) < 2.5 )
	{
	    FLAG_cut["lepEta"] = true;
	}

	// lepton Pt > 7Gev
	if( pt_l[0]>7 && pt_l[1]>7 && pt_l[2]>7 && pt_l[3]>7 )
	{
	    FLAG_cut["lepPt"] = true;
	}

	// Z mass
	if( abs(M_Z[0] - Z_MASS) < 25 )
	{
	    FLAG_cut["mZ1"] = true;
	}
	if( abs(M_Z[1] - Z_MASS) < 25 )
	{
	    FLAG_cut["mZ2"] = true;
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
	if (FLAG_cut["lepEta"] && FLAG_cut["lepPt"] && FLAG_cut["mZ1"]  && FLAG_cut["mZ2"] && FLAG_cut["numJet"] && FLAG_cut["mjj"] && FLAG_cut["deltaJetEta"] && FLAG_cut["centrality"] )
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





