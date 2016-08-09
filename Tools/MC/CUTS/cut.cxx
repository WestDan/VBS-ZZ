#include <iostream>
#include <string>
#include <time.h>
using namespace std;

typedef map<string, bool> MapType1_Bool;
typedef map<string, map<string, double> > MapType2_Double;
typedef map<string, map<string, TH1D*> > MapType2_TH1D;

time_t start, finish;

MapType1_Bool FLAG_cut;
MapType2_Double HistVar;
MapType2_TH1D histo;

vector<string> CUT;

int convert_1(Double_t px, Double_t py, Double_t pz, Double_t & p, Double_t & pt, Double_t & eta, Double_t & phi)
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
    cout << px << "\t" << py << "\t" << pz << "\t" << p << "\t" << pt << endl;
    eta = -log(sin_theta/( 1+cos_theta) );
    cout << eta << endl;
    cout << endl;
    // the return value range of atan() is [-π, π]
    phi = atan(py/px);
    if(abs(pt*cos(phi) - px) < 1 ) ;
    else 
      if( phi > 0) phi -= TMath::Pi();
      else phi += TMath::Pi();
    return 0;
}

void InitStrVec(vector<string>& out, string in, string de)
{
    int pos=0, pos_pre=0;
    while(true)
    {
	pos = in.find(de, pos_pre);
	if(pos == -1) 
	{
	    out.push_back(in.substr(pos_pre, in.size() - pos_pre));
	    break;
	}
	else out.push_back(in.substr(pos_pre, pos-pos_pre));
	pos_pre = pos + 1;
    }
}

void InitHistVar(string varlist, int nbin, double xmin, double xmax, string cutstep)
{
    vector<string> variables;
    InitStrVec(variables, varlist, ",");

    vector<string> cuts;
    if(cutstep == "") cuts.clear();
    else if(cutstep == "All") cuts = CUT; 
    else if(cutstep.compare(0, 1, "-") == 0){
	cutstep = cutstep.erase(0,1);
	for(int i=0; i<CUT.size(); i++){
	    if(cutstep != CUT[i]) cuts.push_back(CUT[i]);
	    else if(cutstep == CUT[i]) {
		cuts.push_back(CUT[i]);
		break;
	    }
	}
    } // subtract some cut form CUT
    else { InitStrVec(cuts, cutstep, ","); }

    for(int i=0; i<cuts.size(); i++)
    {
	string cut = cuts[i];
	for(int j=0; j<variables.size(); j++) 
	{
	    string var = variables[j];

	    HistVar[var]["Value"] = -9999.0;
	    HistVar[var][cut] = 1;

	    string histo_name = cut + "_" + var;
	    if(nbin == 0 ) continue;
	    if(xmin > xmax)
		cout << "Error, in InitStrVec, for histogram: " << var << "xmin should less than xmax." << endl;
	    TH1D * histo_pointer = new TH1D(histo_name.c_str(), histo_name.c_str(), nbin, xmin, xmax);
	    histo[cut][var] = histo_pointer;
	}
    }
}

void FillHistograms(float weight=1)
{
    MapType2_Double::iterator it;
    for(it=HistVar.begin(); it!=HistVar.end(); ++it)
    {
	string var = (*it).first;
	for(int i=0; i<CUT.size(); i++)
	{
	    string cut = CUT[i];
	    if(int((*it).second[cut]) == 1 && FLAG_cut[cut])
		histo[cut][var]->Fill(HistVar[var]["Value"], weight);
	}
    }

}
void cut() 
{
    InitStrVec(CUT, "all,lepEta,lepPt,numJet,mjj,deltaJetEta,centrality", ",");
    // initialize FLAG_cut
    for(int i=0; i<CUT.size(); i++)
    {
	FLAG_cut[CUT[i]] = false;
    }

    InitHistVar("pt_e1,pt_e2,pt_m1,pt_m2", 50, 0, 500, "All");
//    InitHistVar("pt_max", 50, 0, 500, "All");
    InitHistVar("pt_ll1,pt_ll2", 50, 0, 1000, "All");
    InitHistVar("E_e1,E_e2,E_m1,E_m2", 50, 0, 1000, "All");
    InitHistVar("M_ll1,M_ll2", 50, 66, 116, "All");
    InitHistVar("M_llll", 90, 100, 1000, "All");
    InitHistVar("eta_e1,eta_e2,eta_m1,eta_m2", 35, -3.5, 3.5, "All");
    InitHistVar("eta_ll1,eta_ll2", 35, -3.5, 3.5, "All");
    InitHistVar("phi_e1,phi_e2,phi_m1,phi_m2", 35, -3.5, 3.5, "All");
    InitHistVar("phi_ll1,phi_ll2", 35, -3.5, 3.5, "All");
    // Jets
    InitHistVar("pt_j1,pt_j2", 50, 0, 2000, "All");
    InitHistVar("pt_jj", 100, 0, 5000, "All");
    InitHistVar("E_j1,E_j2", 60, 0, 6000, "All");
    InitHistVar("M_jj", 60, 0, 6000, "All");
    InitHistVar("eta_j1,eta_j2", 35, -3.5, 3.5, "All");
    InitHistVar("phi_j1,phi_j2", 35, -3.5, 3.5, "All");
    InitHistVar("deltaEta_jj", 35, 0, 7, "All");
    InitHistVar("centrality", 40, -10, 10, "All");

    string file_name = "output_tree_1092.root";
    TFile * file = new TFile(file_name.c_str(), "read");
    if( ! file ) 
    {
	cout << "Error, file " << file_name << "doesn't exit" << endl;
	return ;
    }

    TTree * tree = (TTree *) file->Get("vbfnlo");
    if ( ! tree ) 
    {
	cout << "tree vbfnlo doesn't exist in file " << file_name << endl;
	return ;
    }

    TFile * out_file = new TFile("output.root", "recreate");
    if(! out_file)
    {
	cout << "Error! unable to create file: output.root" << endl;
	return;
    }

    TCanvas * cc = new TCanvas("cc", "cc", 1400, 900);

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
    double m_origin=0, m_lepEta=0, m_lepPt=0, m_numJet=0, m_mjj=0, m_deltaJetEta=0, m_centrality=0;

    Double_t px[10];
    Double_t py[10];
    Double_t pz[10];
    Double_t E[10];
    Double_t M[10];
    Int_t pdgID[10];
    Int_t status[10];
    tree->SetBranchAddress("mc_px", px);
    tree->SetBranchAddress("mc_py", py);
    tree->SetBranchAddress("mc_pz", pz);
    tree->SetBranchAddress("mc_E", E);
    tree->SetBranchAddress("mc_m", M);
    tree->SetBranchAddress("mc_pdgID", pdgID);
    tree->SetBranchAddress("mc_status", status);

    for(int i=0; i<(int)tree->GetEntries(); i++)
    {
	tree->GetEntry(i);

	int index_e=0, index_m=0, index_j=0;
	for (int j=0; j<10; j++)
	{
	    if(status[j] != 1) continue;
	    else if (abs(pdgID[j]) == 1 || abs(pdgID[j]) == 2)
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
	    convert_1(px_e[i], py_e[i], pz_e[i], p_e[i], pt_e[i], eta_e[i], phi_e[i]);
	}
	for (int i=0; i<2; i++)
	{
	    convert_1(px_m[i], py_m[i], pz_m[i], p_m[i], pt_m[i], eta_m[i], phi_m[i]);
	}
	for (int i=0; i<2; i++)
	{
	    convert_1(px_j[i], py_j[i], pz_j[i], p_j[i], pt_j[i], eta_j[i], phi_j[i]);
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
	convert_1(px_ll[0], py_ll[0], pz_ll[0], p_ll[0], pt_ll[0], eta_ll[0], phi_ll[0]);
	M_ll[0] = sqrt(pow(E_ll[0],2) - pow(p_ll[0], 2));
	
	px_ll[1] = px_m[0] + px_m[1];
	py_ll[1] = py_m[0] + py_m[1];
	pz_ll[1] = pz_m[0] + pz_m[1];
	E_ll[1] = E_m[0] + E_m[1];
	convert_1(px_ll[1], py_ll[1], pz_ll[1], p_ll[1], pt_ll[1], eta_ll[1], phi_ll[1]);
	M_ll[1] = sqrt(pow(E_ll[1],2) - pow(p_ll[1], 2));

	px_llll = px_ll[0] + px_ll[1];
	py_llll = py_ll[0] + py_ll[1];
	pz_llll = pz_ll[0] + pz_ll[1];
	E_llll = E_ll[0] + E_ll[1];
	convert_1(px_llll, py_llll, pz_llll, p_llll, pt_llll, eta_llll, phi_llll);
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
	convert_1(px_jj, py_jj, pz_jj, p_jj, pt_jj, eta_jj, phi_jj);
	M_jj = sqrt(pow(E_jj, 2) - pow(p_jj, 2));

	deltaEta_jj = abs(eta_j[0] - eta_j[1]);
	centrality = (rap_llll - (rap_j[0] + rap_j[1])/2 ) / abs(rap_j[0] - rap_j[1]);

	HistVar["pt_e1"]["Value"] = pt_e[0];
	HistVar["pt_e2"]["Value"] = pt_e[1];
	HistVar["pt_m1"]["Value"] = pt_m[0];
	HistVar["pt_m1"]["Value"] = pt_m[1];
	HistVar["eta_e1"]["Value"] = eta_e[0];
	HistVar["eta_e2"]["Value"] = eta_e[1];
	HistVar["eta_m1"]["Value"] = eta_m[0];
	HistVar["eta_m1"]["Value"] = eta_m[1];
	HistVar["phi_e1"]["Value"] = phi_e[0];
	HistVar["phi_e2"]["Value"] = phi_e[1];
	HistVar["phi_m1"]["Value"] = phi_m[0];
	HistVar["phi_m1"]["Value"] = phi_m[1];
	HistVar["E_e1"]["Value"] = E_e[0];
	HistVar["E_e2"]["Value"] = E_e[1];
	HistVar["E_m1"]["Value"] = E_m[0];
	HistVar["E_m1"]["Value"] = E_m[1];
	HistVar["M_ll1"]["Value"] = M_ll[0];
	HistVar["M_ll2"]["Value"] = M_ll[1];
	HistVar["pt_ll1"]["Value"] = pt_ll[0];
	HistVar["pt_ll2"]["Value"] = pt_ll[1];
	HistVar["eta_ll1"]["Value"] = eta_ll[0];
	HistVar["eta_ll2"]["Value"] = eta_ll[1];
	HistVar["phi_ll1"]["Value"] = phi_ll[0];
	HistVar["phi_ll2"]["Value"] = phi_ll[1];
	HistVar["M_llll"]["Value"] = M_llll;

	HistVar["pt_j1"]["Value"] = pt_j[0];
	HistVar["pt_j2"]["Value"] = pt_j[1];
	HistVar["eta_j1"]["Value"] = eta_j[0];
	HistVar["eta_j2"]["Value"] = eta_j[1];
	HistVar["phi_j1"]["Value"] = phi_j[0];
	HistVar["phi_j2"]["Value"] = phi_j[1];
	HistVar["E_j1"]["Value"] = E_j[0];
	HistVar["E_j2"]["Value"] = E_j[1];
	HistVar["pt_jj"]["Value"] = pt_jj;
	HistVar["M_jj"]["Value"] = M_jj;
	HistVar["deltaEta_jj"]["Value"] = deltaEta_jj;
	HistVar["centrality"]["Value"] = centrality;

	FLAG_cut["all"] = true;
	m_origin++;

	// lepton Eta < 2.5
	if(abs(eta_e[0])<2.5 && abs(eta_e[1])<2.5 && abs(eta_m[0])<2.5 && abs(eta_m[1]) < 2.5)
	{
	    FLAG_cut["lepEta"] = true;
	    m_lepEta++;
	}

	// lepton Pt > 7Gev
	if(FLAG_cut["lepEta"] && pt_e[0]>7 && pt_e[1]>7 && pt_m[0]>7 && pt_m[1]>7)
	{
	    FLAG_cut["lepPt"] = true;
	    m_lepPt++;
	}

	// number of Jet (Jet Pt > 25Gev)
	if(FLAG_cut["lepPt"] && pt_j[0]>25 && pt_j[1]>25)
	{
	    FLAG_cut["numJet"] = true;
	    m_numJet++;
	}

	// Mjj(>500GeV)
	if(FLAG_cut["numJet"] && M_jj>500)
	{
	    FLAG_cut["mjj"] = true;
	    m_mjj++;
	}

	// DeltaJetEta( > 3)
	if(FLAG_cut["mjj"] && deltaEta_jj>3)
	{
	    FLAG_cut["deltaJetEta"] = true;
	    m_deltaJetEta++;
	}

	// Centrality(< 3) 
	if(FLAG_cut["deltaJetEta"] && abs(centrality)<3)
	{
	    FLAG_cut["centrality"] = true;
	    m_centrality++;
	}

	FillHistograms();
    }

    time(&finish);
    double timedif = difftime(finish, start);

    m_origin = histo["all"]["pt_e1"]->Integral(0,-1);
    m_lepEta = histo["lepEta"]["pt_e1"]->Integral(0,-1);
    m_lepPt = histo["lepPt"]["pt_e1"]->Integral(0,-1);
//    m_leadingLepPt = histo["leadingLepPt"]["p3"]->Integral(0,-1);
    m_numJet = histo["numJet"]["pt_e1"]->Integral(0,-1);
    m_mjj = histo["mjj"]["pt_e1"]->Integral(0,-1);
    m_deltaJetEta = histo["deltaJetEta"]["pt_e1"]->Integral(0,-1);
    m_centrality = histo["centrality"]["pt_e1"]->Integral(0,-1);

    // Draw Histograms
    MapType2_Double::iterator it;
    for(it=HistVar.begin(); it!=HistVar.end(); ++it)
    {
	string var = (*it).first;
	for(int i=0; i<CUT.size(); i++)
	{
	    string cut = CUT[i];
	    if(HistVar[var][cut])
	    {
		string Histname = cut + "_" + var;
		histo[cut][var]->Write(Histname.c_str());
		histo[cut][var]->Draw();
		string output = "Plot/" + cut + "_" + var + ".png";
		cc->SaveAs(output.c_str());
	    }
	}
    }

    cout << "Finalize: Time Cost: " << timedif << "seconds." << endl;
    cout << "Origin: " << m_origin << endl;
    cout << "lepEta: " << m_lepEta << endl;
    cout << "lepPt: " << m_lepPt << endl;
//    cout << "leadingLepPt: " << m_leadingLepPt << endl;
    cout << "numJet: " << m_numJet << endl;
    cout << "mjj: " << m_mjj << endl;
    cout << "deltaJetEta: " << m_deltaJetEta << endl;
    cout << "centrality: " << m_centrality << endl;

}





