#include <iostream>
#include <string>
#include <time.h>
using namespace std;

// this script extract histogram from trees of input files
void hist_extraction()
{
    // input file
    string input_file_name = "output.root";
    TFile * input_file = new TFile(input_file_name.c_str(), "read");
    if (! input_file)
    {
	cout << "Error, input file " << input_file_name << " doesn't exist" << endl;
	return 1;
    }

    // specify variables and cuts 
    // cut order: all -> lepPt -> mZ1 -> mZ2 -> numJet -> mjj -> deltaJetEta -> centrality
    // important variables: 
    //	    lepton: pt(pt_l1..pt_l4), eta, phi
    //      Z: pt, eta, phi, m
    //      ZZ: m
    //      Jet: pt, eta, phi, m
    //      mjj, deltaJetEta, centrality 
    
}
