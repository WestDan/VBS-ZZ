#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

struct entry {
  double valx1;
  double valx2;
  double valy;
  double vale;
};

int main( int argc, const char* argv[] ) {
// read arguments
// arguments are hist-file, output-dir, hist-dir, input-dirs
  if (argc < 5) {
    cerr << "Usage: " << argv[0] << " <hist-file> <output-dir> <hist-dir> <input-dirs> <...> " << endl;
    return 1;
  }
  string histfile = argv[1];
  string outdir = argv[2];
  string histdir = argv[3];
  vector<string> indirs;
  for (int i=4;i<argc;i++) {
    indirs.push_back(argv[i]);
  }

  string dummy;
// LO
  int numreadlo = 0;
  vector<entry> histentrylo;
  for (unsigned int i=0;i<indirs.size();i++) {
    string filenamelo = indirs[i]+"/"+histdir+"/LO/"+histfile;
    ifstream infile(filenamelo.c_str());
    if (! infile.is_open()) {
      cerr << "Cannot open file " << filenamelo << endl;
      continue;
    }
    getline(infile,dummy);

    infile.exceptions(ifstream::badbit);
    numreadlo++;
    try{
      int k=0;
      while (! infile.eof() ) {   
        entry vals;
        infile >> vals.valx1 >> vals.valx2 >> vals.valy >> vals.vale;  
        if (infile.eof()) break;                        
        if (numreadlo==1) {                           
          histentrylo.push_back(vals);               
        } else {
          if (histentrylo[k].valx1 != vals.valx1) {       
            cerr << "Entry A mismatch: " << filenamelo << " " << histentrylo[k].valx1 << " <-> " << vals.valx1 << endl;
            return 1;
	  }
          if (histentrylo[k].valx2 != vals.valx2) {       
            cerr << "Entry B mismatch: " << filenamelo << " " << histentrylo[k].valx2 << " <-> " << vals.valx2 << " <-> " << vals.valx1 << " " << histentrylo[k].valx1 <<endl;
            return 1;
          }
          histentrylo[k].valy += vals.valy;              
          histentrylo[k].vale += vals.vale*vals.vale;	 
	}
        k++;
      }
    } catch (ifstream::failure e) {
      cerr << "Error reading file " << filenamelo << endl;
      return 1;
    }
    infile.close();
  }

  string fileout = outdir+"/"+histdir+"/LO/"+histfile;
  ofstream outfile(fileout.c_str());                        
  if (! outfile.is_open()) {
    cerr << "Cannot open file " << fileout << endl;
    return 1;
  }

  for (unsigned int i=0; i<histentrylo.size(); i++) {
    outfile << histentrylo[i].valx1 << " " << histentrylo[i].valx2 << " " << histentrylo[i].valy/static_cast<double>(numreadlo) << " " << sqrt(histentrylo[i].vale)/static_cast<double>(numreadlo) << endl;
  }                
  cout << "Written " << numreadlo << " 2d LO entries for file " << histfile << endl;
  outfile.close();
 
// NLO
  int numread = 0;
    vector<entry> histentrynlo;
  for (unsigned int i=0;i<indirs.size();i++) {
    string filenamenlo = indirs[i]+"/"+histdir+"/NLO/"+histfile;
    ifstream infile(filenamenlo.c_str());
    if (! infile.is_open()) {
      cerr << "Cannot open file " << filenamenlo << endl;
      continue;
    }
    getline(infile,dummy);

    infile.exceptions(ifstream::badbit);
    numread++;
    try{
      int k=0;
      while (! infile.eof() ) {
        entry vals;
        infile >> vals.valx1 >> vals.valx2 >> vals.valy >> vals.vale;  
        if (infile.eof()) break;
        if (numread==1) {
          histentrynlo.push_back(vals);
        } else {
          if (histentrynlo[k].valx1 != vals.valx1) {
            cerr << "Entry C mismatch: " << filenamenlo << " " << histentrynlo[k].valx1 << " <-> " << vals.valx1 << endl;
            return 1;
	  }
	  if (histentrynlo[k].valx2 != vals.valx2) {
            cerr << "Entry D mismatch: " << filenamenlo << " " << histentrynlo[k].valx2 << " <-> " << vals.valx2 << endl;
            return 1;
          }
          histentrynlo[k].valy += vals.valy;
	  histentrynlo[k].vale += vals.vale*vals.vale;
        }
        k++;
      }
    } catch (ifstream::failure e) {
      cerr << "Error reading file " << filenamenlo << endl;
      return 1;
    }
    infile.close();
  }

  fileout = outdir+"/"+histdir+"/NLO/"+histfile;
  outfile.open(fileout.c_str());
  if (! outfile.is_open()) {
    cerr << "Cannot open file " << fileout << endl;
    return 1;
  }

  for (unsigned int i=0; i<histentrynlo.size(); i++) {
    outfile << histentrynlo[i].valx1 << " " << histentrynlo[i].valx2 << " " << histentrynlo[i].valy/static_cast<double>(numread) << " " << sqrt(histentrynlo[i].vale)/static_cast<double>(numread) << endl;
  }
  cout << "Written " << numread << " 2d NLO entries for file " << histfile << endl;
  outfile.close();
  
// Kfac
  fileout = outdir+"/"+histdir+"/Kfac/"+histfile;
  outfile.open(fileout.c_str());
  if (! outfile.is_open()) {
    cerr << "Cannot open file " << fileout << endl;
    return 1;
  }

  for (unsigned int i=0; i<histentrylo.size(); i++) {
    if ( histentrylo[i].valx1 !=  histentrynlo[i].valx1 ) {
      cerr << "Mismatch in LO vs. NLO bins (x1) for file " << fileout << " " << histentrylo[i].valx1 << " <-> " << histentrynlo[i].valx1 << endl;
      return 1;
    }
    if ( histentrylo[i].valx2 !=  histentrynlo[i].valx2 ) {
      cerr << "Mismatch in LO vs. NLO bins (x2) for file " << fileout << " " << histentrylo[i].valx2 << " <-> " << histentrynlo[i].valx2 << endl;
      return 1;
    }
    if (numreadlo != numread) {
      cerr << "Error in number of entries LO vs. NLO " << numreadlo << " <-> " << numread << endl;
      return 1;
    }
    outfile << histentrynlo[i].valx1 << " " << histentrynlo[i].valx2 << " ";
    if (histentrylo[i].valy == 0.) {
      outfile << -1. << " " << 0;
    } else {
      outfile << histentrynlo[i].valy/histentrylo[i].valy << " " << 0;
    }
    outfile << endl;
  }
  cout << "Written " << numread << " 2d K-factor entries for file " << histfile << endl;
  outfile.close();
  
  return 0;
}
