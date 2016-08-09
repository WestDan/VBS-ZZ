#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

struct entry {
  double valx;
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
      while (! infile.eof() ) {    // as long as we haven't reached the end of the file, do ...
        entry vals;
        infile >> vals.valx >> vals.valy >> vals.vale;   // read in the 3 values (x, y, error)
        if (infile.eof()) break;                        // break when at end of file
        if (numreadlo==1) {                             // 1st order input
          histentrylo.push_back(vals);                 // extend datarray to row of 3 numbers
        } else {
          if (histentrylo[k].valx != vals.valx) {          // nur machen, wenn x-werte der inputfiles übereinstimmen
            cerr << "Entry mismatch: " << filenamelo << " " << histentrylo[k].valx << " <-> " << vals.valx << endl;
            return 1;
          }
          histentrylo[k].valy += vals.valy;              //outputeintrag=summe der inputeinträge
          histentrylo[k].vale += vals.vale*vals.vale;	  //für fehler quadratisch 
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
  ofstream outfile(fileout.c_str());                        //outputfile öffnen 
  if (! outfile.is_open()) {
    cerr << "Cannot open file " << fileout << endl;
    return 1;
  }

  for (unsigned int i=0; i<histentrylo.size(); i++) {
    outfile << histentrylo[i].valx << " " << histentrylo[i].valy/static_cast<double>(numreadlo) << " " << sqrt(histentrylo[i].vale)/static_cast<double>(numreadlo) << endl;
  }                 // schreibe raus: xwert, addierte ywerte / # der inputs  , wurzel der quadrat. addierten fehler/  # der inputs (nicht: /(n-1)?)
  cout << "Written " << numreadlo << " LO entries for file " << histfile << endl;
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
        infile >> vals.valx >> vals.valy >> vals.vale;  
        if (infile.eof()) break;
        if (numread==1) {
          histentrynlo.push_back(vals);
        } else {
          if (histentrynlo[k].valx != vals.valx) {
            cerr << "Entry mismatch: " << filenamenlo << " " << histentrynlo[k].valx << " <-> " << vals.valx << endl;
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
    outfile << histentrynlo[i].valx << " " << histentrynlo[i].valy/static_cast<double>(numread) << " " << sqrt(histentrynlo[i].vale)/static_cast<double>(numread) << endl;
  }
  cout << "Written " << numread << " NLO entries for file " << histfile << endl;
  outfile.close();
  
// Kfac
  fileout = outdir+"/"+histdir+"/Kfac/"+histfile;
  outfile.open(fileout.c_str());
  if (! outfile.is_open()) {
    cerr << "Cannot open file " << fileout << endl;
    return 1;
  }

  for (unsigned int i=0; i<histentrylo.size(); i++) {
    if ( histentrylo[i].valx !=  histentrynlo[i].valx ) {
      cerr << "Mismatch in LO vs. NLO bins for file " << fileout << " " << histentrylo[i].valx << " <-> " << histentrynlo[i].valx << endl;
      return 1;
    }
    if (numreadlo != numread) {
      cerr << "Error in number of entries LO vs. NLO " << numreadlo << " <-> " << numread << endl;
      return 1;
    }
    outfile << histentrynlo[i].valx << " ";
    if (histentrylo[i].valy == 0.) {
      outfile << -1. << " " << 0;
    } else {
      outfile << histentrynlo[i].valy/histentrylo[i].valy << " " << 0;
    }
    outfile << endl;
  }
  cout << "Written " << numread << " K-factor entries for file " << histfile << endl;
  outfile.close();
  
  return 0;
}
