// AnalysisVar.h
// Aim to provide variables and objects for analysis

#ifndef AnalysisVar_h
#define AnalysisVar_h

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <map>
#include <vector>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

using namespace std;

typedef struct {
    float num;
    float wnum;
    float err;
} COUNT;

typedef map<string, COUNT> MapType_Counting;
typedef map<string, map<string, COUNT> > MapType2_Counting;
typedef map<string, map<string, map<string, COUNT> > > MapType3_Counting;
typedef map<string, map<string, map<string, map<string, COUNT> > > > MapType4_Counting;

typedef map<string, map<string, int> > MapType2_Int;
typedef map<string, map<string, map<string, int> > > MapType3_Int;
typedef map<string, map<string, vector<int> > > MapType2_VInt;

typedef map<string, map<string, unsigned long long> > MapType2_ULong64;

typedef map<string, float> MapType_Float;
typedef map<string, map<string, float> > MapType2_Float;
typedef map<string, map<string, vector<float> > > MapType2_VFloat;

typedef map<string, map<string, double> > MapType2_Double;
typedef map<string, map<string, vector<double> > > MapType2_VDouble;

typedef map<string, map<string, string> > MapType2_String;
typedef map<string, vector<string> > MapType_VString;

typedef map<string, TTree*> MapType_TTree;

namespace SB{
    enum {
	_signal,
	_wz,
	_z
    };
}
class AnalysisVar 
{
public:
    vector<string> SETNAME, SYSNAME;
    MapType2_Int SETTING;

    MapType_VString STEP_obj;

    MapType3_Counting CNT_obj;

    MapType2_Int	TreeIntVar;
    MapType2_ULong64	TreeUloVar;
    MapType2_Float      TreeFltVar;
    MapType2_String	TreeStrVar;
    MapType2_VInt	TreeIntVVar;
    MapType2_VFloat	TreeFltVVar;
    MapType2_VDouble	TreeDblVVar;

    MapType_TTree Tree;
    TTree * tree, * tree_bkg;
};

#endif
