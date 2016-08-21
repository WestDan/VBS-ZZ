//<AnalysisVar.h>
//<Aim to provide variables and objects for analysis>
//<Yusheng WU, April 2011, Ann Arbor>

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
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>

//<Use STD map datatypes to manage counting and histograms>
using namespace std;

//<Use this struct for base element in cut flow, {num,error}>
typedef struct{
    float num;
    float wnum;
    float err;
} COUNT; 

//<Map type STRING:COUNT...>
typedef map<string, COUNT> MapType_Counting;
typedef map<string, map<string,COUNT> > MapType2_Counting;
typedef map<string, map<string,map<string,COUNT> > > MapType3_Counting;
typedef map<string, map<string, map<string, map<string,COUNT> > > > MapType4_Counting;

//<Map type STRING:INT...>
typedef map<string, int> MapType_Int;
typedef map<string, map<string,int> > MapType2_Int;
typedef map<string, map<string, vector<int> > > MapType2_VInt;
typedef map<string, map<string, map<string,int> > > MapType3_Int;

//<Map type STRING: UNSIGNED LONG LONG INT...>
typedef map<string, map<string, unsigned long long> > MapType2_ULong64;

//<Map type STRING:STRING...>
typedef map<string, string> MapType_String;
typedef map<string, map<string,string> > MapType2_String;
typedef map<string, vector<string> > MapType_VString;

//<Map type STRING:FLOAT...>
typedef map<string, float> MapType_Float;
typedef map<string, map<string, float> > MapType2_Float;
typedef map<string, map<string, vector<float> > > MapType2_VFloat;

//<Map type STRING:DOUBLE...>
typedef map<string, double> MapType_Double;
typedef map<string, vector<double> > MapType_VDouble;
typedef map<string, vector<double>* > MapType_VDoubleStar;
typedef map<string, map<string, double> > MapType2_Double;
typedef map<string, map<string, vector<double> > > MapType2_VDouble;

typedef map<string, map<string, pair<double,double> > > MapType2_Double2D;
typedef map<string, map<string, vector<pair<double,double> > > > MapType2_V2DDouble;

//<Map type STRING:TH1F...>
typedef map<string, TH1F*> MapType_TH1F;
typedef map<string, map<string,TH1F*> > MapType2_TH1F;
typedef map<string, map<string,map<string,TH1F*> > > MapType3_TH1F;
typedef map<string, map<string,map<string,map<string,TH1F*> > > > MapType4_TH1F;
typedef map<string, vector<TH1F*> > MapType_VTH1F;
typedef map<string, map<string,vector<TH1F*> > > MapType2_VTH1F;
typedef map<string, map<string,map<string,vector<TH1F*> > > > MapType3_VTH1F;

typedef map<string, map<string,map<string,map<string,TH1D*> > > > MapType4_TH1D;

//<Map type STRING:TH2F...>
typedef map<string, TH2F*> MapType_TH2F;
typedef map<string, map<string,TH2F*> > MapType2_TH2F;
typedef map<string, map<string,map<string,TH2F*> > > MapType3_TH2F;
typedef map<string, map<string,map<string,map<string, TH2F*> > > > MapType4_TH2F;

//<Map type STRING:TTree...>
typedef map<string, TTree*> MapType_TTree;

//<Options for Variables to determine using which event weight to fill into histograms>
//<Both means both histograms and trees>
enum {InBoth,InBothNoPileup,InBothNoSF,InBothNoPileupNoSF,InTreeOnly,InHistoOnly};

namespace SB{
   enum {
      _signal,
      _wz,
      _z
   };
}

namespace quadType{
   enum {
      _2mu,
      _2e
   };
}

namespace TriLepTye{
   enum {
      _3mu,
      _2mu1e,
      _3e,
      _2e1mu
   };
}


//<AnalysisVar Class>
class AnalysisVar {

    public:        
    
        //<Settings in the analysis>
        //<SETTING => setting name: {set1:var1, set2:var2}>
        //<CUTVAR => setting name: {set1:var1, set2:var2}>
        vector<string> SETNAME, SYSNAME;
        vector<string> Var_Map; // variation values
        MapType_Double Var_List; // variation values
        MapType2_Int SETTING;
        MapType2_Double CUTVAR;

        //<Define analysis channels>
        vector<string> CHN;       
    
        //<Define analysis steps / cut steps>
        vector<string> STEP_cut; // For event selection
	vector<string> TRUTH_STEP_cut;
        vector<string> STEP_cut_tree; // For filling selection
        //MapType_VString STEP_cut; // For event selection
        MapType_VString STEP_obj; // For object selection
 
        //<Define interested variables, to be filled in trees or histograms>
        //setting name: Name: {NBins, Xmin, Xmax, CutStep, Option,VectorOrNot}
        //NBins, Xmin, Xmax -> used to create histogram for this variable
        //CutStep -> used to decide at which step the variable will be filled in histograms,
        //positive number means this step and after, negative number means exactly at this step
        //Option -> in histogram/in tree/ using which weight...
        //VectorOrNot -> this variable shall be vector or not
        MapType2_Double VarName; 
       
        //<Define Cut Flow Maps>
        //1. Event Counting Map: CNT_XX
        //"setting name" : "Channel" : "CutName" : {number of events, statistic errors}
        //"setting name" : "obj" : "CutName" : {number of events, statistic errors}
        //2. PassFlag for Selection: FLAG_XX
        //"Channel/Obj" : "CutName" : 0 or 1
        //*. Use FLAG_XX_temp first to compute decision
        MapType3_Counting CNT_obj; 
        MapType4_Counting CNT_cut; 
        //MapType2_Counting CNT_cut, CNT_obj; 
        //MapType2_Int     FLAG_cut,FLAG_cut_temp; 
        MapType3_Int     FLAG_cut_temp; 
        MapType3_Int     FLAG_cut;
	MapType3_Int	 FLAG_Truth_cut_temp;
	MapType3_Int	 FLAG_Truth_cut;  // added by weibin
        //MapType2_Counting CNT_obj;  //changed by Cong 
        MapType2_Int     FLAG_obj,FLAG_obj_temp;
                
        //<Define the event weight map>
        //i.e. event weights can be different at different cut stages
        //setting name : CutName : (Generator Weight, Pileup Weight, Reco SF, Trigger SF, Bunch/Train Weight ...)
        //Make sure you understand exactly what each element of the array means
        MapType2_Double Evt_Weight;
                
        //<Define Tree Variables and Histogram Maps>        
        //Fill the variables defined in "VarName" into both trees and histograms
        //1. Tree Variable
        //"Variable" : double
        //2. Histogram Map
        //<- TH1F>
        //"SettingName" : "Channel" : "CutName" : "Variable" : TH1F*
		//<- TH1F slices>
		//"SettingName" : "Channel" : "CutName" : "Var1_Var2Bin#" : TH1F*
        //<- TH2F>
        //"SeetingName" : "Channel" : "CutName" : "Var1_Var2" : TH2F*
        //<- TH2F Helper>
        //"Var1_Var2" : "Var1"{1,xbins,xmin,xmax},"Var2"{2,ybins,ymin,ymax}
        //MapType_Double Var;
        MapType2_Int TreeIntVar;
        MapType2_ULong64 TreeUloVar;
        MapType2_String TreeStrVar;
        MapType2_Float TreeFltVar;
        MapType2_VFloat TreeFltVVar;
        MapType2_VInt TreeIntVVar;
        MapType2_VDouble TreeDblVVar;


        MapType2_Double HistVar;
        MapType2_Double TruthHistVar;  // added by weibin 
        MapType2_Double2D Hist2DVar;
        MapType2_VDouble VVar;
        MapType2_V2DDouble V2DVar;
        //MapType_VDoubleStar VVar;
        //MapType4_TH1F histo, histo_slices;
        MapType4_TH1D histo, histo_slices;
        //MapType3_TH1F histo;

        MapType4_TH2F histo_2D, histo_v2D;
        MapType2_VDouble helper_2D, VarSlices;
        
        //<Eventually, you'd like to create root files to store
        //trees of pre-selected events and histograms on disk>
        TFile *TreeFile, *TreeFile_original, *TreeFile_selection, *HistoFile;      
        MapType_TTree Tree, Tree_original, Tree_selection;
        TTree *tree, *tree_bkg;
        
      
};


#endif
