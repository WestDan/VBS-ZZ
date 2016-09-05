#include <string>
#include <vector>
#include <math.h>
#include <iostream>

using namespace std;

#include "xAODMuon/Muon.h"
#include "xAODEgamma/Electron.h"
#include "MyAnalysis/AnalysisVar.h"

void SetWeight(MapType2_Double& map, string cut, string chn, float value) {

    if(chn=="All" || chn=="ALL" || chn=="all") {
      MapType2_Double::iterator it;
      for(it=map.begin(); it!=map.end(); it++) {
        string chn1=(*it).first;
        map[chn1][cut] *= value;
      }
    }
    else {
       map[chn][cut] *= value;
    }
}


void SetFlag(MapType3_Int& map, string cut, string chn, int value, string opt="default") {

    if(chn=="All" || chn=="ALL" || chn=="all") {
      MapType3_Int::iterator it;
      for(it=map.begin(); it!=map.end(); it++) {
        string chn1=(*it).first;
        map[chn1][cut][opt]=value;
      }
    }
    else {
       map[chn][cut][opt]=value;
    }
}

void InitStrVec2(vector<string>& out, string in, string de) {
    int pos=0, pos_pre=0;
    while(true) {
        pos=in.find(de,pos_pre);
        if(pos==-1) {out.push_back(in.substr(pos_pre,in.size()-pos_pre)); break;}
        else  out.push_back(in.substr(pos_pre,pos-pos_pre));
        pos_pre=pos+1;
    }
}


void DoCounting(string sysname, MapType3_Counting& CNT, string chn, string cut) {
    double num = 1.0;
    CNT[sysname][chn][cut].num += num;
    CNT[sysname][chn][cut].err = sqrt(pow(CNT[sysname][chn][cut].err,2)+num);
}


void DoCounting(string sysname, MapType4_Counting& CNT, string chn, string cut, double w=1.0, string var="default") {
    double num = 1.0;
    double wnum = w*1.0;
    CNT[sysname][chn][cut][var].num += num;
    CNT[sysname][chn][cut][var].wnum += wnum;
    CNT[sysname][chn][cut][var].err = sqrt(pow(CNT[sysname][chn][cut][var].err,2)+wnum);
}


void CountMuObj(xAOD::Muon* muon, MapType_VString STEP_obj, MapType3_Counting& CNT, string sysname, bool& passAll) {

    vector<string> objstr = STEP_obj["mu"];
    
    for(int i=0; i<(int)objstr.size(); i++) {
      if(objstr[i].find("Pt_F") != string::npos) break;
      passAll = passAll && muon->auxdata< char >( objstr[i].c_str() );
      if(passAll) DoCounting(sysname, CNT, "mu", objstr[i]); 
    }
}

void CountEleObj(xAOD::Electron* electron, MapType_VString STEP_obj, MapType3_Counting& CNT, string sysname, bool& passAll) {

    vector<string> objstr = STEP_obj["ele"];

    for(int i=0; i<(int)objstr.size(); i++) {
      if(objstr[i].find("Pt_F") != string::npos) break;
      passAll = passAll && electron->auxdata< char >( objstr[i].c_str() ) ;
      if(passAll) DoCounting(sysname, CNT, "ele", objstr[i]);
    }
}

void CountJetObj(xAOD::Jet* jet, MapType_VString STEP_obj, MapType3_Counting& CNT, string sysname, bool& passAll) {

    vector<string> objstr = STEP_obj["jet"];
    for(int i=0; i<(int)objstr.size(); i++) {
      if(objstr[i].find("Clean") != string::npos) break;
      passAll = passAll && jet->auxdata< char >( objstr[i].c_str() );
      if(passAll) DoCounting(sysname, CNT, "jet", objstr[i]);
    }
}

void CountEvt(string sysname, vector<string> CHN, vector<string> STEP_cut, MapType3_Int& FLAG_cut_temp, MapType3_Int& FLAG_cut, MapType4_Counting& CNT, MapType2_Double& Evt_Weight) {

    MapType_VString::iterator it;
    bool passAll;
    for(int i=0; i<(int)CHN.size(); i++) {
      string chl =CHN[i];
      
      passAll=true;
      for(int j=0; j<(int)STEP_cut.size(); j++) {
        string cut=STEP_cut[j];
        passAll = passAll && FLAG_cut_temp[chl][cut]["default"];
        if(passAll) {
          FLAG_cut[chl][cut]["default"]=1;
          double w = Evt_Weight[chl][cut];
          DoCounting(sysname, CNT, chl, cut, w);
        }
      }
    }
}

void CountEvt(string sysname, vector<string> CHN, vector<string> STEP_cut, MapType3_Int& FLAG_cut_temp, MapType3_Int& FLAG_cut, MapType4_Counting& CNT, MapType2_Double& Evt_Weight, vector<string> Var_Map) {

    MapType_VString::iterator it;
    for(int i=0; i<(int)CHN.size(); i++) {
      string chl =CHN[i];

      for(int n=0; n<(int)Var_Map.size(); n++) {
        string fullvar=Var_Map[n];

        vector<string> splitvar;
        InitStrVec2(splitvar, fullvar, ";");           

        bool passAll=true;
        for(int j=0; j<(int)STEP_cut.size(); j++) {
          string cut=STEP_cut[j];
        
          if(fullvar.find(cut) == string::npos) {
            passAll = passAll && FLAG_cut_temp[chl][cut]["default"];
            if(passAll) {
              FLAG_cut[chl][cut][fullvar]=1;
              double w = Evt_Weight[chl][cut];
              DoCounting(sysname, CNT, chl, cut, w, fullvar);
            }
          }else {

            string var=splitvar[0];
            for(int k=0; k<(int)splitvar.size(); k++) {
              var=splitvar[k];
              if(var.find(cut) != string::npos) break;
            }
              
            passAll = passAll && FLAG_cut_temp[chl][cut][var];
            if(passAll) {
              FLAG_cut[chl][cut][fullvar]=1;
              double w = Evt_Weight[chl][cut];
              DoCounting(sysname, CNT, chl, cut, w, fullvar);
            }
          }
        }
      }
    }
}
