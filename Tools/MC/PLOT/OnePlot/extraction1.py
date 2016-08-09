#!/usr/bin/python

##############
## Make Plots
##############

########################
# Libraries and Includes
########################

# Python System Modules
import os,sys,glob
import logging

import array
# Python Paths
ROOTLIB, PATH_OnePlot_MODULE1, PATH_OnePlot_MODULE2 = "", "", ""

if os.getenv("ROOTLIB"): ROOTLIB=os.getenv("ROOTLIB")
if os.getenv("PATH_OnePlot_MODULE1"): PATH_OnePlot_MODULE1=os.getenv("PATH_OnePlot_MODULE1")
if os.getenv("PATH_OnePlot_MODULE2"): PATH_OnePlot_MODULE2=os.getenv("PATH_OnePlot_MODULE2")

sys.path.append(ROOTLIB)
sys.path.append(PATH_OnePlot_MODULE1)
sys.path.append(PATH_OnePlot_MODULE2)

# Python Math Modules
import array
from math import sqrt,fabs,sin,pow

# ROOT Modules
from ROOT import TFile,TTree,TChain,TBranch,TH1F,TH2F
from ROOT import TCanvas,TPad
from ROOT import TLorentzVector
from ROOT import gROOT, gDirectory,gStyle
from ROOT import TPaveText, TLegend
import ROOT

# user modules
from oneplot import oneplot
from module_syst import systPDFCT10, systEnvelope, syst1v1, systCombine
from module_style import atlas_style
    
###############################
# User defined drawing function
###############################

if __name__ == "__main__":

    ###############################################
    # Simple Example for how to use OnePlot modules
    ###############################################
    
    # specify file name
    file="./output.root"
    Channel=["eeee","eemm","mmmm"]
    for channel in Channel:
        #Cut=["4Leptons","2Pairs","Mll1","Mll2","NumJets"]
        Cut=["4Leptons"]
        for cut in Cut:
            Lep_Var=["Pt","Eta","Phi"]
            Bin={"Pt" :range(0,241000,15000),"Eta":[x * 0.1 for x in range(-30,31,4)],"Phi":[x * 0.01 for x in range(-350,351,35) ] }
            for var in Lep_Var:
                for i in range(1,5):
                    hist1="sig"+ "_" + channel + "_" + cut + "_" + var + str(i)
                    hist2="bg"+ "_" + channel + "_" + cut + "_" + var + str(i)
                    # specify the histogram name
                    names=[hist1,hist2]
                    # set drawing options
                    legends=["VBS ZZ","QCD ZZ"]
                    #options=["LPE", "HIST"]
                    options=["HIST", "HIST"]
                    options_ratio=["LPE", "LPE"]
                    opt_legends=["F", "F"]
                    marker_types=[22, 0]            
                    marker_sizes=line_sizes=[2, 2]  
                    fill_colors=marker_colors=line_colors=[2,4] 
                    fill_types=[3001, 3001]
                    binning=Bin[var]
                    # set fig name
                    figname=channel + "_" + cut + "_" + var + str(i) 
                    # draw the histograms                          
                    theone = oneplot()                                                                                 
                    theone.initialize(file=file,names=names,legends=legends, opt_legends=opt_legends, 
                        options = options,                 
                        ratio=0.2,
                        options_ratio = options_ratio,
                        xtitle="leptons_"+var,ytitle="Events",
                        figname=figname, binning=binning,  
                        ratiotitle="VBS/QCD+1", 
                        marker_types=marker_types, marker_sizes=marker_sizes, marker_colors=marker_colors, 
                        fill_colors=fill_colors,
                        fill_types=fill_types,
                        line_sizes=line_sizes, line_colors=line_colors,
                        ratio_denominater=hist2)                                    
                    theone.plot1DHistogram()           
                    theone.finish()                         

                    #end


