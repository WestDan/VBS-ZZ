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
    file="./examplev1.root"
    
    # specify the histogram name
    names=["histo1", "histo2"]
    
    # set drawing options
    legends=["Test Histo 1", "Test Histo 2"]
    options=["LPE", "HIST"]
    opt_legends=["LPE", "L"]
    marker_types=[22, 0]
    marker_sizes=line_sizes=[2, 2]
    marker_colors=line_colors=[2,4]
    binning=[0,10,20,30,40,50,60,70,80,90,100]
    
    # set fig name
    figname="examplev1"
    
    # draw the histograms
    theone = oneplot()
    theone.initialize(file=file,names=names,legends=legends, opt_legends=opt_legends,
                    ratio=0.2,xtitle="Arbitrary",ytitle="Arbitrary",
                    figname=figname, binning=binning,
                    ratiotitle="Histo2/Histo1",
                    marker_types=marker_types, marker_sizes=marker_sizes, marker_colors=marker_colors,
                    line_sizes=line_sizes, line_colors=line_colors)
    theone.plot1DHistogram()
    theone.finish()    
    
    # end
    
    
    
    
    
    
