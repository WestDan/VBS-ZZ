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
from ROOT import TFile,TTree,TChain,TBranch,TH1F,TH2F,TH1D
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

    # plot for QBH v.s. ttbar: files & histo names
    files = [
        "../rootfiles/hist_all/ggH300.root",
        "../rootfiles/hist_all/ggH600.root",
        "../rootfiles/hist_all/ggH1000.root",
        "../rootfiles/hist_all/qqZZ.root",
        "../rootfiles/hist_all/Z.root"
        ]
        
    #hname = "NOMINAL_ee_MetCut120_DLepR100_DMetZ0_FracPt100_Njets100_Nbjets100_dphijet25met0_Metsig0_dMetZPhi"
    hname = "NOMINAL_ee_MetCut120_DLepR1.7_DMetZ0_FracPt0.2_Njets100_Nbjets0_dphijet25met0.7_Metsig0_dMetZPhi"
    histos = []
    names = [] 
    for file in files:
            fin_temp = TFile(file)
            hist_temp = fin_temp.Get(hname)
            hist_temp.SetDirectory(0)
            hist_temp.Scale(1./hist_temp.Integral())
            histos.append(hist_temp)
            names.append(file.split("/")[-1].split(".")[0])

    # set up plotting attributes
    legends=["ggH300", "ggH600", "ggH1000", "qqZZ", "Z+jets"]
    options=["HISTL", "HISTL", "HISTL", "HISTL", "HISTL"]
    opt_legends=["L", "L", "L", "L", "L"]
    #marker_types=[22, 0]
    #marker_sizes=line_sizes=[2, 2]
    #marker_colors=line_colors=[2,4]
    line_sizes=[2,2,2,2,2]
    line_colors=[2,3,4,5,6]
    binning=[0., 1., 1.2, 1.4, 1.6, 1.8, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.1, 3.2, 3.3, 3.4, 3.5]

    # set fig name
    figname="test_ee_post"

    # draw the histograms
    theone = oneplot()
    theone.initialize(#file=file,names=names,legends=legends, opt_legends=opt_legends,
                    list_histo=histos,names=names,legends=legends, opt_legends=opt_legends,
                    ratio=0.2,xtitle="dPhi(Z,MET)",ytitle="Arbitrary",
                    figname=figname, binning=binning,
                    ratiotitle="Ratio",
                    line_sizes=line_sizes, line_colors=line_colors)
    theone.plot1DHistogram()
    theone.finish()

        
