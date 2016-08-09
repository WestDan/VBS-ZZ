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

    ############################
    # Plot for Unfolding Results
    ############################
    
    files=[
        "/atlas/data18a/wyusheng/atlas/analysis/workdir/2014-WorkDir/20140403-ZZUnfolding/Aug1-Unfolding/unfolding/20141009-EWKUpdate/M4lData_FV/Results/ZZIncl8TeV_Unfolding.root",
        "/atlas/data18a/wyusheng/atlas/analysis/workdir/2014-WorkDir/20140403-ZZUnfolding/Aug1-Unfolding/unfolding/20141009-EWKUpdate/M4lData_PS/Results/ZZIncl8TeV_Unfolding.root",
        "/atlas/data18a/wyusheng/atlas/analysis/workdir/2014-WorkDir/20140403-ZZUnfolding/Aug1-Unfolding/unfolding/20141009-EWKUpdate/MZ2Data_FV/Results/ZZIncl8TeV_Unfolding.root",
        "/atlas/data18a/wyusheng/atlas/analysis/workdir/2014-WorkDir/20140403-ZZUnfolding/Aug1-Unfolding/unfolding/20141009-EWKUpdate/MZ2Data_PS/Results/ZZIncl8TeV_Unfolding.root",
        "/atlas/data18a/wyusheng/atlas/analysis/workdir/2014-WorkDir/20140403-ZZUnfolding/Aug1-Unfolding/unfolding/20141009-EWKUpdate/PtZZData_FV/Results/ZZIncl8TeV_Unfolding.root",
        "/atlas/data18a/wyusheng/atlas/analysis/workdir/2014-WorkDir/20140403-ZZUnfolding/Aug1-Unfolding/unfolding/20141009-EWKUpdate/PtZZData_PS/Results/ZZIncl8TeV_Unfolding.root",
    ]
    
    keynames=[
        "MZZ_FiducialVolume",
        "MZZ_PhaseSpace",
        "MZ2_FiducialVolume",
        "MZ2_PhaseSpace",        
        "PtZZ_FiducialVolume",
        "PtZZ_PhaseSpace",    
    ]
    
    xnames=[
        "M_{4l}",
        "M_{4l}",
        "M_{34}",
        "M_{34}",
        "p_{T}^{4l}",
        "p_{T}^{4l}",
    ]
    
    for findex in range(len(files)):

        # open file
        file= files[findex]
        keyname=keynames[findex]
        xname=xnames[findex]
        fin = TFile(file)
        
        # load basic histograms
        histo_unfold = fin.Get("Final/Final_Nominal_Unfolded_Normalized")
        histo_fulluncert = fin.Get("Final/Final_CombinedUncertainty")
        histo_mc = fin.Get("Final/Final_MCTruth_Normalized")
        
        # load the histograms to calculate MC uncertainty
        histo_mc_psup = fin.Get("Systematic/Systematic_PSUp_MCTruth")
        histo_mc_psdn = fin.Get("Systematic/Systematic_PSDown_MCTruth")
        histo_mc_scaleup = fin.Get("Systematic/Systematic_ScaleUp_MCTruth") 
        histo_mc_scaledn = fin.Get("Systematic/Systematic_ScaleDown_MCTruth") 
        histo_mc_pdfup = fin.Get("Systematic/Systematic_PDFUp_MCTruth")
        histo_mc_pdfdn = fin.Get("Systematic/Systematic_PDFDown_MCTruth")
        histo_mc_ewkup = fin.Get("Systematic/Systematic_EWCorrUp_MCTruth")
        histo_mc_ewkdn = fin.Get("Systematic/Systematic_EWCorrDown_MCTruth")
        
        dict_ps = systEnvelope([histo_mc,histo_mc_psup, histo_mc_psdn])
        dict_scale = systEnvelope([histo_mc,histo_mc_scaleup, histo_mc_scaledn])
        dict_pdf = systEnvelope([histo_mc,histo_mc_pdfup, histo_mc_pdfdn])
        dict_ewk = systEnvelope([histo_mc,histo_mc_ewkup, histo_mc_ewkdn])
        dict_all = systCombine(list_histo_sysup=[dict_ps["NOM"], dict_ps["UP"], dict_scale["UP"], dict_pdf["UP"], dict_ewk["UP"]], 
                            list_histo_sysdn=[dict_ps["NOM"], dict_ps["DN"], dict_scale["DN"], dict_pdf["DN"], dict_ewk["DN"]])
        
        # get final histograms
        nbins=histo_unfold.GetNbinsX()
        histo_unfold_final = histo_unfold.Clone("Unfolded Final+Sys")
        histo_mc_final = histo_mc.Clone("Truth Final+Sys")
        for index in range(nbins):
            histo_unfold.SetBinError(index+1, 0.00001)
            histo_mc.SetBinError(index+1, 0.00001)
            histo_unfold_final.SetBinError(index+1, histo_fulluncert.GetBinContent(index+1))
            histo_mc_final.SetBinError(index+1, fabs(dict_all["UP"].GetBinContent(index+1)-dict_all["NOM"].GetBinContent(index+1)))
        
        # Draw histograms
        names=["Unfolded Data", "Unfolded Uncertainty", "MC Truth", "Truth Uncertainty"]
        legends=["Unfolded Data", "Uncert. (Unfolded Data)", "MC Truth", "Uncert. (MC Truth)"]
        histos=[histo_unfold, histo_unfold_final, histo_mc, histo_mc_final]
        figname="final_unfolded_{0}".format(keyname)
        options=["LPE", "E2", "LPE", "E2"]
        opt_legends=["LP", "F", "LP", "F"]
        marker_types=[26, 26, 32, 32]
        marker_sizes=[1, 1, 1, 1]
        line_sizes=[2, 2, 2, 2]
        marker_colors=line_colors=fill_colors=[4, 4, 2, 2]
        fill_types=[0, 3001, 0, 3001]
        
        # call one plot
        theone = oneplot()
        theone._figtype=["png", "pdf"]
        theone.initialize(list_histo=histos,names=names,legends=legends, opt_legends=opt_legends,
                        ratio=0.2,xtitle="{0} [GeV]".format(xname),ytitle="#Delta#sigma / Bin [fb]",
                        figname=figname, ratiotitle="Data / MC", ratio_denominater="MC Truth",
                        equalbin=True, options=options, marker_types=marker_types, marker_colors=marker_colors, marker_sizes=marker_sizes,
                        line_colors=line_colors, line_sizes=line_sizes,
                        fill_colors=fill_colors, fill_types=fill_types)
        theone.plot1DHistogram()
        theone.finish()
        
        # plot test 1 (plot all truth distribution together)
        #names=["nom", "PSUp", "PSDown", "ScaleUp", "ScaleDown", "PDFUp", "PDFDown", "EWKUp", "EWKDown"]
        #legends=["nom", "PSUp", "PSDown", "ScaleUp", "ScaleDown", "PDFUp", "PDFDown", "EWKUp", "EWKDown"]
        #histos=[histo_mc, histo_mc_psup, histo_mc_psdn, histo_mc_scaleup, histo_mc_scaledn, histo_mc_pdfup, histo_mc_pdfdn, histo_mc_ewkup, histo_mc_ewkdn]
        #figname="compare_mctruth_{0}".format(keyname)
        #options=["HIST", "LPE", "LPE", "LPE", "LPE", "LPE", "LPE", "LPE", "LPE"]
        #marker_types=[1, 1, 1, 1, 1, 1, 1, 1, 1]
        #line_colors=[1, 2, 2, 3, 3, 4, 4, 5, 5]
        #line_sizes=[2, 2, 2, 2, 2, 2, 2, 2, 2]
        #theone = oneplot()
        #theone.initialize(list_histo=histos,names=names,legends=legends,
        #                ratio=0.2,xtitle="{0} [GeV]".format(xname),ytitle="#Delta#sigma / Bin [fb]",
        #                figname=figname, ratiotitle="Ratio",
        #                equalbin=True, options=options, line_colors=line_colors, line_sizes=line_sizes)
        #theone.plot1DHistogram()
        #theone.finish()
        
        fin.Close()
        
        # save the files into histograms
        #fout = TFile(figname+".root", "recreate")
        #histo_unfold_final.Write()
        #histo_mc_final.Write()
        #fout.Close()
