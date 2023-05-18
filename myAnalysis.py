#!/usr/bin/env python
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True


def get_vec(obj):
  return ROOT.Math.PtEtaPhiMVector(obj.pt, obj.eta, obj.phi, obj.mass)


def dR(eta1, eta2, phi1, phi2):
    dEta = abs(eta1 - eta2)
    dPhi = abs(phi1 - phi2)  # abs() is fine because delta_phi will be squared eventually 
    # ensure that delta_phi is between 0 and pi (not -pi and pi due to abs())
    if dPhi > ROOT.TMath.Pi(): dPhi -= ROOT.TMath.Pi()
    return ROOT.TMath.Sqrt(dEta**2 + dPhi**2)


class MyAnalysis(Module):
    def __init__(self, datamc, lumi=1.0, dict_xs=None, dict_ngen=None):
        self.writeHistFile = True
        self.lumi = lumi
        self.dict_xs = dict_xs
        self.dict_ngen = dict_ngen
        self.datamc = datamc

    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)
        
        self.twoprong_chargedIso_iso_sym = ROOT.TH1F('twoprong_chargedIso_iso_sym', '; Twoprong_chargedIso', 2000, 0, 100)
        self.twoprong_chargedIso_iso_asym = ROOT.TH1F('twoprong_chargedIso_iso_asym', '; Twoprong_chargedIso', 2000, 0, 100)
        self.twoprong_chargedIso_noniso_sym = ROOT.TH1F('twoprong_chargedIso_noniso_sym', '; Twoprong_chargedIso', 2000, 0, 100)
        self.twoprong_chargedIso_noniso_asym = ROOT.TH1F('twoprong_chargedIso_noniso_asym', '; Twoprong_chargedIso', 2000, 0, 100)
        self.addObject(self.twoprong_chargedIso_iso_sym)
        self.addObject(self.twoprong_chargedIso_iso_asym)
        self.addObject(self.twoprong_chargedIso_noniso_sym)
        self.addObject(self.twoprong_chargedIso_noniso_asym)

        self.twoprong_phi_noniso_asym = ROOT.TH1F('twoprong_phi_noniso_asym', ';Twoprong_phi', 70, -3.5, 3.5)
        self.addObject(self.twoprong_phi_noniso_asym)
        self.twoprong_phi_noniso_sym = ROOT.TH1F('twoprong_phi_noniso_sym', ';Twoprong_phi', 70, -3.5, 3.5)
        self.addObject(self.twoprong_phi_noniso_sym)
        self.twoprong_phi_iso_asym = ROOT.TH1F('twoprong_phi_iso_asym', ';Twoprong_phi', 70, -3.5, 3.5)
        self.addObject(self.twoprong_phi_iso_asym)
        self.twoprong_phi_iso_sym = ROOT.TH1F('twoprong_phi_iso_sym', ';Twoprong_phi', 70, -3.5, 3.5)
        self.addObject(self.twoprong_phi_iso_sym)

        self.twoprong_pt_noniso_asym = ROOT.TH1F('twoprong_pt_noniso_asym', '; Twoprong_pt', 150, 0, 1500)  # loose-loose twoprong
        self.addObject(self.twoprong_pt_noniso_asym)
        self.twoprong_pt_noniso_sym = ROOT.TH1F('twoprong_pt_noniso_sym', '; Twoprong_pt', 150, 0, 1500)  # loose-tight twoprong
        self.addObject(self.twoprong_pt_noniso_sym)
        self.twoprong_pt_iso_asym = ROOT.TH1F('twoprong_pt_iso_asym', '; Twoprong_pt', 150, 0, 1500)  # tight-loose twoprong
        self.addObject(self.twoprong_pt_iso_asym)
        self.twoprong_pt_iso_sym = ROOT.TH1F('twoprong_pt_iso_sym', '; Twoprong_pt', 150, 0, 1500)  # tight-tight twoprong
        self.addObject(self.twoprong_pt_iso_sym)
        
        self.twoprong_eta_noniso_asym = ROOT.TH1F('twoprong_eta_noniso_asym', '; Twoprong_Eta', 100, -5, 5)  # loose-loose twoprong
        self.addObject(self.twoprong_eta_noniso_asym)
        self.twoprong_eta_noniso_sym = ROOT.TH1F('twoprong_eta_noniso_sym', '; Twoprong_Eta', 100, -5, 5)  # loose-tight twoprong
        self.addObject(self.twoprong_eta_noniso_sym)
        self.twoprong_eta_iso_asym = ROOT.TH1F('twoprong_eta_iso_asym', '; Twoprong_Eta', 100, -5, 5)  # tight-loose twoprong
        self.addObject(self.twoprong_eta_iso_asym)
        self.twoprong_eta_iso_sym = ROOT.TH1F('twoprong_eta_iso_sym', '; Twoprong_Eta', 100, -5, 5)  # tight-tight twoprong
        self.addObject(self.twoprong_eta_iso_sym)
        
        self.ntwoprong_regions = ROOT.TH1F('ntwoprong_regions', 'Twoporng type', 4, 0, 4)
        self.addObject(self.ntwoprong_regions)

        self.twoprong_masspi0_noniso_asym = ROOT.TH1F('twoprong_masspi0_noniso_asym', '; Twoprong_massPi0', 300, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_noniso_asym)
        self.twoprong_masspi0_noniso_sym = ROOT.TH1F('twoprong_masspi0_noniso_sym', '; Twoprong_massPi0', 300, 0, 50)  # loose-tight twoprong
        self.addObject(self.twoprong_masspi0_noniso_sym)
        self.twoprong_masspi0_iso_asym = ROOT.TH1F('twoprong_masspi0_iso_asym', '; Twoprong_massPi0', 300, 0, 50)  # tight-loose twoprong
        self.addObject(self.twoprong_masspi0_iso_asym)
        self.twoprong_masspi0_iso_sym = ROOT.TH1F('twoprong_masspi0_iso_sym', '; Twoprong_massPi0', 300, 0, 50)  # tight-tight twoprong
        self.addObject(self.twoprong_masspi0_iso_sym)
        
        self.twoprong_masspi0_noniso_asym_barrel = ROOT.TH1F('twoprong_masspi0_noniso_asym_barrel', '; Twoprong_massPi0', 300, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel)
        self.twoprong_masspi0_noniso_sym_barrel = ROOT.TH1F('twoprong_masspi0_noniso_sym_barrel', '; Twoprong_massPi0', 300, 0, 50)  # loose-tight twoprong
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel)
        self.twoprong_masspi0_iso_asym_barrel = ROOT.TH1F('twoprong_masspi0_iso_asym_barrel', '; Twoprong_massPi0', 300, 0, 50)  # tight-loose twoprong
        self.addObject(self.twoprong_masspi0_iso_asym_barrel)
        self.twoprong_masspi0_iso_sym_barrel = ROOT.TH1F('twoprong_masspi0_iso_sym_barrel', '; Twoprong_massPi0', 300, 0, 50)  # tight-tight twoprong
        self.addObject(self.twoprong_masspi0_iso_sym_barrel)
        
        self.twoprong_masspi0_noniso_asym_endcap = ROOT.TH1F('twoprong_masspi0_noniso_asym_endcap', '; Twoprong_massPi0', 300, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap)
        self.twoprong_masspi0_noniso_sym_endcap = ROOT.TH1F('twoprong_masspi0_noniso_sym_endcap', '; Twoprong_massPi0', 300, 0, 50)  # loose-tight twoprong
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap)
        self.twoprong_masspi0_iso_asym_endcap = ROOT.TH1F('twoprong_masspi0_iso_asym_endcap', '; Twoprong_massPi0', 300, 0, 50)  # tight-loose twoprong
        self.addObject(self.twoprong_masspi0_iso_asym_endcap)
        self.twoprong_masspi0_iso_sym_endcap = ROOT.TH1F('twoprong_masspi0_iso_sym_endcap', '; Twoprong_massPi0', 300, 0, 50)  # tight-tight twoprong
        self.addObject(self.twoprong_masspi0_iso_sym_endcap)
        
        self.hthat_gjets = ROOT.TH1F('GJETS_hthat_lhe', '; htHat_lhe', 150, 0, 1500)
        self.addObject(self.hthat_gjets)

        # TwoProng Pt Binning
        self.pt_bins = [0,20,40,60,80,120,160,240,300,460,700]
        
        masspi0_iso_sym_barrel_hists = []
        masspi0_iso_sym_endcap_hists = []
        masspi0_iso_asym_barrel_hists = []
        masspi0_iso_asym_endcap_hists = []
        masspi0_noniso_sym_barrel_hists = []
        masspi0_noniso_sym_endcap_hists = []
        masspi0_noniso_asym_barrel_hists = []
        masspi0_noniso_asym_endcap_hists = []

        bin0 = self.pt_bins[0] 
        bin1 = self.pt_bins[1] 
        bin2 = self.pt_bins[2] 
        bin3 = self.pt_bins[3] 
        bin4 = self.pt_bins[4] 
        bin5 = self.pt_bins[5] 
        bin6 = self.pt_bins[6] 
        bin7 = self.pt_bins[7] 
        bin8 = self.pt_bins[8] 
        bin9 = self.pt_bins[9] 
        bin10 = self.pt_bins[10] 

        # bin0 - bin1
        self.twoprong_masspi0_noniso_asym_barrel_bin0 = ROOT.TH1F('twoprong_masspi0_noniso_asym_barrel_'+str(bin0)+'_'+str(bin1), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin0)
        self.twoprong_masspi0_noniso_asym_endcap_bin0 = ROOT.TH1F('twoprong_masspi0_noniso_asym_endcap_'+str(bin0)+'_'+str(bin1), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin0)
        self.twoprong_masspi0_iso_sym_barrel_bin0 = ROOT.TH1F('twoprong_masspi0_iso_sym_barrel_'+str(bin0)+'_'+str(bin1), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin0)
        self.twoprong_masspi0_iso_sym_endcap_bin0 = ROOT.TH1F('twoprong_masspi0_iso_sym_endcap_'+str(bin0)+'_'+str(bin1), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin0)
        self.twoprong_masspi0_iso_asym_barrel_bin0 = ROOT.TH1F('twoprong_masspi0_iso_asym_barrel_'+str(bin0)+'_'+str(bin1), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin0)
        self.twoprong_masspi0_iso_asym_endcap_bin0 = ROOT.TH1F('twoprong_masspi0_iso_asym_endcap_'+str(bin0)+'_'+str(bin1), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin0)
        self.twoprong_masspi0_noniso_sym_barrel_bin0 = ROOT.TH1F('twoprong_masspi0_noniso_sym_barrel_'+str(bin0)+'_'+str(bin1), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin0)
        self.twoprong_masspi0_noniso_sym_endcap_bin0 = ROOT.TH1F('twoprong_masspi0_noniso_sym_endcap_'+str(bin0)+'_'+str(bin1), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin0)

        # bin1 - bin2
        self.twoprong_masspi0_noniso_asym_barrel_bin1 = ROOT.TH1F('twoprong_masspi0_noniso_asym_barrel_'+str(bin1)+'_'+str(bin2), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin1)
        self.twoprong_masspi0_noniso_asym_endcap_bin1 = ROOT.TH1F('twoprong_masspi0_noniso_asym_endcap_'+str(bin1)+'_'+str(bin2), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin1)
        self.twoprong_masspi0_iso_sym_barrel_bin1 = ROOT.TH1F('twoprong_masspi0_iso_sym_barrel_'+str(bin1)+'_'+str(bin2), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin1)
        self.twoprong_masspi0_iso_sym_endcap_bin1 = ROOT.TH1F('twoprong_masspi0_iso_sym_endcap_'+str(bin1)+'_'+str(bin2), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin1)
        self.twoprong_masspi0_iso_asym_barrel_bin1 = ROOT.TH1F('twoprong_masspi0_iso_asym_barrel_'+str(bin1)+'_'+str(bin2), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin1)
        self.twoprong_masspi0_iso_asym_endcap_bin1 = ROOT.TH1F('twoprong_masspi0_iso_asym_endcap_'+str(bin1)+'_'+str(bin2), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin1)
        self.twoprong_masspi0_noniso_sym_barrel_bin1 = ROOT.TH1F('twoprong_masspi0_noniso_sym_barrel_'+str(bin1)+'_'+str(bin2), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin1)
        self.twoprong_masspi0_noniso_sym_endcap_bin1 = ROOT.TH1F('twoprong_masspi0_noniso_sym_endcap_'+str(bin1)+'_'+str(bin2), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin1)
        
        # bin2 - bin3
        self.twoprong_masspi0_noniso_asym_barrel_bin2 = ROOT.TH1F('twoprong_masspi0_noniso_asym_barrel_'+str(bin2)+'_'+str(bin3), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin2)
        self.twoprong_masspi0_noniso_asym_endcap_bin2 = ROOT.TH1F('twoprong_masspi0_noniso_asym_endcap_'+str(bin2)+'_'+str(bin3), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin2)
        self.twoprong_masspi0_iso_sym_barrel_bin2 = ROOT.TH1F('twoprong_masspi0_iso_sym_barrel_'+str(bin2)+'_'+str(bin3), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin2)
        self.twoprong_masspi0_iso_sym_endcap_bin2 = ROOT.TH1F('twoprong_masspi0_iso_sym_endcap_'+str(bin2)+'_'+str(bin3), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin2)
        self.twoprong_masspi0_iso_asym_barrel_bin2 = ROOT.TH1F('twoprong_masspi0_iso_asym_barrel_'+str(bin2)+'_'+str(bin3), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin2)
        self.twoprong_masspi0_iso_asym_endcap_bin2 = ROOT.TH1F('twoprong_masspi0_iso_asym_endcap_'+str(bin2)+'_'+str(bin3), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin2)
        self.twoprong_masspi0_noniso_sym_barrel_bin2 = ROOT.TH1F('twoprong_masspi0_noniso_sym_barrel_'+str(bin2)+'_'+str(bin3), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin2)
        self.twoprong_masspi0_noniso_sym_endcap_bin2 = ROOT.TH1F('twoprong_masspi0_noniso_sym_endcap_'+str(bin2)+'_'+str(bin3), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin2)

        # bin3 - bin4
        self.twoprong_masspi0_noniso_asym_barrel_bin3 = ROOT.TH1F('twoprong_masspi0_noniso_asym_barrel_'+str(bin3)+'_'+str(bin4), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin3)
        self.twoprong_masspi0_noniso_asym_endcap_bin3 = ROOT.TH1F('twoprong_masspi0_noniso_asym_endcap_'+str(bin3)+'_'+str(bin4), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin3)
        self.twoprong_masspi0_iso_sym_barrel_bin3 = ROOT.TH1F('twoprong_masspi0_iso_sym_barrel_'+str(bin3)+'_'+str(bin4), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin3)
        self.twoprong_masspi0_iso_sym_endcap_bin3 = ROOT.TH1F('twoprong_masspi0_iso_sym_endcap_'+str(bin3)+'_'+str(bin4), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin3)
        self.twoprong_masspi0_iso_asym_barrel_bin3 = ROOT.TH1F('twoprong_masspi0_iso_asym_barrel_'+str(bin3)+'_'+str(bin4), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin3)
        self.twoprong_masspi0_iso_asym_endcap_bin3 = ROOT.TH1F('twoprong_masspi0_iso_asym_endcap_'+str(bin3)+'_'+str(bin4), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin3)
        self.twoprong_masspi0_noniso_sym_barrel_bin3 = ROOT.TH1F('twoprong_masspi0_noniso_sym_barrel_'+str(bin3)+'_'+str(bin4), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin3)
        self.twoprong_masspi0_noniso_sym_endcap_bin3 = ROOT.TH1F('twoprong_masspi0_noniso_sym_endcap_'+str(bin3)+'_'+str(bin4), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin3)
        
        # bin4 - bin5
        self.twoprong_masspi0_noniso_asym_barrel_bin4 = ROOT.TH1F('twoprong_masspi0_noniso_asym_barrel_'+str(bin4)+'_'+str(bin5), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin4)
        self.twoprong_masspi0_noniso_asym_endcap_bin4 = ROOT.TH1F('twoprong_masspi0_noniso_asym_endcap_'+str(bin4)+'_'+str(bin5), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin4)
        self.twoprong_masspi0_iso_sym_barrel_bin4 = ROOT.TH1F('twoprong_masspi0_iso_sym_barrel_'+str(bin4)+'_'+str(bin5), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin4)
        self.twoprong_masspi0_iso_sym_endcap_bin4 = ROOT.TH1F('twoprong_masspi0_iso_sym_endcap_'+str(bin4)+'_'+str(bin5), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin4)
        self.twoprong_masspi0_iso_asym_barrel_bin4 = ROOT.TH1F('twoprong_masspi0_iso_asym_barrel_'+str(bin4)+'_'+str(bin5), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin4)
        self.twoprong_masspi0_iso_asym_endcap_bin4 = ROOT.TH1F('twoprong_masspi0_iso_asym_endcap_'+str(bin4)+'_'+str(bin5), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin4)
        self.twoprong_masspi0_noniso_sym_barrel_bin4 = ROOT.TH1F('twoprong_masspi0_noniso_sym_barrel_'+str(bin4)+'_'+str(bin5), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin4)
        self.twoprong_masspi0_noniso_sym_endcap_bin4 = ROOT.TH1F('twoprong_masspi0_noniso_sym_endcap_'+str(bin4)+'_'+str(bin5), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin4)

        # bin5 - bin6
        self.twoprong_masspi0_noniso_asym_barrel_bin5 = ROOT.TH1F('twoprong_masspi0_noniso_asym_barrel_'+str(bin5)+'_'+str(bin6), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin5)
        self.twoprong_masspi0_noniso_asym_endcap_bin5 = ROOT.TH1F('twoprong_masspi0_noniso_asym_endcap_'+str(bin5)+'_'+str(bin6), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin5)
        self.twoprong_masspi0_iso_sym_barrel_bin5 = ROOT.TH1F('twoprong_masspi0_iso_sym_barrel_'+str(bin5)+'_'+str(bin6), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin5)
        self.twoprong_masspi0_iso_sym_endcap_bin5 = ROOT.TH1F('twoprong_masspi0_iso_sym_endcap_'+str(bin5)+'_'+str(bin6), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin5)
        self.twoprong_masspi0_iso_asym_barrel_bin5 = ROOT.TH1F('twoprong_masspi0_iso_asym_barrel_'+str(bin5)+'_'+str(bin6), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin5)
        self.twoprong_masspi0_iso_asym_endcap_bin5 = ROOT.TH1F('twoprong_masspi0_iso_asym_endcap_'+str(bin5)+'_'+str(bin6), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin5)
        self.twoprong_masspi0_noniso_sym_barrel_bin5 = ROOT.TH1F('twoprong_masspi0_noniso_sym_barrel_'+str(bin5)+'_'+str(bin6), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin5)
        self.twoprong_masspi0_noniso_sym_endcap_bin5 = ROOT.TH1F('twoprong_masspi0_noniso_sym_endcap_'+str(bin5)+'_'+str(bin6), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin5)
    
        # bin6 - bin7
        self.twoprong_masspi0_noniso_asym_barrel_bin6 = ROOT.TH1F('twoprong_masspi0_noniso_asym_barrel_'+str(bin6)+'_'+str(bin7), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin6)
        self.twoprong_masspi0_noniso_asym_endcap_bin6 = ROOT.TH1F('twoprong_masspi0_noniso_asym_endcap_'+str(bin6)+'_'+str(bin7), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin6)
        self.twoprong_masspi0_iso_sym_barrel_bin6 = ROOT.TH1F('twoprong_masspi0_iso_sym_barrel_'+str(bin6)+'_'+str(bin7), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin6)
        self.twoprong_masspi0_iso_sym_endcap_bin6 = ROOT.TH1F('twoprong_masspi0_iso_sym_endcap_'+str(bin6)+'_'+str(bin7), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin6)
        self.twoprong_masspi0_iso_asym_barrel_bin6 = ROOT.TH1F('twoprong_masspi0_iso_asym_barrel_'+str(bin6)+'_'+str(bin7), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin6)
        self.twoprong_masspi0_iso_asym_endcap_bin6 = ROOT.TH1F('twoprong_masspi0_iso_asym_endcap_'+str(bin6)+'_'+str(bin7), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin6)
        self.twoprong_masspi0_noniso_sym_barrel_bin6 = ROOT.TH1F('twoprong_masspi0_noniso_sym_barrel_'+str(bin6)+'_'+str(bin7), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin6)
        self.twoprong_masspi0_noniso_sym_endcap_bin6 = ROOT.TH1F('twoprong_masspi0_noniso_sym_endcap_'+str(bin6)+'_'+str(bin7), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin6)
        
        # bin7 - bin8
        self.twoprong_masspi0_noniso_asym_barrel_bin7 = ROOT.TH1F('twoprong_masspi0_noniso_asym_barrel_'+str(bin7)+'_'+str(bin8), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin7)
        self.twoprong_masspi0_noniso_asym_endcap_bin7 = ROOT.TH1F('twoprong_masspi0_noniso_asym_endcap_'+str(bin7)+'_'+str(bin8), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin7)
        self.twoprong_masspi0_iso_sym_barrel_bin7 = ROOT.TH1F('twoprong_masspi0_iso_sym_barrel_'+str(bin7)+'_'+str(bin8), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin7)
        self.twoprong_masspi0_iso_sym_endcap_bin7 = ROOT.TH1F('twoprong_masspi0_iso_sym_endcap_'+str(bin7)+'_'+str(bin8), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin7)
        self.twoprong_masspi0_iso_asym_barrel_bin7 = ROOT.TH1F('twoprong_masspi0_iso_asym_barrel_'+str(bin7)+'_'+str(bin8), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin7)
        self.twoprong_masspi0_iso_asym_endcap_bin7 = ROOT.TH1F('twoprong_masspi0_iso_asym_endcap_'+str(bin7)+'_'+str(bin8), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin7)
        self.twoprong_masspi0_noniso_sym_barrel_bin7 = ROOT.TH1F('twoprong_masspi0_noniso_sym_barrel_'+str(bin7)+'_'+str(bin8), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin7)
        self.twoprong_masspi0_noniso_sym_endcap_bin7 = ROOT.TH1F('twoprong_masspi0_noniso_sym_endcap_'+str(bin7)+'_'+str(bin8), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin7)
        
        # bin8 - bin9
        self.twoprong_masspi0_noniso_asym_barrel_bin8 = ROOT.TH1F('twoprong_masspi0_noniso_asym_barrel_'+str(bin8)+'_'+str(bin9), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin8)
        self.twoprong_masspi0_noniso_asym_endcap_bin8 = ROOT.TH1F('twoprong_masspi0_noniso_asym_endcap_'+str(bin8)+'_'+str(bin9), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin8)
        self.twoprong_masspi0_iso_sym_barrel_bin8 = ROOT.TH1F('twoprong_masspi0_iso_sym_barrel_'+str(bin8)+'_'+str(bin9), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin8)
        self.twoprong_masspi0_iso_sym_endcap_bin8 = ROOT.TH1F('twoprong_masspi0_iso_sym_endcap_'+str(bin8)+'_'+str(bin9), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin8)
        self.twoprong_masspi0_iso_asym_barrel_bin8 = ROOT.TH1F('twoprong_masspi0_iso_asym_barrel_'+str(bin8)+'_'+str(bin9), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin8)
        self.twoprong_masspi0_iso_asym_endcap_bin8 = ROOT.TH1F('twoprong_masspi0_iso_asym_endcap_'+str(bin8)+'_'+str(bin9), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin8)
        self.twoprong_masspi0_noniso_sym_barrel_bin8 = ROOT.TH1F('twoprong_masspi0_noniso_sym_barrel_'+str(bin8)+'_'+str(bin9), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin8)
        self.twoprong_masspi0_noniso_sym_endcap_bin8 = ROOT.TH1F('twoprong_masspi0_noniso_sym_endcap_'+str(bin8)+'_'+str(bin9), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin8)
        
        # bin9 - bin10
        self.twoprong_masspi0_noniso_asym_barrel_bin9 = ROOT.TH1F('twoprong_masspi0_noniso_asym_barrel_'+str(bin9)+'_'+str(bin10), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin9)
        self.twoprong_masspi0_noniso_asym_endcap_bin9 = ROOT.TH1F('twoprong_masspi0_noniso_asym_endcap_'+str(bin9)+'_'+str(bin10), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin9)
        self.twoprong_masspi0_iso_sym_barrel_bin9 = ROOT.TH1F('twoprong_masspi0_iso_sym_barrel_'+str(bin9)+'_'+str(bin10), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin9)
        self.twoprong_masspi0_iso_sym_endcap_bin9 = ROOT.TH1F('twoprong_masspi0_iso_sym_endcap_'+str(bin9)+'_'+str(bin10), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin9)
        self.twoprong_masspi0_iso_asym_barrel_bin9 = ROOT.TH1F('twoprong_masspi0_iso_asym_barrel_'+str(bin9)+'_'+str(bin10), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin9)
        self.twoprong_masspi0_iso_asym_endcap_bin9 = ROOT.TH1F('twoprong_masspi0_iso_asym_endcap_'+str(bin9)+'_'+str(bin10), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin9)
        self.twoprong_masspi0_noniso_sym_barrel_bin9 = ROOT.TH1F('twoprong_masspi0_noniso_sym_barrel_'+str(bin9)+'_'+str(bin10), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin9)
        self.twoprong_masspi0_noniso_sym_endcap_bin9 = ROOT.TH1F('twoprong_masspi0_noniso_sym_endcap_'+str(bin9)+'_'+str(bin10), '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin9)
        
        # bin10 +
        self.twoprong_masspi0_noniso_asym_barrel_bin10 = ROOT.TH1F('twoprong_masspi0_noniso_asym_barrel_'+str(bin10)+'+', '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin10)
        self.twoprong_masspi0_noniso_asym_endcap_bin10 = ROOT.TH1F('twoprong_masspi0_noniso_asym_endcap_'+str(bin10)+'+', '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin10)
        self.twoprong_masspi0_iso_sym_barrel_bin10 = ROOT.TH1F('twoprong_masspi0_iso_sym_barrel_'+str(bin10)+'+', '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin10)
        self.twoprong_masspi0_iso_sym_endcap_bin10 = ROOT.TH1F('twoprong_masspi0_iso_sym_endcap_'+str(bin10)+'+', '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin10)
        self.twoprong_masspi0_iso_asym_barrel_bin10 = ROOT.TH1F('twoprong_masspi0_iso_asym_barrel_'+str(bin10)+'+', '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin10)
        self.twoprong_masspi0_iso_asym_endcap_bin10 = ROOT.TH1F('twoprong_masspi0_iso_asym_endcap_'+str(bin10)+'+', '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin10)
        self.twoprong_masspi0_noniso_sym_barrel_bin10 = ROOT.TH1F('twoprong_masspi0_noniso_sym_barrel_'+str(bin10)+'+', '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin10)
        self.twoprong_masspi0_noniso_sym_endcap_bin10 = ROOT.TH1F('twoprong_masspi0_noniso_sym_endcap_'+str(bin10)+'+', '; Twoprong_massPi0', 300, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin10)
   

    def analyze(self, event):
        if self.dict_xs and self.dict_ngen:
            dataset_id = event.dataset_id
            xs = self.dict_xs[dataset_id]
            Ngen = self.dict_ngen[dataset_id]
            weight = xs * self.lumi / Ngen
        else:
            weight = 1.0
    
        #hthat = Object(event, "htHat")
        twoprongs = Collection(event, "TwoProng")
        photons = Collection(event, "Photon")
        pass_trigger = event.HLT_Photon200

        # Immediately fill lhe histogram
        #if len(photons) > 0: self.hthat_gjets.Fill(hthat.lhe, weight)
        
        bins = []
        for i in range(len(self.pt_bins)):
            bins.append(self.pt_bins[i])
        

        sel_photons = []
        for photon in photons:
            if photon.cutBased >= 1 and photon.pt > 220 and abs(photon.scEta) < 1.4442: sel_photons.append(photon)
            break

        self.ntwoprong_regions.Fill(len(sel_photons))
        if len(sel_photons) == 0 or not pass_trigger: return True
        sel_photon = sel_photons[0]
        
        iso_sym_tp = []
        iso_asym_tp = []
        noniso_sym_tp = []
        noniso_asym_tp = []
        for twoprong in twoprongs:
            if twoprong.pt < 20 or dR(twoprong.eta, sel_photon.eta, twoprong.phi, sel_photon.phi) < 0.1 or twoprong.chargedIso > 1: continue
            if twoprong.passIso and twoprong.passSym: 
                iso_sym_tp.append(twoprong)
                break
            if twoprong.passIso and not twoprong.passSym: iso_asym_tp.append(twoprong)
            if not twoprong.passIso and twoprong.passSym: noniso_sym_tp.append(twoprong)
            if not twoprong.passIso and not twoprong.passSym: noniso_asym_tp.append(twoprong)
        
        if len(iso_sym_tp) != 0: sel_tp = iso_sym_tp[0]
        elif len(iso_asym_tp) != 0: sel_tp = iso_asym_tp[0]
        elif len(noniso_sym_tp) != 0: sel_tp = noniso_sym_tp[0]
        elif len(noniso_asym_tp) != 0: sel_tp = noniso_asym_tp[0]
        else: return True

        if len(iso_sym_tp) != 0: # tight tight
            self.twoprong_phi_iso_sym.Fill(sel_tp.phi, weight)
            self.twoprong_pt_iso_sym.Fill(sel_tp.pt, weight)
            self.twoprong_eta_iso_sym.Fill(sel_tp.eta, weight)
            self.ntwoprong_regions.Fill(0)
            self.twoprong_chargedIso_iso_sym.Fill(sel_tp.chargedIso, weight)
            self.twoprong_masspi0_iso_sym.Fill(sel_tp.massPi0, weight)
            if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_iso_sym_barrel.Fill(sel_tp.massPi0, weight)
            else: self.twoprong_masspi0_iso_sym_endcap.Fill(sel_tp.massPi0, weight)
            if bins[0] < sel_tp.pt < bins[1]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_iso_sym_barrel_bin0.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_sym_endcap_bin0.Fill(sel_tp.massPi0, weight)
            elif bins[1] < sel_tp.pt < bins[2]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_iso_sym_barrel_bin1.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_sym_endcap_bin1.Fill(sel_tp.massPi0, weight)
            elif bins[2] < sel_tp.pt < bins[3]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_iso_sym_barrel_bin2.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_sym_endcap_bin2.Fill(sel_tp.massPi0, weight)
            elif bins[3] < sel_tp.pt < bins[4]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_iso_sym_barrel_bin3.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_sym_endcap_bin3.Fill(sel_tp.massPi0, weight)
            elif bins[4] < sel_tp.pt < bins[5]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_iso_sym_barrel_bin4.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_sym_endcap_bin4.Fill(sel_tp.massPi0, weight)
            elif bins[5] < sel_tp.pt < bins[6]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_iso_sym_barrel_bin5.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_sym_endcap_bin5.Fill(sel_tp.massPi0, weight)
            elif bins[6] < sel_tp.pt < bins[7]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_iso_sym_barrel_bin6.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_sym_endcap_bin6.Fill(sel_tp.massPi0, weight)
            elif bins[7] < sel_tp.pt < bins[8]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_iso_sym_barrel_bin7.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_sym_endcap_bin7.Fill(sel_tp.massPi0, weight)
            elif bins[8] < sel_tp.pt < bins[9]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_iso_sym_barrel_bin8.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_sym_endcap_bin8.Fill(sel_tp.massPi0, weight)
            elif bins[9] < sel_tp.pt < bins[10]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_iso_sym_barrel_bin9.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_sym_endcap_bin9.Fill(sel_tp.massPi0, weight)
            else:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_iso_sym_barrel_bin10.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_sym_endcap_bin10.Fill(sel_tp.massPi0, weight)
        elif len(iso_asym_tp) != 0:  # tight loose
            self.twoprong_phi_iso_asym.Fill(sel_tp.phi, weight)
            self.twoprong_pt_iso_asym.Fill(sel_tp.pt, weight)
            self.twoprong_eta_iso_asym.Fill(sel_tp.eta, weight)
            self.twoprong_chargedIso_iso_asym.Fill(sel_tp.chargedIso, weight)
            self.ntwoprong_regions.Fill(1)
            self.twoprong_masspi0_iso_asym.Fill(sel_tp.massPi0, weight)
            if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_iso_asym_barrel.Fill(sel_tp.massPi0, weight)
            else: self.twoprong_masspi0_iso_asym_endcap.Fill(sel_tp.massPi0, weight)
            if bins[0] < sel_tp.pt < bins[1]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_iso_asym_barrel_bin0.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_asym_endcap_bin0.Fill(sel_tp.massPi0, weight)
            elif bins[1] < sel_tp.pt < bins[2]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_iso_asym_barrel_bin1.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_asym_endcap_bin1.Fill(sel_tp.massPi0, weight)
            elif bins[2] < sel_tp.pt < bins[3]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_iso_asym_barrel_bin2.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_asym_endcap_bin2.Fill(sel_tp.massPi0, weight)
            elif bins[3] < sel_tp.pt < bins[4]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_iso_asym_barrel_bin3.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_asym_endcap_bin3.Fill(sel_tp.massPi0, weight)
            elif bins[4] < sel_tp.pt < bins[5]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_iso_asym_barrel_bin4.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_asym_endcap_bin4.Fill(sel_tp.massPi0, weight)
            elif bins[5] < sel_tp.pt < bins[6]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_iso_asym_barrel_bin5.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_asym_endcap_bin5.Fill(sel_tp.massPi0, weight)
            elif bins[6] < sel_tp.pt < bins[7]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_iso_asym_barrel_bin6.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_asym_endcap_bin6.Fill(sel_tp.massPi0, weight)
            elif bins[7] < sel_tp.pt < bins[8]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_iso_asym_barrel_bin7.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_asym_endcap_bin7.Fill(sel_tp.massPi0, weight)
            elif bins[8] < sel_tp.pt < bins[9]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_iso_asym_barrel_bin8.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_asym_endcap_bin8.Fill(sel_tp.massPi0, weight)
            elif bins[9] < sel_tp.pt < bins[10]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_iso_asym_barrel_bin9.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_asym_endcap_bin9.Fill(sel_tp.massPi0, weight)
            else:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_iso_asym_barrel_bin10.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_asym_endcap_bin10.Fill(sel_tp.massPi0, weight)
        elif len(noniso_sym_tp) != 0:  # loose tight
            self.twoprong_phi_noniso_sym.Fill(sel_tp.phi, weight)
            self.twoprong_pt_noniso_sym.Fill(sel_tp.pt, weight)
            self.twoprong_eta_noniso_sym.Fill(sel_tp.eta, weight)
            self.ntwoprong_regions.Fill(2)
            self.twoprong_chargedIso_noniso_sym.Fill(sel_tp.chargedIso, weight)
            self.twoprong_masspi0_noniso_sym.Fill(sel_tp.massPi0, weight)
            if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_noniso_sym_barrel.Fill(sel_tp.massPi0, weight)
            else: self.twoprong_masspi0_noniso_sym_endcap.Fill(sel_tp.massPi0, weight)
            if bins[0] < sel_tp.pt < bins[1]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_noniso_sym_barrel_bin0.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_sym_endcap_bin0.Fill(sel_tp.massPi0, weight)
            elif bins[1] < sel_tp.pt < bins[2]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_noniso_sym_barrel_bin1.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_sym_endcap_bin1.Fill(sel_tp.massPi0, weight)
            elif bins[2] < sel_tp.pt < bins[3]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_noniso_sym_barrel_bin2.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_sym_endcap_bin2.Fill(sel_tp.massPi0, weight)
            elif bins[3] < sel_tp.pt < bins[4]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_noniso_sym_barrel_bin3.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_sym_endcap_bin3.Fill(sel_tp.massPi0, weight)
            elif bins[4] < sel_tp.pt < bins[5]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_noniso_sym_barrel_bin4.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_sym_endcap_bin4.Fill(sel_tp.massPi0, weight)
            elif bins[5] < sel_tp.pt < bins[6]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_noniso_sym_barrel_bin5.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_sym_endcap_bin5.Fill(sel_tp.massPi0, weight)
            elif bins[6] < sel_tp.pt < bins[7]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_noniso_sym_barrel_bin6.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_sym_endcap_bin6.Fill(sel_tp.massPi0, weight)
            elif bins[7] < sel_tp.pt < bins[8]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_noniso_sym_barrel_bin7.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_sym_endcap_bin7.Fill(sel_tp.massPi0, weight)
            elif bins[8] < sel_tp.pt < bins[9]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_noniso_sym_barrel_bin8.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_sym_endcap_bin8.Fill(sel_tp.massPi0, weight)
            elif bins[9] < sel_tp.pt < bins[10]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_noniso_sym_barrel_bin9.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_sym_endcap_bin9.Fill(sel_tp.massPi0, weight)
            else:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_noniso_sym_barrel_bin10.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_sym_endcap_bin10.Fill(sel_tp.massPi0, weight)
        elif len(noniso_asym_tp) != 0:  # loose loose
            self.twoprong_phi_noniso_asym.Fill(sel_tp.phi, weight)
            self.twoprong_pt_noniso_asym.Fill(sel_tp.pt, weight)
            self.twoprong_eta_noniso_asym.Fill(sel_tp.eta, weight)
            self.ntwoprong_regions.Fill(3)
            self.twoprong_chargedIso_noniso_asym.Fill(sel_tp.chargedIso, weight)
            self.twoprong_masspi0_noniso_asym.Fill(sel_tp.massPi0, weight)
            if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_noniso_asym_barrel.Fill(sel_tp.massPi0, weight)
            else: self.twoprong_masspi0_noniso_asym_endcap.Fill(sel_tp.massPi0, weight)
            if bins[0] < sel_tp.pt < bins[1]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_noniso_asym_barrel_bin0.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_asym_endcap_bin0.Fill(sel_tp.massPi0, weight)
            elif bins[1] < sel_tp.pt < bins[2]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_noniso_asym_barrel_bin1.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_asym_endcap_bin1.Fill(sel_tp.massPi0, weight)
            elif bins[2] < sel_tp.pt < bins[3]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_noniso_asym_barrel_bin2.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_asym_endcap_bin2.Fill(sel_tp.massPi0, weight)
            elif bins[3] < sel_tp.pt < bins[4]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_noniso_asym_barrel_bin3.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_asym_endcap_bin3.Fill(sel_tp.massPi0, weight)
            elif bins[4] < sel_tp.pt < bins[5]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_noniso_asym_barrel_bin4.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_asym_endcap_bin4.Fill(sel_tp.massPi0, weight)
            elif bins[5] < sel_tp.pt < bins[6]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_noniso_asym_barrel_bin5.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_asym_endcap_bin5.Fill(sel_tp.massPi0, weight)
            elif bins[6] < sel_tp.pt < bins[7]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_noniso_asym_barrel_bin6.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_asym_endcap_bin6.Fill(sel_tp.massPi0, weight)
            elif bins[7] < sel_tp.pt < bins[8]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_noniso_asym_barrel_bin7.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_asym_endcap_bin7.Fill(sel_tp.massPi0, weight)
            elif bins[8] < sel_tp.pt < bins[9]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_noniso_asym_barrel_bin8.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_asym_endcap_bin8.Fill(sel_tp.massPi0, weight)
            elif bins[9] < sel_tp.pt < bins[10]:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_noniso_asym_barrel_bin9.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_asym_endcap_bin9.Fill(sel_tp.massPi0, weight)
            else:
                if abs(sel_tp.eta) < 1.4442: self.twoprong_masspi0_noniso_asym_barrel_bin10.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_asym_endcap_bin10.Fill(sel_tp.massPi0, weight)
    
        return True


