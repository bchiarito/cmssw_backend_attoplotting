#!/usr/bin/env python
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True


def dR(eta1, eta2, phi1, phi2):
    pi = ROOT.TMath.Pi()
    dEta = abs(eta1 - eta2)
    dPhi = phi1 - phi2
    if dPhi > pi: dPhi -= 2.0*pi
    elif dPhi <= -pi: dPhi += 2.0*pi
    return ROOT.TMath.Sqrt(dEta**2 + dPhi**2)

class MyAnalysis(Module):
    def __init__(self, datamc, lumi=1.0, dict_xs=None, dict_ngen=None, phislice=0):
        self.writeHistFile = True
        self.lumi = lumi
        self.dict_xs = dict_xs
        self.dict_ngen = dict_ngen
        self.datamc = datamc
        self.phislice = phislice

    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)
        
        self.recophi_mass = ROOT.TH1F('recophi_mass', '; RecoPhi_mass', 600, 0, 6000)
        self.addObject(self.recophi_mass)

        self.twoprong_chargedIso_iso_sym = ROOT.TH1F('twoprong_chargedIso_iso_sym', '; Twoprong_chargedIso', 2000, 0, 100)
        self.twoprong_chargedIso_iso_asym = ROOT.TH1F('twoprong_chargedIso_iso_asym', '; Twoprong_chargedIso', 2000, 0, 100)
        self.twoprong_chargedIso_noniso_sym = ROOT.TH1F('twoprong_chargedIso_noniso_sym', '; Twoprong_chargedIso', 2000, 0, 100)
        self.twoprong_chargedIso_noniso_asym = ROOT.TH1F('twoprong_chargedIso_noniso_asym', '; Twoprong_chargedIso', 2000, 0, 100)
        self.addObject(self.twoprong_chargedIso_iso_sym)
        self.addObject(self.twoprong_chargedIso_iso_asym)
        self.addObject(self.twoprong_chargedIso_noniso_sym)
        self.addObject(self.twoprong_chargedIso_noniso_asym)

        self.photon_tight_hoe_iso_sym = ROOT.TH1F('photon_tight_hoe_iso_sym', '; Photon_hoe', 400, 0, 2)
        self.addObject(self.photon_tight_hoe_iso_sym)
        self.photon_loose_hoe_iso_sym = ROOT.TH1F('photon_loose_hoe_iso_sym', '; Photon_hoe', 400, 0, 2)
        self.addObject(self.photon_loose_hoe_iso_sym)
        self.photon_tight_sieie_iso_sym = ROOT.TH1F('photon_tight_sieie_iso_sym', '; Photon_sieie', 2000, 0, 2)
        self.addObject(self.photon_tight_sieie_iso_sym)
        self.photon_loose_sieie_iso_sym = ROOT.TH1F('photon_loose_sieie_iso_sym', '; Photon_sieie', 2000, 0, 2)
        self.addObject(self.photon_loose_sieie_iso_sym)
        self.photon_tight_chgIso_iso_sym = ROOT.TH1F('photon_tight_pfRelIso03_chg_iso_sym', '; Photon_pfRelIso03_chg', 100, 0, 0.01)
        self.addObject(self.photon_tight_chgIso_iso_sym)
        self.photon_loose_chgIso_iso_sym = ROOT.TH1F('photon_loose_pfRelIso03_chg_iso_sym', '; Photon_pfRelIso03_chg', 500, 0, 5)
        self.addObject(self.photon_loose_chgIso_iso_sym)
        self.photon_tight_hadTow_iso_sym = ROOT.TH1F('photon_tight_hadTow_iso_sym', '; Photon_hadTowOverEm', 400, 0, 2)
        self.addObject(self.photon_tight_hadTow_iso_sym)
        self.photon_loose_hadTow_iso_sym = ROOT.TH1F('photon_loose_hadTow_iso_sym', '; Photon_hadTowOverEm', 400, 0, 2)
        self.addObject(self.photon_loose_hadTow_iso_sym)

        self.photon_tight_hoe_iso_asym = ROOT.TH1F('photon_tight_hoe_iso_asym', '; Photon_hoe', 400, 0, 2)
        self.addObject(self.photon_tight_hoe_iso_asym)
        self.photon_loose_hoe_iso_asym = ROOT.TH1F('photon_loose_hoe_iso_asym', '; Photon_hoe', 400, 0, 2)
        self.addObject(self.photon_loose_hoe_iso_asym)
        self.photon_tight_sieie_iso_asym = ROOT.TH1F('photon_tight_sieie_iso_asym', '; Photon_sieie', 2000, 0, 2)
        self.addObject(self.photon_tight_sieie_iso_asym)
        self.photon_loose_sieie_iso_asym = ROOT.TH1F('photon_loose_sieie_iso_asym', '; Photon_sieie', 2000, 0, 2)
        self.addObject(self.photon_loose_sieie_iso_asym)
        self.photon_tight_chgIso_iso_asym = ROOT.TH1F('photon_tight_pfRelIso03_chg_iso_asym', '; Photon_pfRelIso03_chg', 100, 0, 0.01)
        self.addObject(self.photon_tight_chgIso_iso_asym)
        self.photon_loose_chgIso_iso_asym = ROOT.TH1F('photon_loose_pfRelIso03_chg_iso_asym', '; Photon_pfRelIso03_chg', 500, 0, 5)
        self.addObject(self.photon_loose_chgIso_iso_asym)
        self.photon_tight_hadTow_iso_asym = ROOT.TH1F('photon_tight_hadTow_iso_asym', '; Photon_hadTowOverEm', 400, 0, 2)
        self.addObject(self.photon_tight_hadTow_iso_asym)
        self.photon_loose_hadTow_iso_asym = ROOT.TH1F('photon_loose_hadTow_iso_asym', '; Photon_hadTowOverEm', 400, 0, 2)
        self.addObject(self.photon_loose_hadTow_iso_asym)
        
        self.photon_tight_hoe_noniso_sym = ROOT.TH1F('photon_tight_hoe_noniso_sym', '; Photon_hoe', 400, 0, 2)
        self.addObject(self.photon_tight_hoe_noniso_sym)
        self.photon_loose_hoe_noniso_sym = ROOT.TH1F('photon_loose_hoe_noniso_sym', '; Photon_hoe', 400, 0, 2)
        self.addObject(self.photon_loose_hoe_noniso_sym)
        self.photon_tight_sieie_noniso_sym = ROOT.TH1F('photon_tight_sieie_noniso_sym', '; Photon_sieie', 2000, 0, 2)
        self.addObject(self.photon_tight_sieie_noniso_sym)
        self.photon_loose_sieie_noniso_sym = ROOT.TH1F('photon_loose_sieie_noniso_sym', '; Photon_sieie', 2000, 0, 2)
        self.addObject(self.photon_loose_sieie_noniso_sym)
        self.photon_tight_chgIso_noniso_sym = ROOT.TH1F('photon_tight_pfRelIso03_chg_noniso_sym', '; Photon_pfRelIso03_chg', 100, 0, 0.01)
        self.addObject(self.photon_tight_chgIso_noniso_sym)
        self.photon_loose_chgIso_noniso_sym = ROOT.TH1F('photon_loose_pfRelIso03_chg_noniso_sym', '; Photon_pfRelIso03_chg', 500, 0, 5)
        self.addObject(self.photon_loose_chgIso_noniso_sym)
        self.photon_tight_hadTow_noniso_sym = ROOT.TH1F('photon_tight_hadTow_noniso_sym', '; Photon_hadTowOverEm', 400, 0, 2)
        self.addObject(self.photon_tight_hadTow_noniso_sym)
        self.photon_loose_hadTow_noniso_sym = ROOT.TH1F('photon_loose_hadTow_noniso_sym', '; Photon_hadTowOverEm', 400, 0, 2)
        self.addObject(self.photon_loose_hadTow_noniso_sym)

        self.photon_tight_hoe_noniso_asym = ROOT.TH1F('photon_tight_hoe_noniso_asym', '; Photon_hoe', 400, 0, 2)
        self.addObject(self.photon_tight_hoe_noniso_asym)
        self.photon_loose_hoe_noniso_asym = ROOT.TH1F('photon_loose_hoe_noniso_asym', '; Photon_hoe', 400, 0, 2)
        self.addObject(self.photon_loose_hoe_noniso_asym)
        self.photon_tight_sieie_noniso_asym = ROOT.TH1F('photon_tight_sieie_noniso_asym', '; Photon_sieie', 2000, 0, 2)
        self.addObject(self.photon_tight_sieie_noniso_asym)
        self.photon_loose_sieie_noniso_asym = ROOT.TH1F('photon_loose_sieie_noniso_asym', '; Photon_sieie', 2000, 0, 2)
        self.addObject(self.photon_loose_sieie_noniso_asym)
        self.photon_tight_chgIso_noniso_asym = ROOT.TH1F('photon_tight_pfRelIso03_chg_noniso_asym', '; Photon_pfRelIso03_chg', 100, 0, 0.01)
        self.addObject(self.photon_tight_chgIso_noniso_asym)
        self.photon_loose_chgIso_noniso_asym = ROOT.TH1F('photon_loose_pfRelIso03_chg_noniso_asym', '; Photon_pfRelIso03_chg', 500, 0, 5)
        self.addObject(self.photon_loose_chgIso_noniso_asym)
        self.photon_tight_hadTow_noniso_asym = ROOT.TH1F('photon_tight_hadTow_noniso_asym', '; Photon_hadTowOverEm', 400, 0, 2)
        self.addObject(self.photon_tight_hadTow_noniso_asym)
        self.photon_loose_hadTow_noniso_asym = ROOT.TH1F('photon_loose_hadTow_noniso_asym', '; Photon_hadTowOverEm', 400, 0, 2)
        self.addObject(self.photon_loose_hadTow_noniso_asym)
        
        self.twoprong_phi_noniso_asym = ROOT.TH1F('twoprong_phi_noniso_asym', ';Twoprong_phi', 70, -3.5, 3.5)
        self.addObject(self.twoprong_phi_noniso_asym)
        self.twoprong_phi_noniso_sym = ROOT.TH1F('twoprong_phi_noniso_sym', ';Twoprong_phi', 70, -3.5, 3.5)
        self.addObject(self.twoprong_phi_noniso_sym)
        self.twoprong_phi_iso_asym = ROOT.TH1F('twoprong_phi_iso_asym', ';Twoprong_phi', 70, -3.5, 3.5)
        self.addObject(self.twoprong_phi_iso_asym)
        self.twoprong_phi_iso_sym = ROOT.TH1F('twoprong_phi_iso_sym', ';Twoprong_phi', 70, -3.5, 3.5)
        self.addObject(self.twoprong_phi_iso_sym)

        # tight photon twoprong pT
        self.twoprong_pt_noniso_asym_tight = ROOT.TH1F('twoprong_pt_noniso_asym_tight', '; Twoprong_pt', 150, 0, 1500)  # loose-loose twoprong
        self.addObject(self.twoprong_pt_noniso_asym_tight)
        self.twoprong_pt_noniso_sym_tight = ROOT.TH1F('twoprong_pt_noniso_sym_tight', '; Twoprong_pt', 150, 0, 1500)  # loose-tight twoprong
        self.addObject(self.twoprong_pt_noniso_sym_tight)
        self.twoprong_pt_iso_asym_tight = ROOT.TH1F('twoprong_pt_iso_asym_tight', '; Twoprong_pt', 150, 0, 1500)  # tight-loose twoprong
        self.addObject(self.twoprong_pt_iso_asym_tight)
        self.twoprong_pt_iso_sym_tight = ROOT.TH1F('twoprong_pt_iso_sym_tight', '; Twoprong_pt', 150, 0, 1500)  # tight-tight twoprong
        self.addObject(self.twoprong_pt_iso_sym_tight)
        
        # loose photon twoprong pT
        self.twoprong_pt_noniso_asym_loose = ROOT.TH1F('twoprong_pt_noniso_asym_loose', '; Twoprong_pt', 150, 0, 1500)  # loose-loose twoprong
        self.addObject(self.twoprong_pt_noniso_asym_loose)
        self.twoprong_pt_noniso_sym_loose = ROOT.TH1F('twoprong_pt_noniso_sym_loose', '; Twoprong_pt', 150, 0, 1500)  # loose-loose twoprong
        self.addObject(self.twoprong_pt_noniso_sym_loose)
        self.twoprong_pt_iso_asym_loose = ROOT.TH1F('twoprong_pt_iso_asym_loose', '; Twoprong_pt', 150, 0, 1500)  # loose-loose twoprong
        self.addObject(self.twoprong_pt_iso_asym_loose)
        self.twoprong_pt_iso_sym_loose = ROOT.TH1F('twoprong_pt_iso_sym_loose', '; Twoprong_pt', 150, 0, 1500)  # loose-loose twoprong
        self.addObject(self.twoprong_pt_iso_sym_loose)
        
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

        self.twoprong_masspi0_noniso_asym_tight = ROOT.TH1D('twoprong_masspi0_noniso_asym_tight', '; Twoprong_massPi0', 1000, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_noniso_asym_tight)
        self.twoprong_masspi0_noniso_asym_loose = ROOT.TH1D('twoprong_masspi0_noniso_asym_loose', '; Twoprong_massPi0', 1000, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_noniso_asym_loose)
        self.twoprong_masspi0_iso_sym_tight = ROOT.TH1D('twoprong_masspi0_iso_sym_tight', '; Twoprong_massPi0', 1000, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_iso_sym_tight)
        self.twoprong_masspi0_iso_sym_loose = ROOT.TH1D('twoprong_masspi0_iso_sym_loose', '; Twoprong_massPi0', 1000, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_iso_sym_loose)
        self.twoprong_masspi0_noniso_sym_tight = ROOT.TH1D('twoprong_masspi0_noniso_sym_tight', '; Twoprong_massPi0', 1000, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_noniso_sym_tight)
        self.twoprong_masspi0_noniso_sym_loose = ROOT.TH1D('twoprong_masspi0_noniso_sym_loose', '; Twoprong_massPi0', 1000, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_noniso_sym_loose)
        self.twoprong_masspi0_iso_asym_tight = ROOT.TH1D('twoprong_masspi0_iso_asym_tight', '; Twoprong_massPi0', 1000, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_iso_asym_tight)
        self.twoprong_masspi0_iso_asym_loose = ROOT.TH1D('twoprong_masspi0_iso_asym_loose', '; Twoprong_massPi0', 1000, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_iso_asym_loose)
        
        self.twoprong_masspi0_noniso_asym_barrel_tight = ROOT.TH1D('twoprong_masspi0_noniso_asym_barrel_tight', '; Twoprong_massPi0', 1000, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_tight)
        self.twoprong_masspi0_noniso_asym_barrel_loose = ROOT.TH1D('twoprong_masspi0_noniso_asym_barrel_loose', '; Twoprong_massPi0', 1000, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_loose)
        self.twoprong_masspi0_iso_sym_barrel_tight = ROOT.TH1D('twoprong_masspi0_iso_sym_barrel_tight', '; Twoprong_massPi0', 1000, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_tight)
        self.twoprong_masspi0_iso_sym_barrel_loose = ROOT.TH1D('twoprong_masspi0_iso_sym_barrel_loose', '; Twoprong_massPi0', 1000, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_loose)
        self.twoprong_masspi0_noniso_sym_barrel_tight = ROOT.TH1D('twoprong_masspi0_noniso_sym_barrel_tight', '; Twoprong_massPi0', 1000, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_tight)
        self.twoprong_masspi0_noniso_sym_barrel_loose = ROOT.TH1D('twoprong_masspi0_noniso_sym_barrel_loose', '; Twoprong_massPi0', 1000, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_loose)
        self.twoprong_masspi0_iso_asym_barrel_tight = ROOT.TH1D('twoprong_masspi0_iso_asym_barrel_tight', '; Twoprong_massPi0', 1000, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_tight)
        self.twoprong_masspi0_iso_asym_barrel_loose = ROOT.TH1D('twoprong_masspi0_iso_asym_barrel_loose', '; Twoprong_massPi0', 1000, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_loose)
        
        self.twoprong_masspi0_noniso_asym_endcap_tight = ROOT.TH1D('twoprong_masspi0_noniso_asym_endcap_tight', '; Twoprong_massPi0', 1000, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_tight)
        self.twoprong_masspi0_noniso_asym_endcap_loose = ROOT.TH1D('twoprong_masspi0_noniso_asym_endcap_loose', '; Twoprong_massPi0', 1000, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_loose)
        self.twoprong_masspi0_iso_sym_endcap_tight = ROOT.TH1D('twoprong_masspi0_iso_sym_endcap_tight', '; Twoprong_massPi0', 1000, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_tight)
        self.twoprong_masspi0_iso_sym_endcap_loose = ROOT.TH1D('twoprong_masspi0_iso_sym_endcap_loose', '; Twoprong_massPi0', 1000, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_loose)
        self.twoprong_masspi0_noniso_sym_endcap_tight = ROOT.TH1D('twoprong_masspi0_noniso_sym_endcap_tight', '; Twoprong_massPi0', 1000, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_tight)
        self.twoprong_masspi0_noniso_sym_endcap_loose = ROOT.TH1D('twoprong_masspi0_noniso_sym_endcap_loose', '; Twoprong_massPi0', 1000, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_loose)
        self.twoprong_masspi0_iso_asym_endcap_tight = ROOT.TH1D('twoprong_masspi0_iso_asym_endcap_tight', '; Twoprong_massPi0', 1000, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_tight)
        self.twoprong_masspi0_iso_asym_endcap_loose = ROOT.TH1D('twoprong_masspi0_iso_asym_endcap_loose', '; Twoprong_massPi0', 1000, 0, 50)  # loose-loose twoprong
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_loose)

        self.hthat_gjets = ROOT.TH1F('GJETS_hthat_lhe', '; htHat_lhe', 150, 0, 1500)
        self.addObject(self.hthat_gjets)

        # TwoProng Pt Binning
        self.pt_bins = [20,40,60,80,100,140,180,220,300,380]
        
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

        # bin0 - bin1
        self.twoprong_masspi0_noniso_asym_barrel_bin0_tight = ROOT.TH1D('twoprong_masspi0_noniso_asym_barrel_'+str(bin0)+'_'+str(bin1)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin0_tight)
        self.twoprong_masspi0_noniso_asym_barrel_bin0_loose = ROOT.TH1D('twoprong_masspi0_noniso_asym_barrel_'+str(bin0)+'_'+str(bin1)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin0_loose)
        self.twoprong_masspi0_noniso_asym_endcap_bin0_tight = ROOT.TH1D('twoprong_masspi0_noniso_asym_endcap_'+str(bin0)+'_'+str(bin1)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin0_tight)
        self.twoprong_masspi0_noniso_asym_endcap_bin0_loose = ROOT.TH1D('twoprong_masspi0_noniso_asym_endcap_'+str(bin0)+'_'+str(bin1)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin0_loose)
        self.twoprong_masspi0_iso_sym_barrel_bin0_tight = ROOT.TH1D('twoprong_masspi0_iso_sym_barrel_'+str(bin0)+'_'+str(bin1)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin0_tight)
        self.twoprong_masspi0_iso_sym_barrel_bin0_loose = ROOT.TH1D('twoprong_masspi0_iso_sym_barrel_'+str(bin0)+'_'+str(bin1)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin0_loose)
        self.twoprong_masspi0_iso_sym_endcap_bin0_tight = ROOT.TH1D('twoprong_masspi0_iso_sym_endcap_'+str(bin0)+'_'+str(bin1)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin0_tight)
        self.twoprong_masspi0_iso_sym_endcap_bin0_loose = ROOT.TH1D('twoprong_masspi0_iso_sym_endcap_'+str(bin0)+'_'+str(bin1)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin0_loose)
        self.twoprong_masspi0_noniso_sym_barrel_bin0_tight = ROOT.TH1D('twoprong_masspi0_noniso_sym_barrel_'+str(bin0)+'_'+str(bin1)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin0_tight)
        self.twoprong_masspi0_noniso_sym_barrel_bin0_loose = ROOT.TH1D('twoprong_masspi0_noniso_sym_barrel_'+str(bin0)+'_'+str(bin1)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin0_loose)
        self.twoprong_masspi0_noniso_sym_endcap_bin0_tight = ROOT.TH1D('twoprong_masspi0_noniso_sym_endcap_'+str(bin0)+'_'+str(bin1)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin0_tight)
        self.twoprong_masspi0_noniso_sym_endcap_bin0_loose = ROOT.TH1D('twoprong_masspi0_noniso_sym_endcap_'+str(bin0)+'_'+str(bin1)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin0_loose)
        self.twoprong_masspi0_iso_asym_barrel_bin0_tight = ROOT.TH1D('twoprong_masspi0_iso_asym_barrel_'+str(bin0)+'_'+str(bin1)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin0_tight)
        self.twoprong_masspi0_iso_asym_barrel_bin0_loose = ROOT.TH1D('twoprong_masspi0_iso_asym_barrel_'+str(bin0)+'_'+str(bin1)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin0_loose)
        self.twoprong_masspi0_iso_asym_endcap_bin0_tight = ROOT.TH1D('twoprong_masspi0_iso_asym_endcap_'+str(bin0)+'_'+str(bin1)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin0_tight)
        self.twoprong_masspi0_iso_asym_endcap_bin0_loose = ROOT.TH1D('twoprong_masspi0_iso_asym_endcap_'+str(bin0)+'_'+str(bin1)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin0_loose)

        # bin1 - bin2
        self.twoprong_masspi0_noniso_asym_barrel_bin1_tight = ROOT.TH1D('twoprong_masspi0_noniso_asym_barrel_'+str(bin1)+'_'+str(bin2)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin1_tight)
        self.twoprong_masspi0_noniso_asym_barrel_bin1_loose = ROOT.TH1D('twoprong_masspi0_noniso_asym_barrel_'+str(bin1)+'_'+str(bin2)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin1_loose)
        self.twoprong_masspi0_noniso_asym_endcap_bin1_tight = ROOT.TH1D('twoprong_masspi0_noniso_asym_endcap_'+str(bin1)+'_'+str(bin2)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin1_tight)
        self.twoprong_masspi0_noniso_asym_endcap_bin1_loose = ROOT.TH1D('twoprong_masspi0_noniso_asym_endcap_'+str(bin1)+'_'+str(bin2)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin1_loose)
        self.twoprong_masspi0_iso_sym_barrel_bin1_tight = ROOT.TH1D('twoprong_masspi0_iso_sym_barrel_'+str(bin1)+'_'+str(bin2)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin1_tight)
        self.twoprong_masspi0_iso_sym_barrel_bin1_loose = ROOT.TH1D('twoprong_masspi0_iso_sym_barrel_'+str(bin1)+'_'+str(bin2)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin1_loose)
        self.twoprong_masspi0_iso_sym_endcap_bin1_tight = ROOT.TH1D('twoprong_masspi0_iso_sym_endcap_'+str(bin1)+'_'+str(bin2)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin1_tight)
        self.twoprong_masspi0_iso_sym_endcap_bin1_loose = ROOT.TH1D('twoprong_masspi0_iso_sym_endcap_'+str(bin1)+'_'+str(bin2)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin1_loose)
        self.twoprong_masspi0_noniso_sym_barrel_bin1_tight = ROOT.TH1D('twoprong_masspi0_noniso_sym_barrel_'+str(bin1)+'_'+str(bin2)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin1_tight)
        self.twoprong_masspi0_noniso_sym_barrel_bin1_loose = ROOT.TH1D('twoprong_masspi0_noniso_sym_barrel_'+str(bin1)+'_'+str(bin2)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin1_loose)
        self.twoprong_masspi0_noniso_sym_endcap_bin1_tight = ROOT.TH1D('twoprong_masspi0_noniso_sym_endcap_'+str(bin1)+'_'+str(bin2)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin1_tight)
        self.twoprong_masspi0_noniso_sym_endcap_bin1_loose = ROOT.TH1D('twoprong_masspi0_noniso_sym_endcap_'+str(bin1)+'_'+str(bin2)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin1_loose)
        self.twoprong_masspi0_iso_asym_barrel_bin1_tight = ROOT.TH1D('twoprong_masspi0_iso_asym_barrel_'+str(bin1)+'_'+str(bin2)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin1_tight)
        self.twoprong_masspi0_iso_asym_barrel_bin1_loose = ROOT.TH1D('twoprong_masspi0_iso_asym_barrel_'+str(bin1)+'_'+str(bin2)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin1_loose)
        self.twoprong_masspi0_iso_asym_endcap_bin1_tight = ROOT.TH1D('twoprong_masspi0_iso_asym_endcap_'+str(bin1)+'_'+str(bin2)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin1_tight)
        self.twoprong_masspi0_iso_asym_endcap_bin1_loose = ROOT.TH1D('twoprong_masspi0_iso_asym_endcap_'+str(bin1)+'_'+str(bin2)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin1_loose)

        # bin2 - bin3
        self.twoprong_masspi0_noniso_asym_barrel_bin2_tight = ROOT.TH1D('twoprong_masspi0_noniso_asym_barrel_'+str(bin2)+'_'+str(bin3)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin2_tight)
        self.twoprong_masspi0_noniso_asym_barrel_bin2_loose = ROOT.TH1D('twoprong_masspi0_noniso_asym_barrel_'+str(bin2)+'_'+str(bin3)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin2_loose)
        self.twoprong_masspi0_noniso_asym_endcap_bin2_tight = ROOT.TH1D('twoprong_masspi0_noniso_asym_endcap_'+str(bin2)+'_'+str(bin3)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin2_tight)
        self.twoprong_masspi0_noniso_asym_endcap_bin2_loose = ROOT.TH1D('twoprong_masspi0_noniso_asym_endcap_'+str(bin2)+'_'+str(bin3)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin2_loose)
        self.twoprong_masspi0_iso_sym_barrel_bin2_tight = ROOT.TH1D('twoprong_masspi0_iso_sym_barrel_'+str(bin2)+'_'+str(bin3)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin2_tight)
        self.twoprong_masspi0_iso_sym_barrel_bin2_loose = ROOT.TH1D('twoprong_masspi0_iso_sym_barrel_'+str(bin2)+'_'+str(bin3)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin2_loose)
        self.twoprong_masspi0_iso_sym_endcap_bin2_tight = ROOT.TH1D('twoprong_masspi0_iso_sym_endcap_'+str(bin2)+'_'+str(bin3)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin2_tight)
        self.twoprong_masspi0_iso_sym_endcap_bin2_loose = ROOT.TH1D('twoprong_masspi0_iso_sym_endcap_'+str(bin2)+'_'+str(bin3)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin2_loose)
        self.twoprong_masspi0_noniso_sym_barrel_bin2_tight = ROOT.TH1D('twoprong_masspi0_noniso_sym_barrel_'+str(bin2)+'_'+str(bin3)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin2_tight)
        self.twoprong_masspi0_noniso_sym_barrel_bin2_loose = ROOT.TH1D('twoprong_masspi0_noniso_sym_barrel_'+str(bin2)+'_'+str(bin3)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin2_loose)
        self.twoprong_masspi0_noniso_sym_endcap_bin2_tight = ROOT.TH1D('twoprong_masspi0_noniso_sym_endcap_'+str(bin2)+'_'+str(bin3)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin2_tight)
        self.twoprong_masspi0_noniso_sym_endcap_bin2_loose = ROOT.TH1D('twoprong_masspi0_noniso_sym_endcap_'+str(bin2)+'_'+str(bin3)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin2_loose)
        self.twoprong_masspi0_iso_asym_barrel_bin2_tight = ROOT.TH1D('twoprong_masspi0_iso_asym_barrel_'+str(bin2)+'_'+str(bin3)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin2_tight)
        self.twoprong_masspi0_iso_asym_barrel_bin2_loose = ROOT.TH1D('twoprong_masspi0_iso_asym_barrel_'+str(bin2)+'_'+str(bin3)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin2_loose)
        self.twoprong_masspi0_iso_asym_endcap_bin2_tight = ROOT.TH1D('twoprong_masspi0_iso_asym_endcap_'+str(bin2)+'_'+str(bin3)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin2_tight)
        self.twoprong_masspi0_iso_asym_endcap_bin2_loose = ROOT.TH1D('twoprong_masspi0_iso_asym_endcap_'+str(bin2)+'_'+str(bin3)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin2_loose)

        # bin3 - bin4
        self.twoprong_masspi0_noniso_asym_barrel_bin3_tight = ROOT.TH1D('twoprong_masspi0_noniso_asym_barrel_'+str(bin3)+'_'+str(bin4)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin3_tight)
        self.twoprong_masspi0_noniso_asym_barrel_bin3_loose = ROOT.TH1D('twoprong_masspi0_noniso_asym_barrel_'+str(bin3)+'_'+str(bin4)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin3_loose)
        self.twoprong_masspi0_noniso_asym_endcap_bin3_tight = ROOT.TH1D('twoprong_masspi0_noniso_asym_endcap_'+str(bin3)+'_'+str(bin4)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin3_tight)
        self.twoprong_masspi0_noniso_asym_endcap_bin3_loose = ROOT.TH1D('twoprong_masspi0_noniso_asym_endcap_'+str(bin3)+'_'+str(bin4)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin3_loose)
        self.twoprong_masspi0_iso_sym_barrel_bin3_tight = ROOT.TH1D('twoprong_masspi0_iso_sym_barrel_'+str(bin3)+'_'+str(bin4)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin3_tight)
        self.twoprong_masspi0_iso_sym_barrel_bin3_loose = ROOT.TH1D('twoprong_masspi0_iso_sym_barrel_'+str(bin3)+'_'+str(bin4)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin3_loose)
        self.twoprong_masspi0_iso_sym_endcap_bin3_tight = ROOT.TH1D('twoprong_masspi0_iso_sym_endcap_'+str(bin3)+'_'+str(bin4)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin3_tight)
        self.twoprong_masspi0_iso_sym_endcap_bin3_loose = ROOT.TH1D('twoprong_masspi0_iso_sym_endcap_'+str(bin3)+'_'+str(bin4)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin3_loose)
        self.twoprong_masspi0_noniso_sym_barrel_bin3_tight = ROOT.TH1D('twoprong_masspi0_noniso_sym_barrel_'+str(bin3)+'_'+str(bin4)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin3_tight)
        self.twoprong_masspi0_noniso_sym_barrel_bin3_loose = ROOT.TH1D('twoprong_masspi0_noniso_sym_barrel_'+str(bin3)+'_'+str(bin4)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin3_loose)
        self.twoprong_masspi0_noniso_sym_endcap_bin3_tight = ROOT.TH1D('twoprong_masspi0_noniso_sym_endcap_'+str(bin3)+'_'+str(bin4)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin3_tight)
        self.twoprong_masspi0_noniso_sym_endcap_bin3_loose = ROOT.TH1D('twoprong_masspi0_noniso_sym_endcap_'+str(bin3)+'_'+str(bin4)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin3_loose)
        self.twoprong_masspi0_iso_asym_barrel_bin3_tight = ROOT.TH1D('twoprong_masspi0_iso_asym_barrel_'+str(bin3)+'_'+str(bin4)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin3_tight)
        self.twoprong_masspi0_iso_asym_barrel_bin3_loose = ROOT.TH1D('twoprong_masspi0_iso_asym_barrel_'+str(bin3)+'_'+str(bin4)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin3_loose)
        self.twoprong_masspi0_iso_asym_endcap_bin3_tight = ROOT.TH1D('twoprong_masspi0_iso_asym_endcap_'+str(bin3)+'_'+str(bin4)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin3_tight)
        self.twoprong_masspi0_iso_asym_endcap_bin3_loose = ROOT.TH1D('twoprong_masspi0_iso_asym_endcap_'+str(bin3)+'_'+str(bin4)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin3_loose)
       
        # bin4 - bin5
        self.twoprong_masspi0_noniso_asym_barrel_bin4_tight = ROOT.TH1D('twoprong_masspi0_noniso_asym_barrel_'+str(bin4)+'_'+str(bin5)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin4_tight)
        self.twoprong_masspi0_noniso_asym_barrel_bin4_loose = ROOT.TH1D('twoprong_masspi0_noniso_asym_barrel_'+str(bin4)+'_'+str(bin5)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin4_loose)
        self.twoprong_masspi0_noniso_asym_endcap_bin4_tight = ROOT.TH1D('twoprong_masspi0_noniso_asym_endcap_'+str(bin4)+'_'+str(bin5)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin4_tight)
        self.twoprong_masspi0_noniso_asym_endcap_bin4_loose = ROOT.TH1D('twoprong_masspi0_noniso_asym_endcap_'+str(bin4)+'_'+str(bin5)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin4_loose)
        self.twoprong_masspi0_iso_sym_barrel_bin4_tight = ROOT.TH1D('twoprong_masspi0_iso_sym_barrel_'+str(bin4)+'_'+str(bin5)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin4_tight)
        self.twoprong_masspi0_iso_sym_barrel_bin4_loose = ROOT.TH1D('twoprong_masspi0_iso_sym_barrel_'+str(bin4)+'_'+str(bin5)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin4_loose)
        self.twoprong_masspi0_iso_sym_endcap_bin4_tight = ROOT.TH1D('twoprong_masspi0_iso_sym_endcap_'+str(bin4)+'_'+str(bin5)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin4_tight)
        self.twoprong_masspi0_iso_sym_endcap_bin4_loose = ROOT.TH1D('twoprong_masspi0_iso_sym_endcap_'+str(bin4)+'_'+str(bin5)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin4_loose)
        self.twoprong_masspi0_noniso_sym_barrel_bin4_tight = ROOT.TH1D('twoprong_masspi0_noniso_sym_barrel_'+str(bin4)+'_'+str(bin5)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin4_tight)
        self.twoprong_masspi0_noniso_sym_barrel_bin4_loose = ROOT.TH1D('twoprong_masspi0_noniso_sym_barrel_'+str(bin4)+'_'+str(bin5)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin4_loose)
        self.twoprong_masspi0_noniso_sym_endcap_bin4_tight = ROOT.TH1D('twoprong_masspi0_noniso_sym_endcap_'+str(bin4)+'_'+str(bin5)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin4_tight)
        self.twoprong_masspi0_noniso_sym_endcap_bin4_loose = ROOT.TH1D('twoprong_masspi0_noniso_sym_endcap_'+str(bin4)+'_'+str(bin5)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin4_loose)
        self.twoprong_masspi0_iso_asym_barrel_bin4_tight = ROOT.TH1D('twoprong_masspi0_iso_asym_barrel_'+str(bin4)+'_'+str(bin5)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin4_tight)
        self.twoprong_masspi0_iso_asym_barrel_bin4_loose = ROOT.TH1D('twoprong_masspi0_iso_asym_barrel_'+str(bin4)+'_'+str(bin5)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin4_loose)
        self.twoprong_masspi0_iso_asym_endcap_bin4_tight = ROOT.TH1D('twoprong_masspi0_iso_asym_endcap_'+str(bin4)+'_'+str(bin5)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin4_tight)
        self.twoprong_masspi0_iso_asym_endcap_bin4_loose = ROOT.TH1D('twoprong_masspi0_iso_asym_endcap_'+str(bin4)+'_'+str(bin5)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin4_loose)
        
        # bin5 - bin6
        self.twoprong_masspi0_noniso_asym_barrel_bin5_tight = ROOT.TH1D('twoprong_masspi0_noniso_asym_barrel_'+str(bin5)+'_'+str(bin6)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin5_tight)
        self.twoprong_masspi0_noniso_asym_barrel_bin5_loose = ROOT.TH1D('twoprong_masspi0_noniso_asym_barrel_'+str(bin5)+'_'+str(bin6)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin5_loose)
        self.twoprong_masspi0_noniso_asym_endcap_bin5_tight = ROOT.TH1D('twoprong_masspi0_noniso_asym_endcap_'+str(bin5)+'_'+str(bin6)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin5_tight)
        self.twoprong_masspi0_noniso_asym_endcap_bin5_loose = ROOT.TH1D('twoprong_masspi0_noniso_asym_endcap_'+str(bin5)+'_'+str(bin6)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin5_loose)
        self.twoprong_masspi0_iso_sym_barrel_bin5_tight = ROOT.TH1D('twoprong_masspi0_iso_sym_barrel_'+str(bin5)+'_'+str(bin6)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin5_tight)
        self.twoprong_masspi0_iso_sym_barrel_bin5_loose = ROOT.TH1D('twoprong_masspi0_iso_sym_barrel_'+str(bin5)+'_'+str(bin6)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin5_loose)
        self.twoprong_masspi0_iso_sym_endcap_bin5_tight = ROOT.TH1D('twoprong_masspi0_iso_sym_endcap_'+str(bin5)+'_'+str(bin6)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin5_tight)
        self.twoprong_masspi0_iso_sym_endcap_bin5_loose = ROOT.TH1D('twoprong_masspi0_iso_sym_endcap_'+str(bin5)+'_'+str(bin6)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin5_loose)
        self.twoprong_masspi0_noniso_sym_barrel_bin5_tight = ROOT.TH1D('twoprong_masspi0_noniso_sym_barrel_'+str(bin5)+'_'+str(bin6)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin5_tight)
        self.twoprong_masspi0_noniso_sym_barrel_bin5_loose = ROOT.TH1D('twoprong_masspi0_noniso_sym_barrel_'+str(bin5)+'_'+str(bin6)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin5_loose)
        self.twoprong_masspi0_noniso_sym_endcap_bin5_tight = ROOT.TH1D('twoprong_masspi0_noniso_sym_endcap_'+str(bin5)+'_'+str(bin6)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin5_tight)
        self.twoprong_masspi0_noniso_sym_endcap_bin5_loose = ROOT.TH1D('twoprong_masspi0_noniso_sym_endcap_'+str(bin5)+'_'+str(bin6)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin5_loose)
        self.twoprong_masspi0_iso_asym_barrel_bin5_tight = ROOT.TH1D('twoprong_masspi0_iso_asym_barrel_'+str(bin5)+'_'+str(bin6)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin5_tight)
        self.twoprong_masspi0_iso_asym_barrel_bin5_loose = ROOT.TH1D('twoprong_masspi0_iso_asym_barrel_'+str(bin5)+'_'+str(bin6)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin5_loose)
        self.twoprong_masspi0_iso_asym_endcap_bin5_tight = ROOT.TH1D('twoprong_masspi0_iso_asym_endcap_'+str(bin5)+'_'+str(bin6)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin5_tight)
        self.twoprong_masspi0_iso_asym_endcap_bin5_loose = ROOT.TH1D('twoprong_masspi0_iso_asym_endcap_'+str(bin5)+'_'+str(bin6)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin5_loose)
        
        # bin6 - bin7
        self.twoprong_masspi0_noniso_asym_barrel_bin6_tight = ROOT.TH1D('twoprong_masspi0_noniso_asym_barrel_'+str(bin6)+'_'+str(bin7)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin6_tight)
        self.twoprong_masspi0_noniso_asym_barrel_bin6_loose = ROOT.TH1D('twoprong_masspi0_noniso_asym_barrel_'+str(bin6)+'_'+str(bin7)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin6_loose)
        self.twoprong_masspi0_noniso_asym_endcap_bin6_tight = ROOT.TH1D('twoprong_masspi0_noniso_asym_endcap_'+str(bin6)+'_'+str(bin7)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin6_tight)
        self.twoprong_masspi0_noniso_asym_endcap_bin6_loose = ROOT.TH1D('twoprong_masspi0_noniso_asym_endcap_'+str(bin6)+'_'+str(bin7)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin6_loose)
        self.twoprong_masspi0_iso_sym_barrel_bin6_tight = ROOT.TH1D('twoprong_masspi0_iso_sym_barrel_'+str(bin6)+'_'+str(bin7)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin6_tight)
        self.twoprong_masspi0_iso_sym_barrel_bin6_loose = ROOT.TH1D('twoprong_masspi0_iso_sym_barrel_'+str(bin6)+'_'+str(bin7)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin6_loose)
        self.twoprong_masspi0_iso_sym_endcap_bin6_tight = ROOT.TH1D('twoprong_masspi0_iso_sym_endcap_'+str(bin6)+'_'+str(bin7)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin6_tight)
        self.twoprong_masspi0_iso_sym_endcap_bin6_loose = ROOT.TH1D('twoprong_masspi0_iso_sym_endcap_'+str(bin6)+'_'+str(bin7)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin6_loose)
        self.twoprong_masspi0_noniso_sym_barrel_bin6_tight = ROOT.TH1D('twoprong_masspi0_noniso_sym_barrel_'+str(bin6)+'_'+str(bin7)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin6_tight)
        self.twoprong_masspi0_noniso_sym_barrel_bin6_loose = ROOT.TH1D('twoprong_masspi0_noniso_sym_barrel_'+str(bin6)+'_'+str(bin7)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin6_loose)
        self.twoprong_masspi0_noniso_sym_endcap_bin6_tight = ROOT.TH1D('twoprong_masspi0_noniso_sym_endcap_'+str(bin6)+'_'+str(bin7)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin6_tight)
        self.twoprong_masspi0_noniso_sym_endcap_bin6_loose = ROOT.TH1D('twoprong_masspi0_noniso_sym_endcap_'+str(bin6)+'_'+str(bin7)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin6_loose)
        self.twoprong_masspi0_iso_asym_barrel_bin6_tight = ROOT.TH1D('twoprong_masspi0_iso_asym_barrel_'+str(bin6)+'_'+str(bin7)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin6_tight)
        self.twoprong_masspi0_iso_asym_barrel_bin6_loose = ROOT.TH1D('twoprong_masspi0_iso_asym_barrel_'+str(bin6)+'_'+str(bin7)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin6_loose)
        self.twoprong_masspi0_iso_asym_endcap_bin6_tight = ROOT.TH1D('twoprong_masspi0_iso_asym_endcap_'+str(bin6)+'_'+str(bin7)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin6_tight)
        self.twoprong_masspi0_iso_asym_endcap_bin6_loose = ROOT.TH1D('twoprong_masspi0_iso_asym_endcap_'+str(bin6)+'_'+str(bin7)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin6_loose)
        
        # bin7 - bin8
        self.twoprong_masspi0_noniso_asym_barrel_bin7_tight = ROOT.TH1D('twoprong_masspi0_noniso_asym_barrel_'+str(bin7)+'_'+str(bin8)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin7_tight)
        self.twoprong_masspi0_noniso_asym_barrel_bin7_loose = ROOT.TH1D('twoprong_masspi0_noniso_asym_barrel_'+str(bin7)+'_'+str(bin8)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin7_loose)
        self.twoprong_masspi0_noniso_asym_endcap_bin7_tight = ROOT.TH1D('twoprong_masspi0_noniso_asym_endcap_'+str(bin7)+'_'+str(bin8)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin7_tight)
        self.twoprong_masspi0_noniso_asym_endcap_bin7_loose = ROOT.TH1D('twoprong_masspi0_noniso_asym_endcap_'+str(bin7)+'_'+str(bin8)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin7_loose)
        self.twoprong_masspi0_iso_sym_barrel_bin7_tight = ROOT.TH1D('twoprong_masspi0_iso_sym_barrel_'+str(bin7)+'_'+str(bin8)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin7_tight)
        self.twoprong_masspi0_iso_sym_barrel_bin7_loose = ROOT.TH1D('twoprong_masspi0_iso_sym_barrel_'+str(bin7)+'_'+str(bin8)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin7_loose)
        self.twoprong_masspi0_iso_sym_endcap_bin7_tight = ROOT.TH1D('twoprong_masspi0_iso_sym_endcap_'+str(bin7)+'_'+str(bin8)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin7_tight)
        self.twoprong_masspi0_iso_sym_endcap_bin7_loose = ROOT.TH1D('twoprong_masspi0_iso_sym_endcap_'+str(bin7)+'_'+str(bin8)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin7_loose)
        self.twoprong_masspi0_noniso_sym_barrel_bin7_tight = ROOT.TH1D('twoprong_masspi0_noniso_sym_barrel_'+str(bin7)+'_'+str(bin8)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin7_tight)
        self.twoprong_masspi0_noniso_sym_barrel_bin7_loose = ROOT.TH1D('twoprong_masspi0_noniso_sym_barrel_'+str(bin7)+'_'+str(bin8)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin7_loose)
        self.twoprong_masspi0_noniso_sym_endcap_bin7_tight = ROOT.TH1D('twoprong_masspi0_noniso_sym_endcap_'+str(bin7)+'_'+str(bin8)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin7_tight)
        self.twoprong_masspi0_noniso_sym_endcap_bin7_loose = ROOT.TH1D('twoprong_masspi0_noniso_sym_endcap_'+str(bin7)+'_'+str(bin8)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin7_loose)
        self.twoprong_masspi0_iso_asym_barrel_bin7_tight = ROOT.TH1D('twoprong_masspi0_iso_asym_barrel_'+str(bin7)+'_'+str(bin8)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin7_tight)
        self.twoprong_masspi0_iso_asym_barrel_bin7_loose = ROOT.TH1D('twoprong_masspi0_iso_asym_barrel_'+str(bin7)+'_'+str(bin8)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin7_loose)
        self.twoprong_masspi0_iso_asym_endcap_bin7_tight = ROOT.TH1D('twoprong_masspi0_iso_asym_endcap_'+str(bin7)+'_'+str(bin8)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin7_tight)
        self.twoprong_masspi0_iso_asym_endcap_bin7_loose = ROOT.TH1D('twoprong_masspi0_iso_asym_endcap_'+str(bin7)+'_'+str(bin8)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin7_loose)
        
        # bin8 - bin9
        self.twoprong_masspi0_noniso_asym_barrel_bin8_tight = ROOT.TH1D('twoprong_masspi0_noniso_asym_barrel_'+str(bin8)+'_'+str(bin9)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin8_tight)
        self.twoprong_masspi0_noniso_asym_barrel_bin8_loose = ROOT.TH1D('twoprong_masspi0_noniso_asym_barrel_'+str(bin8)+'_'+str(bin9)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin8_loose)
        self.twoprong_masspi0_noniso_asym_endcap_bin8_tight = ROOT.TH1D('twoprong_masspi0_noniso_asym_endcap_'+str(bin8)+'_'+str(bin9)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin8_tight)
        self.twoprong_masspi0_noniso_asym_endcap_bin8_loose = ROOT.TH1D('twoprong_masspi0_noniso_asym_endcap_'+str(bin8)+'_'+str(bin9)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin8_loose)
        self.twoprong_masspi0_iso_sym_barrel_bin8_tight = ROOT.TH1D('twoprong_masspi0_iso_sym_barrel_'+str(bin8)+'_'+str(bin9)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin8_tight)
        self.twoprong_masspi0_iso_sym_barrel_bin8_loose = ROOT.TH1D('twoprong_masspi0_iso_sym_barrel_'+str(bin8)+'_'+str(bin9)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin8_loose)
        self.twoprong_masspi0_iso_sym_endcap_bin8_tight = ROOT.TH1D('twoprong_masspi0_iso_sym_endcap_'+str(bin8)+'_'+str(bin9)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin8_tight)
        self.twoprong_masspi0_iso_sym_endcap_bin8_loose = ROOT.TH1D('twoprong_masspi0_iso_sym_endcap_'+str(bin8)+'_'+str(bin9)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin8_loose)
        self.twoprong_masspi0_noniso_sym_barrel_bin8_tight = ROOT.TH1D('twoprong_masspi0_noniso_sym_barrel_'+str(bin8)+'_'+str(bin9)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin8_tight)
        self.twoprong_masspi0_noniso_sym_barrel_bin8_loose = ROOT.TH1D('twoprong_masspi0_noniso_sym_barrel_'+str(bin8)+'_'+str(bin9)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin8_loose)
        self.twoprong_masspi0_noniso_sym_endcap_bin8_tight = ROOT.TH1D('twoprong_masspi0_noniso_sym_endcap_'+str(bin8)+'_'+str(bin9)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin8_tight)
        self.twoprong_masspi0_noniso_sym_endcap_bin8_loose = ROOT.TH1D('twoprong_masspi0_noniso_sym_endcap_'+str(bin8)+'_'+str(bin9)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin8_loose)
        self.twoprong_masspi0_iso_asym_barrel_bin8_tight = ROOT.TH1D('twoprong_masspi0_iso_asym_barrel_'+str(bin8)+'_'+str(bin9)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin8_tight)
        self.twoprong_masspi0_iso_asym_barrel_bin8_loose = ROOT.TH1D('twoprong_masspi0_iso_asym_barrel_'+str(bin8)+'_'+str(bin9)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin8_loose)
        self.twoprong_masspi0_iso_asym_endcap_bin8_tight = ROOT.TH1D('twoprong_masspi0_iso_asym_endcap_'+str(bin8)+'_'+str(bin9)+'_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin8_tight)
        self.twoprong_masspi0_iso_asym_endcap_bin8_loose = ROOT.TH1D('twoprong_masspi0_iso_asym_endcap_'+str(bin8)+'_'+str(bin9)+'_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin8_loose)
        
        # bin9 +
        self.twoprong_masspi0_noniso_asym_barrel_bin9_tight = ROOT.TH1D('twoprong_masspi0_noniso_asym_barrel_'+str(bin9)+'+_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin9_tight)
        self.twoprong_masspi0_noniso_asym_barrel_bin9_loose = ROOT.TH1D('twoprong_masspi0_noniso_asym_barrel_'+str(bin9)+'+_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_barrel_bin9_loose)
        self.twoprong_masspi0_noniso_asym_endcap_bin9_tight = ROOT.TH1D('twoprong_masspi0_noniso_asym_endcap_'+str(bin9)+'+_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin9_tight)
        self.twoprong_masspi0_noniso_asym_endcap_bin9_loose = ROOT.TH1D('twoprong_masspi0_noniso_asym_endcap_'+str(bin9)+'+_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_asym_endcap_bin9_loose)
        self.twoprong_masspi0_iso_sym_barrel_bin9_tight = ROOT.TH1D('twoprong_masspi0_iso_sym_barrel_'+str(bin9)+'+_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin9_tight)
        self.twoprong_masspi0_iso_sym_barrel_bin9_loose = ROOT.TH1D('twoprong_masspi0_iso_sym_barrel_'+str(bin9)+'+_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_barrel_bin9_loose)
        self.twoprong_masspi0_iso_sym_endcap_bin9_tight = ROOT.TH1D('twoprong_masspi0_iso_sym_endcap_'+str(bin9)+'+_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin9_tight)
        self.twoprong_masspi0_iso_sym_endcap_bin9_loose = ROOT.TH1D('twoprong_masspi0_iso_sym_endcap_'+str(bin9)+'+_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_sym_endcap_bin9_loose)
        self.twoprong_masspi0_noniso_sym_barrel_bin9_tight = ROOT.TH1D('twoprong_masspi0_noniso_sym_barrel_'+str(bin9)+'+_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin9_tight)
        self.twoprong_masspi0_noniso_sym_barrel_bin9_loose = ROOT.TH1D('twoprong_masspi0_noniso_sym_barrel_'+str(bin9)+'+_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_barrel_bin9_loose)
        self.twoprong_masspi0_noniso_sym_endcap_bin9_tight = ROOT.TH1D('twoprong_masspi0_noniso_sym_endcap_'+str(bin9)+'+_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin9_tight)
        self.twoprong_masspi0_noniso_sym_endcap_bin9_loose = ROOT.TH1D('twoprong_masspi0_noniso_sym_endcap_'+str(bin9)+'+_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_noniso_sym_endcap_bin9_loose)
        self.twoprong_masspi0_iso_asym_barrel_bin9_tight = ROOT.TH1D('twoprong_masspi0_iso_asym_barrel_'+str(bin9)+'+_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin9_tight)
        self.twoprong_masspi0_iso_asym_barrel_bin9_loose = ROOT.TH1D('twoprong_masspi0_iso_asym_barrel_'+str(bin9)+'+_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_barrel_bin9_loose)
        self.twoprong_masspi0_iso_asym_endcap_bin9_tight = ROOT.TH1D('twoprong_masspi0_iso_asym_endcap_'+str(bin9)+'+_tight', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin9_tight)
        self.twoprong_masspi0_iso_asym_endcap_bin9_loose = ROOT.TH1D('twoprong_masspi0_iso_asym_endcap_'+str(bin9)+'+_loose', '; Twoprong_massPi0', 1000, 0, 50)  
        self.addObject(self.twoprong_masspi0_iso_asym_endcap_bin9_loose)
        

    def analyze(self, event):
        if self.dict_xs and self.dict_ngen:
            dataset_id = event.dataset_id
            xs = self.dict_xs[dataset_id]
            Ngen = self.dict_ngen[dataset_id]
            weight = xs * self.lumi / Ngen
        else:
            weight = 1.0
    
        twoprongs = Collection(event, "TwoProng")
        photons = Collection(event, "Photon")
        recophi = Object(event, "CBL_RecoPhi")
        PVs = Object(event, "PV")
        #pass_trigger = event.HLT_Photon35_TwoProngs35
        pass_trigger = True

        bins = []
        for i in range(len(self.pt_bins)):
            bins.append(self.pt_bins[i])

        tight_photons = []
        loose_photons = []
        for photon in photons:
            if photon.pt > 220 and photon.isScEtaEB and photon.hadTowOverEm < 0.04596: # photon cuts
                if photon.cutBased >= 1: 
                    tight_photons.append(photon)
                    break
                else: loose_photons.append(photon)
        
        if len(loose_photons) == 0 and len(tight_photons) == 0 or not pass_trigger: return True
        tight_photon = False

        if not len(tight_photons) == 0: 
            sel_photon = tight_photons[0]
            tight_photon = True
        else: sel_photon = loose_photons[0]

        if self.phislice == 1:
          phi_window = [600, 700]
          if recophi.mass < phi_window[0] or recophi.mass > phi_window[1]: return False

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

        self.recophi_mass.Fill(recophi.mass)
        if len(iso_sym_tp) != 0: # tight tight
            self.twoprong_phi_iso_sym.Fill(sel_tp.phi, weight)
            self.twoprong_eta_iso_sym.Fill(sel_tp.eta, weight)
            self.twoprong_chargedIso_iso_sym.Fill(sel_tp.chargedIso, weight)
            if not tight_photon:
                self.twoprong_pt_iso_sym_loose.Fill(sel_tp.pt, weight)
                self.twoprong_masspi0_iso_sym_loose.Fill(sel_tp.massPi0, weight)
                self.photon_loose_chgIso_iso_sym.Fill(sel_photon.pfRelIso03_chg)
                self.photon_loose_hoe_iso_sym.Fill(sel_photon.hoe, weight)
                self.photon_loose_sieie_iso_sym.Fill(sel_photon.sieie, weight)
                self.photon_loose_hadTow_iso_sym.Fill(sel_photon.hadTowOverEm, weight)
            else:
                self.twoprong_pt_iso_sym_tight.Fill(sel_tp.pt, weight)
                self.twoprong_masspi0_iso_sym_tight.Fill(sel_tp.massPi0, weight)
                self.photon_tight_chgIso_iso_sym.Fill(sel_photon.pfRelIso03_chg)
                self.photon_tight_hoe_iso_sym.Fill(sel_photon.hoe, weight)
                self.photon_tight_sieie_iso_sym.Fill(sel_photon.sieie, weight)
                self.photon_tight_hadTow_iso_sym.Fill(sel_photon.hadTowOverEm, weight)
            if abs(sel_tp.eta) < 1.4442:
                if tight_photon: self.twoprong_masspi0_iso_sym_barrel_tight.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_sym_barrel_loose.Fill(sel_tp.massPi0, weight)
            else: 
                if tight_photon: self.twoprong_masspi0_iso_sym_endcap_tight.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_sym_endcap_loose.Fill(sel_tp.massPi0, weight)
            if bins[0] < sel_tp.pt < bins[1]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_iso_sym_barrel_bin0_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_sym_barrel_bin0_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_iso_sym_endcap_bin0_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_sym_endcap_bin0_loose.Fill(sel_tp.massPi0, weight)
            elif bins[1] < sel_tp.pt < bins[2]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_iso_sym_barrel_bin1_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_sym_barrel_bin1_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_iso_sym_endcap_bin1_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_sym_endcap_bin1_loose.Fill(sel_tp.massPi0, weight)
            elif bins[2] < sel_tp.pt < bins[3]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_iso_sym_barrel_bin2_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_sym_barrel_bin2_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_iso_sym_endcap_bin2_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_sym_endcap_bin2_loose.Fill(sel_tp.massPi0, weight)
            elif bins[3] < sel_tp.pt < bins[4]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_iso_sym_barrel_bin3_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_sym_barrel_bin3_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_iso_sym_endcap_bin3_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_sym_endcap_bin3_loose.Fill(sel_tp.massPi0, weight)
            elif bins[4] < sel_tp.pt < bins[5]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_iso_sym_barrel_bin4_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_sym_barrel_bin4_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_iso_sym_endcap_bin4_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_sym_endcap_bin4_loose.Fill(sel_tp.massPi0, weight)
            elif bins[5] < sel_tp.pt < bins[6]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_iso_sym_barrel_bin5_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_sym_barrel_bin5_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_iso_sym_endcap_bin5_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_sym_endcap_bin5_loose.Fill(sel_tp.massPi0, weight)
            elif bins[6] < sel_tp.pt < bins[7]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_iso_sym_barrel_bin6_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_sym_barrel_bin6_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_iso_sym_endcap_bin6_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_sym_endcap_bin6_loose.Fill(sel_tp.massPi0, weight)
            elif bins[7] < sel_tp.pt < bins[8]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_iso_sym_barrel_bin7_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_sym_barrel_bin7_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_iso_sym_endcap_bin7_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_sym_endcap_bin7_loose.Fill(sel_tp.massPi0, weight)
            elif bins[8] < sel_tp.pt < bins[9]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_iso_sym_barrel_bin8_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_sym_barrel_bin8_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_iso_sym_endcap_bin8_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_sym_endcap_bin8_loose.Fill(sel_tp.massPi0, weight)
            else:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_iso_sym_barrel_bin9_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_sym_barrel_bin9_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_iso_sym_endcap_bin9_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_sym_endcap_bin9_loose.Fill(sel_tp.massPi0, weight)
        elif len(iso_asym_tp) != 0:  # tight loose
            self.twoprong_phi_iso_asym.Fill(sel_tp.phi, weight)
            self.twoprong_eta_iso_asym.Fill(sel_tp.eta, weight)
            self.twoprong_chargedIso_iso_asym.Fill(sel_tp.chargedIso, weight)
            if not tight_photon:
                self.twoprong_pt_iso_asym_loose.Fill(sel_tp.pt, weight)
                self.twoprong_masspi0_iso_asym_loose.Fill(sel_tp.massPi0, weight)
                self.photon_loose_chgIso_iso_asym.Fill(sel_photon.pfRelIso03_chg)
                self.photon_loose_hoe_iso_asym.Fill(sel_photon.hoe, weight)
                self.photon_loose_sieie_iso_asym.Fill(sel_photon.sieie, weight)
                self.photon_loose_hadTow_iso_asym.Fill(sel_photon.hadTowOverEm, weight)
            else:
                self.twoprong_pt_iso_asym_tight.Fill(sel_tp.pt, weight)
                self.twoprong_masspi0_iso_asym_tight.Fill(sel_tp.massPi0, weight)
                self.photon_tight_chgIso_iso_asym.Fill(sel_photon.pfRelIso03_chg)
                self.photon_tight_hoe_iso_asym.Fill(sel_photon.hoe, weight)
                self.photon_tight_sieie_iso_asym.Fill(sel_photon.sieie, weight)
                self.photon_tight_hadTow_iso_asym.Fill(sel_photon.hadTowOverEm, weight)
            if abs(sel_tp.eta) < 1.4442:
                if tight_photon: self.twoprong_masspi0_iso_asym_barrel_tight.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_asym_barrel_loose.Fill(sel_tp.massPi0, weight)
            else: 
                if tight_photon: self.twoprong_masspi0_iso_asym_endcap_tight.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_iso_asym_endcap_loose.Fill(sel_tp.massPi0, weight)
            if bins[0] < sel_tp.pt < bins[1]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_iso_asym_barrel_bin0_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_asym_barrel_bin0_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_iso_asym_endcap_bin0_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_asym_endcap_bin0_loose.Fill(sel_tp.massPi0, weight)
            elif bins[1] < sel_tp.pt < bins[2]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_iso_asym_barrel_bin1_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_asym_barrel_bin1_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_iso_asym_endcap_bin1_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_asym_endcap_bin1_loose.Fill(sel_tp.massPi0, weight)
            elif bins[2] < sel_tp.pt < bins[3]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_iso_asym_barrel_bin2_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_asym_barrel_bin2_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_iso_asym_endcap_bin2_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_asym_endcap_bin2_loose.Fill(sel_tp.massPi0, weight)
            elif bins[3] < sel_tp.pt < bins[4]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_iso_asym_barrel_bin3_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_asym_barrel_bin3_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_iso_asym_endcap_bin3_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_asym_endcap_bin3_loose.Fill(sel_tp.massPi0, weight)
            elif bins[4] < sel_tp.pt < bins[5]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_iso_asym_barrel_bin4_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_asym_barrel_bin4_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_iso_asym_endcap_bin4_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_asym_endcap_bin4_loose.Fill(sel_tp.massPi0, weight)
            elif bins[5] < sel_tp.pt < bins[6]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_iso_asym_barrel_bin5_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_asym_barrel_bin5_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_iso_asym_endcap_bin5_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_asym_endcap_bin5_loose.Fill(sel_tp.massPi0, weight)
            elif bins[6] < sel_tp.pt < bins[7]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_iso_asym_barrel_bin6_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_asym_barrel_bin6_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_iso_asym_endcap_bin6_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_asym_endcap_bin6_loose.Fill(sel_tp.massPi0, weight)
            elif bins[7] < sel_tp.pt < bins[8]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_iso_asym_barrel_bin7_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_asym_barrel_bin7_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_iso_asym_endcap_bin7_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_asym_endcap_bin7_loose.Fill(sel_tp.massPi0, weight)
            elif bins[8] < sel_tp.pt < bins[9]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_iso_asym_barrel_bin8_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_asym_barrel_bin8_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_iso_asym_endcap_bin8_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_asym_endcap_bin8_loose.Fill(sel_tp.massPi0, weight)
            else:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_iso_asym_barrel_bin9_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_asym_barrel_bin9_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_iso_asym_endcap_bin9_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_iso_asym_endcap_bin9_loose.Fill(sel_tp.massPi0, weight)
        elif len(noniso_sym_tp) != 0:  # loose tight
            self.twoprong_phi_noniso_sym.Fill(sel_tp.phi, weight)
            self.twoprong_eta_noniso_sym.Fill(sel_tp.eta, weight)
            self.twoprong_chargedIso_noniso_sym.Fill(sel_tp.chargedIso, weight)
            if not tight_photon:
                self.twoprong_pt_noniso_sym_loose.Fill(sel_tp.pt, weight)
                self.twoprong_masspi0_noniso_sym_loose.Fill(sel_tp.massPi0, weight)
                self.photon_loose_chgIso_noniso_sym.Fill(sel_photon.pfRelIso03_chg)
                self.photon_loose_hoe_noniso_sym.Fill(sel_photon.hoe, weight)
                self.photon_loose_sieie_noniso_sym.Fill(sel_photon.sieie, weight)
                self.photon_loose_hadTow_noniso_sym.Fill(sel_photon.hadTowOverEm, weight)
            else:
                self.twoprong_pt_noniso_sym_tight.Fill(sel_tp.pt, weight)
                self.twoprong_masspi0_noniso_sym_tight.Fill(sel_tp.massPi0, weight)
                self.photon_tight_chgIso_noniso_sym.Fill(sel_photon.pfRelIso03_chg)
                self.photon_tight_hoe_noniso_sym.Fill(sel_photon.hoe, weight)
                self.photon_tight_sieie_noniso_sym.Fill(sel_photon.sieie, weight)
                self.photon_tight_hadTow_noniso_sym.Fill(sel_photon.hadTowOverEm, weight)
            if abs(sel_tp.eta) < 1.4442:
                if tight_photon: self.twoprong_masspi0_noniso_sym_barrel_tight.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_sym_barrel_loose.Fill(sel_tp.massPi0, weight)
            else: 
                if tight_photon: self.twoprong_masspi0_noniso_sym_endcap_tight.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_sym_endcap_loose.Fill(sel_tp.massPi0, weight)
            if bins[0] < sel_tp.pt < bins[1]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_noniso_sym_barrel_bin0_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_sym_barrel_bin0_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_noniso_sym_endcap_bin0_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_sym_endcap_bin0_loose.Fill(sel_tp.massPi0, weight)
            elif bins[1] < sel_tp.pt < bins[2]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_noniso_sym_barrel_bin1_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_sym_barrel_bin1_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_noniso_sym_endcap_bin1_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_sym_endcap_bin1_loose.Fill(sel_tp.massPi0, weight)
            elif bins[2] < sel_tp.pt < bins[3]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_noniso_sym_barrel_bin2_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_sym_barrel_bin2_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_noniso_sym_endcap_bin2_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_sym_endcap_bin2_loose.Fill(sel_tp.massPi0, weight)
            elif bins[3] < sel_tp.pt < bins[4]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_noniso_sym_barrel_bin3_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_sym_barrel_bin3_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_noniso_sym_endcap_bin3_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_sym_endcap_bin3_loose.Fill(sel_tp.massPi0, weight)
            elif bins[4] < sel_tp.pt < bins[5]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_noniso_sym_barrel_bin4_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_sym_barrel_bin4_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_noniso_sym_endcap_bin4_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_sym_endcap_bin4_loose.Fill(sel_tp.massPi0, weight)
            elif bins[5] < sel_tp.pt < bins[6]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_noniso_sym_barrel_bin5_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_sym_barrel_bin5_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_noniso_sym_endcap_bin5_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_sym_endcap_bin5_loose.Fill(sel_tp.massPi0, weight)
            elif bins[6] < sel_tp.pt < bins[7]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_noniso_sym_barrel_bin6_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_sym_barrel_bin6_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_noniso_sym_endcap_bin6_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_sym_endcap_bin6_loose.Fill(sel_tp.massPi0, weight)
            elif bins[7] < sel_tp.pt < bins[8]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_noniso_sym_barrel_bin7_tight.Fill(sel_tp.massPi0, weight)
                    else: 
                        self.twoprong_masspi0_noniso_sym_barrel_bin7_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_noniso_sym_endcap_bin7_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_sym_endcap_bin7_loose.Fill(sel_tp.massPi0, weight)
            elif bins[8] < sel_tp.pt < bins[9]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_noniso_sym_barrel_bin8_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_sym_barrel_bin8_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_noniso_sym_endcap_bin8_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_sym_endcap_bin8_loose.Fill(sel_tp.massPi0, weight)
            else:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_noniso_sym_barrel_bin9_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_sym_barrel_bin9_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_noniso_sym_endcap_bin9_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_sym_endcap_bin9_loose.Fill(sel_tp.massPi0, weight)
        elif len(noniso_asym_tp) != 0:  # loose loose
            self.twoprong_phi_noniso_asym.Fill(sel_tp.phi, weight)
            self.twoprong_eta_noniso_asym.Fill(sel_tp.eta, weight)
            self.twoprong_chargedIso_noniso_asym.Fill(sel_tp.chargedIso, weight)
            if not tight_photon:
                self.twoprong_pt_noniso_asym_loose.Fill(sel_tp.pt, weight)
                self.twoprong_masspi0_noniso_asym_loose.Fill(sel_tp.massPi0, weight)
                self.photon_loose_chgIso_noniso_asym.Fill(sel_photon.pfRelIso03_chg)
                self.photon_loose_hoe_noniso_asym.Fill(sel_photon.hoe, weight)
                self.photon_loose_sieie_noniso_asym.Fill(sel_photon.sieie, weight)
                self.photon_loose_hadTow_noniso_asym.Fill(sel_photon.hadTowOverEm, weight)
            else:
                self.twoprong_pt_noniso_asym_tight.Fill(sel_tp.pt, weight)
                self.twoprong_masspi0_noniso_asym_tight.Fill(sel_tp.massPi0, weight)
                self.photon_tight_chgIso_noniso_asym.Fill(sel_photon.pfRelIso03_chg)
                self.photon_tight_hoe_noniso_asym.Fill(sel_photon.hoe, weight)
                self.photon_tight_sieie_noniso_asym.Fill(sel_photon.sieie, weight)
                self.photon_tight_hadTow_noniso_asym.Fill(sel_photon.hadTowOverEm, weight)
            if abs(sel_tp.eta) < 1.4442:
                if tight_photon: self.twoprong_masspi0_noniso_asym_barrel_tight.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_asym_barrel_loose.Fill(sel_tp.massPi0, weight)
            else: 
                if tight_photon: self.twoprong_masspi0_noniso_asym_endcap_tight.Fill(sel_tp.massPi0, weight)
                else: self.twoprong_masspi0_noniso_asym_endcap_loose.Fill(sel_tp.massPi0, weight)
            if bins[0] < sel_tp.pt < bins[1]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_noniso_asym_barrel_bin0_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_asym_barrel_bin0_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_noniso_asym_endcap_bin0_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_asym_endcap_bin0_loose.Fill(sel_tp.massPi0, weight)
            elif bins[1] < sel_tp.pt < bins[2]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_noniso_asym_barrel_bin1_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_asym_barrel_bin1_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_noniso_asym_endcap_bin1_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_asym_endcap_bin1_loose.Fill(sel_tp.massPi0, weight)
            elif bins[2] < sel_tp.pt < bins[3]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_noniso_asym_barrel_bin2_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_asym_barrel_bin2_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_noniso_asym_endcap_bin2_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_asym_endcap_bin2_loose.Fill(sel_tp.massPi0, weight)
            elif bins[3] < sel_tp.pt < bins[4]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_noniso_asym_barrel_bin3_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_asym_barrel_bin3_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_noniso_asym_endcap_bin3_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_asym_endcap_bin3_loose.Fill(sel_tp.massPi0, weight)
            elif bins[4] < sel_tp.pt < bins[5]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_noniso_asym_barrel_bin4_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_asym_barrel_bin4_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_noniso_asym_endcap_bin4_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_asym_endcap_bin4_loose.Fill(sel_tp.massPi0, weight)
            elif bins[5] < sel_tp.pt < bins[6]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_noniso_asym_barrel_bin5_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_asym_barrel_bin5_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_noniso_asym_endcap_bin5_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_asym_endcap_bin5_loose.Fill(sel_tp.massPi0, weight)
            elif bins[6] < sel_tp.pt < bins[7]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_noniso_asym_barrel_bin6_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_asym_barrel_bin6_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_noniso_asym_endcap_bin6_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_asym_endcap_bin6_loose.Fill(sel_tp.massPi0, weight)
            elif bins[7] < sel_tp.pt < bins[8]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_noniso_asym_barrel_bin7_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_asym_barrel_bin7_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_noniso_asym_endcap_bin7_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_asym_endcap_bin7_loose.Fill(sel_tp.massPi0, weight)
            elif bins[8] < sel_tp.pt < bins[9]:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_noniso_asym_barrel_bin8_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_asym_barrel_bin8_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_noniso_asym_endcap_bin8_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_asym_endcap_bin8_loose.Fill(sel_tp.massPi0, weight)
            else:
                if abs(sel_tp.eta) < 1.4442: 
                    if tight_photon: self.twoprong_masspi0_noniso_asym_barrel_bin9_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_asym_barrel_bin9_loose.Fill(sel_tp.massPi0, weight)
                else:
                    if tight_photon: self.twoprong_masspi0_noniso_asym_endcap_bin9_tight.Fill(sel_tp.massPi0, weight)
                    else: self.twoprong_masspi0_noniso_asym_endcap_bin9_loose.Fill(sel_tp.massPi0, weight)
        
        return True


