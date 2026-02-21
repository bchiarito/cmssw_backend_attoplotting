#!/usr/bin/env python
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
from array import array
import os
import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

PHOTON_CUTBASED_ID = 1 # loose
PHOTON_BARREL_ETA = 1.4442
PHOTON_HoverE_CUT = 0.04596

PHI_BINS = [500, 520, 541, 563, 586, 609, 633, 658, 684, 711, 739, 769, 800, 832, 865, 900, 936, 973, 1012, 1052, 1094, 1138, 1184, 1231, 1280, 1331, 1384, 1439, 1497, 1557, 1619, 1684, 1751, 1821, 1894, 1970, 2049, 2131, 2216, 2305, 2397, 2493, 2593, 2697, 2805, 2917, 3034, 3155, 3281, 3412, 3548, 3690, 3838, 3998]
OMEGA_BINS = [0.40, 0.53, 0.58, 0.64, 0.70, 0.77, 0.85, 0.94, 1.03, 1.13, 1.24, 1.36, 1.50, 1.65, 1.81, 1.99, 2.19, 2.41, 2.65, 2.92, 3.21, 3.5, 3.85, 4.24, 4.66, 5.33]

M2P_GEN = 3.0

def dR(eta1, eta2, phi1, phi2):
    pi = ROOT.TMath.Pi()
    dEta = abs(eta1 - eta2)
    dPhi = phi1 - phi2
    if dPhi > pi: dPhi -= 2.0*pi
    elif dPhi <= -pi: dPhi += 2.0*pi
    return ROOT.TMath.Sqrt(dEta**2 + dPhi**2)

def dPhi(vec1, vec2):
  pi = ROOT.TMath.Pi()
  dphi = vec1.Phi() - vec2.Phi()
  if dphi > pi: dphi = dphi - 2.0*pi
  elif dphi <= -pi: dphi = dphi + 2*pi
  return dphi

def get_vec(obj):
  return ROOT.Math.PtEtaPhiMVector(obj.pt, obj.eta, obj.phi, obj.mass)

def clone_vec(vec):
  return ROOT.Math.PtEtaPhiMVector(vec.Pt(), vec.Eta(), vec.Phi(), vec.M())

class SanityAnalysis(Module):
    def __init__(self, datamc, lumi=1.0, dict_xs=None, dict_ngen=None, cut='None', photon='HPID'):
        self.writeHistFile = True
        self.lumi = lumi
        self.dict_xs = dict_xs
        self.dict_ngen = dict_ngen
        self.cut = cut

        self.seed = 1337
        self.rand = ROOT.TRandom(self.seed)

        self.datamc = datamc
        self.year = "18" # default
        if "data" in self.datamc:
            if not len(self.datamc) == 4: 
                self.year = self.datamc[4:]
                self.datamc = self.datamc[:4]
        elif "mc" in self.datamc:
            if not len(self.datamc) == 2:
                self.year = self.datamc[2:]
                self.datamc = self.datamc[:2]
        elif "sigRes" in self.datamc:
            if not len(self.datamc) == 6:
                self.year = self.datamc[6:]
                self.datamc = self.datamc[:6]
        elif "sigNonRes" in self.datamc:
            if not len(self.datamc) == 9:
                self.year = self.datamc[9:]
                self.datamc = self.datamc[:9]

        self.photon = photon
        self.photon_min_pt = 220  # default
        if "CBL" in self.photon: 
            if not len(self.photon) == 3: 
                self.photon_min_pt = int(self.photon[3:])
                self.photon = self.photon[:3]
        elif "HPID" in self.photon: 
            if not len(self.photon) == 4: 
                self.photon_min_pt = int(self.photon[4:])
                self.photon = self.photon[:4]
        print('WWW')
        print(self.photon_min_pt)
        print('WWW')

        '''
        if self.dict_xs and self.dict_ngen:
          dataset_id = event.dataset_id
          xs = self.dict_xs[dataset_id]
          Ngen = self.dict_ngen[dataset_id]
          self.weight = xs * self.lumi / Ngen
        else:
          self.weight = 1.0
        '''
        self.weight = 1.0
        print("weight", self.weight)

    def book_histo(self, name, bins, low, high, title=None):
        if not title: title = name
        exec("self."+name+" = ROOT.TH1F(name, title, bins, low, high)")
        exec("self.addObject(self."+name+")")

    def book_manyEBEB(self, name, var, bins, low, high):
        self.book_histo(name+'_'+var, bins, low, high)
        self.book_histo(name+'_B_B_'+var, bins, low, high)
        self.book_histo(name+'_B_E_'+var, bins, low, high)
        #self.book_histo(name+'_E_B_'+var, bins, low, high)
        #self.book_histo(name+'_E_E_'+var, bins, low, high)

    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)

        self.njets = ROOT.TH1F('njets', 'njets', 20, 0, 20)
        self.addObject(self.njets)
        self.ht = ROOT.TH1F('ht', 'ht', 100, 0, 5000)
        self.addObject(self.ht)
        self.met = ROOT.TH1F('met', 'met', 80, 0, 400)
        self.addObject(self.met)
        self.met_phi = ROOT.TH1F('met_phi', 'met_phi', 70, -3.5, 3.5)
        self.addObject(self.met_phi)
        self.met_phi_cut = ROOT.TH1F('met_phi_cut', 'met_phi_cut', 70, -3.5, 3.5)
        self.addObject(self.met_phi_cut)
        self.npv = ROOT.TH1F('npv', 'npv', 100, 0, 100)
        self.addObject(self.npv)

        self.book_manyEBEB('recophi', 'dphi', 35, 0, 3.5)
        self.book_manyEBEB('recophi', 'deta', 70, 0, 7)
        self.book_manyEBEB('recophi', 'dr', 70, 0, 7)
        self.book_manyEBEB('recophi', 'm', 200, 0, 6000)
        self.book_manyEBEB('recophi', 'pt', 200, 0, 2000)
        self.book_manyEBEB('recophi', 'eta', 200, -10, 10)
        self.book_manyEBEB('recophi', 'phi', 70, -3.5, 3.5)
        self.book_manyEBEB('recophi', 'dr_0jet', 70, 0, 7)
        self.book_manyEBEB('recophi', 'dr_1jet', 70, 0, 7)
        self.book_manyEBEB('recophi', 'dr_2jet', 70, 0, 7)

        self.photon_pt = ROOT.TH1F('photon_pt', 'photon_pt', 160, 0, 1600)
        self.addObject(self.photon_pt)
        self.photon_pt_up = ROOT.TH1F('photon_pt_up', 'photon_pt_up', 160, 0, 1600)
        self.addObject(self.photon_pt_up)
        self.photon_pt_down = ROOT.TH1F('photon_pt_down', 'photon_pt_down', 160, 0, 1600)
        self.addObject(self.photon_pt_down)
        self.photon_eta = ROOT.TH1F('photon_eta', 'photon_eta', 100, -5, 5)
        self.addObject(self.photon_eta)
        self.photon_phi = ROOT.TH1F('photon_phi', 'photon_phi', 70, -3.5, 3.5)
        self.addObject(self.photon_phi)
        self.nphoton = ROOT.TH1F('nphoton', 'nphoton', 10, 0, 10)
        self.addObject(self.nphoton)

        self.twoprong_pt = ROOT.TH1F('twoprong_pt', 'twoprong_pt', 150, 0, 1500)
        self.addObject(self.twoprong_pt)
        self.twoprong_eta = ROOT.TH1F('twoprong_eta', 'twoprong_eta', 100, -5, 5)
        self.addObject(self.twoprong_eta)

        self.twoprong_phi = ROOT.TH1F('twoprong_phi', 'twoprong_phi', 70, -3.5, 3.5)
        self.addObject(self.twoprong_phi)
        self.twoprong_mass = ROOT.TH1F('twoprong_mass', 'twoprong_mass', 300, 0, 15)
        self.addObject(self.twoprong_mass)
        self.twoprong_masspi0 = ROOT.TH1F('twoprong_masspi0', 'twoprong_masspi0', 300, 0, 15)
        self.addObject(self.twoprong_masspi0)
        self.twoprong_masseta = ROOT.TH1F('twoprong_masseta', 'twoprong_masseta', 300, 0, 15)
        self.addObject(self.twoprong_masseta)
        self.ntwoprong = ROOT.TH1F('ntwoprong', 'ntwoprong', 10, 0, 10)
        self.addObject(self.ntwoprong)
        self.twoprong_eta_barrel = ROOT.TH1F('twoprong_eta_barrel', 'twoprong_eta_barrel', 100, -5, 5)
        self.addObject(self.twoprong_eta_barrel)
        self.twoprong_eta_endcap = ROOT.TH1F('twoprong_eta_endcap', 'twoprong_eta_endcap', 100, -5, 5)
        self.addObject(self.twoprong_eta_endcap)

        self.recomass_bare = ROOT.TH2D('recomass_bare', 'recomass_bare', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_bare)
        self.recomass_bare_barrel = ROOT.TH2D('recomass_bare_barrel', 'recomass_barrel', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_bare_barrel)
        self.recomass_bare_endcap = ROOT.TH2D('recomass_bare_endcap', 'recomass_endcap', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_bare_endcap)

        self.recomass_sideband = ROOT.TH2D('recomass_sideband', 'recomass_sideband', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_sideband)
        self.recomass_sideband_barrel = ROOT.TH2D('recomass_sideband_barrel', 'recomass_sideband_barrel', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_sideband_barrel)
        self.recomass_sideband_endcap = ROOT.TH2D('recomass_sideband_endcap', 'recomass_sideband_endcap', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_sideband_endcap)
        self.recomass_uniformbinning = ROOT.TH2D('recomass_uniformbinning', 'recomass_uniformbinning', 300, 0, 15, 200, 0, 6000)
        self.addObject(self.recomass_uniformbinning)


        self.recomassprime_bare = ROOT.TH2D('recomassprime_bare', 'recomassprime_bare', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_bare)
        self.recomassprime_bare_barrel = ROOT.TH2D('recomassprime_bare_barrel', 'recomassprime_barrel', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_bare_barrel)
        self.recomassprime_bare_endcap = ROOT.TH2D('recomassprime_bare_endcap', 'recomassprime_endcap', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_bare_endcap)

        self.recomassprime_sideband = ROOT.TH2D('recomassprime_sideband', 'recomassprime_sideband', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_sideband)
        self.recomassprime_sideband_barrel = ROOT.TH2D('recomassprime_sideband_barrel', 'recomassprime_sideband_barrel', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_sideband_barrel)
        self.recomassprime_sideband_endcap = ROOT.TH2D('recomassprime_sideband_endcap', 'recomassprime_sideband_endcap', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_sideband_endcap)
        self.recomassprime_uniformbinning = ROOT.TH2D('recomassprime_uniformbinning', 'recomassprime_uniformbinning', 300, 0, 15, 200, 0, 6000)
        self.addObject(self.recomassprime_uniformbinning)


        self.recomass = ROOT.TH2D('recomass', 'recomass', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass)

        self.recomass_shiftUp = ROOT.TH2D('recomass_shiftUp', 'recomass_shiftUp', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_shiftUp)
        self.recomass_shiftDown = ROOT.TH2D('recomass_shiftDown', 'recomass_shiftDown', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_shiftDown)

        self.recomass_stretchUp = ROOT.TH2D('recomass_stretchUp', 'recomass_stretchUp', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_stretchUp)
        self.recomass_stretchDown = ROOT.TH2D('recomass_stretchDown', 'recomass_stretchDown', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_stretchDown)

        self.recomass_scaleUp = ROOT.TH2D('recomass_scaleUp', 'recomass_scaleUp', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_scaleUp)
        self.recomass_scaleDown = ROOT.TH2D('recomass_scaleDown', 'recomass_scaleDown', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_scaleDown)

        self.recomass_resUp = ROOT.TH2D('recomass_resUp', 'recomass_resUp', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_resUp)
        self.recomass_resDown = ROOT.TH2D('recomass_resDown', 'recomass_resDown', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_resDown)

        self.recomass_barrel = ROOT.TH2D('recomass_barrel', 'recomass_barrel', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_barrel)

        self.recomass_barrel_shiftUp = ROOT.TH2D('recomass_barrel_shiftUp', 'recomass_barrel_shiftUp', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_barrel_shiftUp)
        self.recomass_barrel_shiftDown = ROOT.TH2D('recomass_barrel_shiftDown', 'recomass_barrel_shiftDown', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_barrel_shiftDown)

        self.recomass_barrel_stretchUp = ROOT.TH2D('recomass_barrel_stretchUp', 'recomass_barrel_stretchUp', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_barrel_stretchUp)
        self.recomass_barrel_stretchDown = ROOT.TH2D('recomass_barrel_stretchDown', 'recomass_barrel_stretchDown', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_barrel_stretchDown)

        self.recomass_barrel_scaleUp = ROOT.TH2D('recomass_barrel_scaleUp', 'recomass_barrel_scaleUp', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_barrel_scaleUp)
        self.recomass_barrel_scaleDown = ROOT.TH2D('recomass_barrel_scaleDown', 'recomass_barrel_scaleDown', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_barrel_scaleDown)

        self.recomass_barrel_resUp = ROOT.TH2D('recomass_barrel_resUp', 'recomass_barrel_resUp', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_barrel_resUp)
        self.recomass_barrel_resDown = ROOT.TH2D('recomass_barrel_resDown', 'recomass_barrel_resDown', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_barrel_resDown)

        self.recomass_endcap = ROOT.TH2D('recomass_endcap', 'recomass_endcap', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_endcap)

        self.recomass_endcap_shiftUp = ROOT.TH2D('recomass_endcap_shiftUp', 'recomass_endcap_shiftUp', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_endcap_shiftUp)
        self.recomass_endcap_shiftDown = ROOT.TH2D('recomass_endcap_shiftDown', 'recomass_endcap_shiftDown', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_endcap_shiftDown)

        self.recomass_endcap_stretchUp = ROOT.TH2D('recomass_endcap_stretchUp', 'recomass_endcap_stretchUp', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_endcap_stretchUp)
        self.recomass_endcap_stretchDown = ROOT.TH2D('recomass_endcap_stretchDown', 'recomass_endcap_stretchDown', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_endcap_stretchDown)

        self.recomass_endcap_scaleUp = ROOT.TH2D('recomass_endcap_scaleUp', 'recomass_endcap_scaleUp', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_endcap_scaleUp)
        self.recomass_endcap_scaleDown = ROOT.TH2D('recomass_endcap_scaleDown', 'recomass_endcap_scaleDown', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_endcap_scaleDown)

        self.recomass_endcap_resUp = ROOT.TH2D('recomass_endcap_resUp', 'recomass_endcap_resUp', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_endcap_resUp)
        self.recomass_endcap_resDown = ROOT.TH2D('recomass_endcap_resDown', 'recomass_endcap_resDown', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomass_endcap_resDown)



        self.recomassprime = ROOT.TH2D('recomassprime', 'recomassprime', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime)

        self.recomassprime_shiftUp = ROOT.TH2D('recomassprime_shiftUp', 'recomassprime_shiftUp', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_shiftUp)
        self.recomassprime_shiftDown = ROOT.TH2D('recomassprime_shiftDown', 'recomassprime_shiftDown', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_shiftDown)

        self.recomassprime_stretchUp = ROOT.TH2D('recomassprime_stretchUp', 'recomassprime_stretchUp', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_stretchUp)
        self.recomassprime_stretchDown = ROOT.TH2D('recomassprime_stretchDown', 'recomassprime_stretchDown', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_stretchDown)

        self.recomassprime_scaleUp = ROOT.TH2D('recomassprime_scaleUp', 'recomassprime_scaleUp', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_scaleUp)
        self.recomassprime_scaleDown = ROOT.TH2D('recomassprime_scaleDown', 'recomassprime_scaleDown', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_scaleDown)

        self.recomassprime_resUp = ROOT.TH2D('recomassprime_resUp', 'recomassprime_resUp', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_resUp)
        self.recomassprime_resDown = ROOT.TH2D('recomassprime_resDown', 'recomassprime_resDown', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_resDown)

        self.recomassprime_barrel = ROOT.TH2D('recomassprime_barrel', 'recomassprime_barrel', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_barrel)

        self.recomassprime_barrel_shiftUp = ROOT.TH2D('recomassprime_barrel_shiftUp', 'recomassprime_barrel_shiftUp', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_barrel_shiftUp)
        self.recomassprime_barrel_shiftDown = ROOT.TH2D('recomassprime_barrel_shiftDown', 'recomassprime_barrel_shiftDown', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_barrel_shiftDown)

        self.recomassprime_barrel_stretchUp = ROOT.TH2D('recomassprime_barrel_stretchUp', 'recomassprime_barrel_stretchUp', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_barrel_stretchUp)
        self.recomassprime_barrel_stretchDown = ROOT.TH2D('recomassprime_barrel_stretchDown', 'recomassprime_barrel_stretchDown', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_barrel_stretchDown)

        self.recomassprime_barrel_scaleUp = ROOT.TH2D('recomassprime_barrel_scaleUp', 'recomassprime_barrel_scaleUp', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_barrel_scaleUp)
        self.recomassprime_barrel_scaleDown = ROOT.TH2D('recomassprime_barrel_scaleDown', 'recomassprime_barrel_scaleDown', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_barrel_scaleDown)

        self.recomassprime_barrel_resUp = ROOT.TH2D('recomassprime_barrel_resUp', 'recomassprime_barrel_resUp', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_barrel_resUp)
        self.recomassprime_barrel_resDown = ROOT.TH2D('recomassprime_barrel_resDown', 'recomassprime_barrel_resDown', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_barrel_resDown)

        self.recomassprime_endcap = ROOT.TH2D('recomassprime_endcap', 'recomassprime_endcap', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_endcap)

        self.recomassprime_endcap_shiftUp = ROOT.TH2D('recomassprime_endcap_shiftUp', 'recomassprime_endcap_shiftUp', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_endcap_shiftUp)
        self.recomassprime_endcap_shiftDown = ROOT.TH2D('recomassprime_endcap_shiftDown', 'recomassprime_endcap_shiftDown', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_endcap_shiftDown)

        self.recomassprime_endcap_stretchUp = ROOT.TH2D('recomassprime_endcap_stretchUp', 'recomassprime_endcap_stretchUp', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_endcap_stretchUp)
        self.recomassprime_endcap_stretchDown = ROOT.TH2D('recomassprime_endcap_stretchDown', 'recomassprime_endcap_stretchDown', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_endcap_stretchDown)

        self.recomassprime_endcap_scaleUp = ROOT.TH2D('recomassprime_endcap_scaleUp', 'recomassprime_endcap_scaleUp', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_endcap_scaleUp)
        self.recomassprime_endcap_scaleDown = ROOT.TH2D('recomassprime_endcap_scaleDown', 'recomassprime_endcap_scaleDown', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_endcap_scaleDown)

        self.recomassprime_endcap_resUp = ROOT.TH2D('recomassprime_endcap_resUp', 'recomassprime_endcap_resUp', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_endcap_resUp)
        self.recomassprime_endcap_resDown = ROOT.TH2D('recomassprime_endcap_resDown', 'recomassprime_endcap_resDown', len(OMEGA_BINS)-1, array('d', OMEGA_BINS), len(PHI_BINS)-1, array('d', PHI_BINS))
        self.addObject(self.recomassprime_endcap_resDown)


        self.hthat_gjets = ROOT.TH1F('GJETS_hthat_lhe', 'hthat', 150, 0, 1500)
        self.addObject(self.hthat_gjets)
        self.hthat_qcd = ROOT.TH1F('QCD_hthat_lhe', 'hthat', 300, 0, 3000)
        self.addObject(self.hthat_qcd)

        self.book_histo('SIGNAL_decaymode', 21, -1, 20, title='decaymode')
        self.book_histo('SIGNAL_prongs', 5, -1, 4, title='prongs')
        self.book_histo('SIGNAL_tag_pt_NUMER', 150, 0, 1500, title='tag_pt_NUMER')
        self.book_histo('SIGNAL_tag_pt_DENOM', 150, 0, 1500, title='tag_pt_DENOM')
        self.book_histo('SIGNAL_tag_eta_NUMER', 100, -5, 5, title='tag_eta_NUMER')
        self.book_histo('SIGNAL_tag_eta_DENOM', 100, -5, 5, title='tag_eta_DENOM')
        self.book_histo('SIGNAL_tag_phi_NUMER', 70, -3.5, 3.5, title='tag_phi_NUMER')
        self.book_histo('SIGNAL_tag_phi_DENOM', 70, -3.5, 3.5, title='tag_phi_DENOM')
        self.book_histo('SIGNAL_tag_npv_NUMER', 100, 0, 100, title='tag_npv_NUMER')
        self.book_histo('SIGNAL_tag_npv_DENOM', 100, 0, 100, title='tag_npv_DENOM')

        self.book_histo('SIGNAL_photontag_pt_NUMER', 150, 0, 1500, title='photontag_pt_NUMER')
        self.book_histo('SIGNAL_photontag_pt_DENOM', 150, 0, 1500, title='photontag_pt_DENOM')
        self.book_histo('SIGNAL_photontag_eta_NUMER', 100, -5, 5, title='photontag_eta_NUMER')
        self.book_histo('SIGNAL_photontag_eta_DENOM', 100, -5, 5, title='photontag_eta_DENOM')
        self.book_histo('SIGNAL_photontag_phi_NUMER', 70, -3.5, 3.5, title='photontag_phi_NUMER')
        self.book_histo('SIGNAL_photontag_phi_DENOM', 70, -3.5, 3.5, title='photontag_phi_DENOM')
        self.book_histo('SIGNAL_photontag_npv_NUMER', 100, 0, 100, title='photontag_npv_NUMER')
        self.book_histo('SIGNAL_photontag_npv_DENOM', 100, 0, 100, title='photontag_npv_DENOM')

        self.cutflow = ROOT.TH1F('cutflow', 'cutflow', 10, 0, 10)
        self.addObject(self.cutflow)

    def analyze(self, event):
        
        # get collections
        flags = Object(event, "Flag")
        twoprongs = Collection(event, "TwoProng")
        twoprongs = sorted(twoprongs, reverse=True, key=lambda obj : obj.pt) # FIXME sorting should be done upstream
        if self.photon == 'HPID':
          photons = Collection(event, "HighPtIdPhoton")
          recophi = Object(event, "HPID_RecoPhi")
          region = event.HPID_Region
          njets = event.HPID_NJets
          ht = event.HPID_HT
        elif self.photon == 'CBL':
          photons = Collection(event, "Photon")
          recophi = Object(event, "CBL_RecoPhi")
          region = event.CBL_Region
          njets = event.CBL_NJets
          ht = event.CBL_HT

        if self.year == "18": pass_trigger = event.HLT_Photon200
        elif self.year == "17": pass_trigger = event.HLT_Photon200
        elif self.year == "16": pass_trigger = event.HLT_Photon175
        else: pass_trigger = True
        ntwoprong = 0
        for twoprong in twoprongs:
          try:
            tight = twoprong.isTight
          except RuntimeError:
            tight = True
          if tight: ntwoprong += 1
        pass_UL18_HEM = True
        if region == 1:
          the_photon = get_vec(photons[recophi.photonindex])
          the_twoprong = get_vec(twoprongs[recophi.twoprongindex])
          recophi_vec = get_vec(twoprongs[recophi.twoprongindex]) + get_vec(photons[recophi.photonindex])
          # HEM correction
          if self.year == "18":
            if the_twoprong.Phi() > -1.57 and the_twoprong.Phi() < -0.87 and \
              the_twoprong.Eta() > -3.0 and the_twoprong.Eta() < -1.4: pass_UL18_HEM = False

        # filter decision for data
        '''
        if self.datamc == 'data':
          pass_filters = False
          if flags.goodVertices and \
            flags.globalSuperTightHalo2016Filter and \
            flags.HBHENoiseFilter and \
            flags.HBHENoiseIsoFilter and \
            flags.EcalDeadCellTriggerPrimitiveFilter and \
            flags.BadPFMuonFilter and \
            flags.eeBadScFilter and \
            flags.ecalBadCalibFilter:
              pass_filters = True
        '''

        pass_analysis_selection = (region == 1 and the_photon.Pt() > self.photon_min_pt and abs(photons[recophi.photonindex].scEta) < 1.4442 and abs(the_photon.Eta() - the_twoprong.Eta() < 1.4 and pass_UL18_HEM))

        # Fill
        self.cutflow.Fill(0)
        if self.datamc == 'mc':
          try:
            self.hthat_gjets.Fill(event.htHat_lhe, self.weight)
            self.hthat_qcd.Fill(event.htHat_lhe, self.weight)
          except RuntimeError:
            pass

        if region == 1:
          self.cutflow.Fill(1)

        if region == 1 and pass_analysis_selection:
          self.cutflow.Fill(2)

        if region == 1 and pass_analysis_selection and pass_trigger:
          self.cutflow.Fill(3)
          if self.photon == 'HPID':
            if abs(photons[recophi.photonindex].scEta)<1.4442: photon_subdet = 'barrel'
            elif abs(photons[recophi.photonindex].scEta)>1.566 and abs(photons[recophi.photonindex].scEta)<2.5: photon_subdet = 'endcap'
            else: photon_subdet = 'other'
          if self.photon == 'CBL':
            if photons[recophi.photonindex].isScEtaEB: photon_subdet = 'barrel'
            elif photons[recophi.photonindex].isScEtaEE: photon_subdet = 'endcap'
            else: photon_subdet = 'other'
          if abs(twoprongs[recophi.twoprongindex].eta)<1.4442: twoprong_subdet = 'barrel'
          else: twoprong_subdet = 'endcap'
          self.recophi_pt.Fill(recophi.pt, self.weight)
          self.recophi_eta.Fill(recophi.eta, self.weight)
          self.recophi_phi.Fill(recophi.phi, self.weight)
          self.recophi_m.Fill(recophi.mass, self.weight)
          self.recophi_dphi.Fill(abs(dPhi(the_photon,the_twoprong)), self.weight)
          self.recophi_deta.Fill(abs(the_photon.Eta() - the_twoprong.Eta()), self.weight)
          self.recophi_dr.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), self.weight)
          if njets == 0: self.recophi_dr_0jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), self.weight)
          if njets == 1: self.recophi_dr_1jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), self.weight)
          if njets >= 2: self.recophi_dr_2jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), self.weight)
          if photon_subdet == 'barrel' and twoprong_subdet == 'barrel':
            self.recophi_B_B_pt.Fill(recophi.pt, self.weight)
            self.recophi_B_B_eta.Fill(recophi.eta, self.weight)
            self.recophi_B_B_phi.Fill(recophi.phi, self.weight)
            self.recophi_B_B_m.Fill(recophi.mass, self.weight)
            self.recophi_B_B_dphi.Fill(abs(dPhi(the_photon,the_twoprong)), self.weight)
            self.recophi_B_B_deta.Fill(abs(the_photon.Eta() - the_twoprong.Eta()), self.weight)
            self.recophi_B_B_dr.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), self.weight)
            if njets == 0: self.recophi_B_B_dr_0jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), self.weight)
            if njets == 1: self.recophi_B_B_dr_1jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), self.weight)
            if njets >= 2: self.recophi_B_B_dr_2jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), self.weight)
          if photon_subdet == 'barrel' and twoprong_subdet == 'endcap':
            self.recophi_B_E_pt.Fill(recophi.pt, self.weight)
            self.recophi_B_E_eta.Fill(recophi.eta, self.weight)
            self.recophi_B_E_phi.Fill(recophi.phi, self.weight)
            self.recophi_B_E_m.Fill(recophi.mass, self.weight)
            self.recophi_B_E_dphi.Fill(abs(dPhi(the_photon,the_twoprong)), self.weight)
            self.recophi_B_E_deta.Fill(abs(the_photon.Eta() - the_twoprong.Eta()), self.weight)
            self.recophi_B_E_dr.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), self.weight)
            if njets == 0: self.recophi_B_E_dr_0jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), self.weight)
            if njets == 1: self.recophi_B_E_dr_1jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), self.weight)
            if njets >= 2: self.recophi_B_E_dr_2jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), self.weight)
          '''
          if photon_subdet == 'endcap' and twoprong_subdet == 'barrel':
            self.recophi_E_B_pt.Fill(recophi.pt, self.weight)
            self.recophi_E_B_eta.Fill(recophi.eta, self.weight)
            self.recophi_E_B_phi.Fill(recophi.phi, self.weight)
            self.recophi_E_B_m.Fill(recophi.mass, self.weight)
            self.recophi_E_B_dphi.Fill(abs(ROOT.Math.VectorUtil.DeltaPhi(the_photon,the_twoprong)), self.weight)
            self.recophi_E_B_deta.Fill(abs(the_photon.Eta() - the_twoprong.Eta()), self.weight)
            self.recophi_E_B_dr.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), self.weight)
            if event.NJets == 0: self.recophi_E_B_dr_0jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), self.weight)
            if event.NJets == 1: self.recophi_E_B_dr_1jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), self.weight)
            if event.NJets >= 2: self.recophi_E_B_dr_2jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), self.weight)
          if photon_subdet == 'endcap' and twoprong_subdet == 'endcap':
            self.recophi_E_E_pt.Fill(recophi.pt, self.weight)
            self.recophi_E_E_eta.Fill(recophi.eta, self.weight)
            self.recophi_E_E_phi.Fill(recophi.phi, self.weight)
            self.recophi_E_E_m.Fill(recophi.mass, self.weight)
            self.recophi_E_E_dphi.Fill(abs(ROOT.Math.VectorUtil.DeltaPhi(the_photon,the_twoprong)), self.weight)
            self.recophi_E_E_deta.Fill(abs(the_photon.Eta() - the_twoprong.Eta()), self.weight)
            self.recophi_E_E_dr.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), self.weight)
            if event.NJets == 0: self.recophi_E_E_dr_0jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), self.weight)
            if event.NJets == 1: self.recophi_E_E_dr_1jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), self.weight)
            if event.NJets >= 2: self.recophi_E_E_dr_2jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), self.weight)
          '''
          self.photon_pt.Fill(the_photon.Pt(), self.weight)
            
          # scaling and smearing for photon_pt
          #corr_photon = photons[recophi.photonindex]
        
          #corr_up = corr_photon.dEscaleUp + corr_photon.dEsigmaUp * corr_photon.eCorr
          #photon_pt_up = corr_photon.pt / corr_photon.eCorr * corr_up
          #corr_down = corr_photon.dEscaledown + corr_photon.dEsigmadown * corr_photon.eCorr
          #photon_pt_down = corr_photon.pt / corr_photon.eCorr * corr_down
          
          #self.photon_pt_up.Fill(photon_pt_up, self.weight)
          #self.photon_pt_down.Fill(photon_pt_down, self.weight)

          # Create photon pt scaled hists
          #scaleup_pt = sel_photon.pt * sel_photon.dEscaleUp
          #scaledown_pt = sel_photon.pt * sel_photon.dEscaleDown

          #self.photon_scaleup_sigmaup.Fill(scaleup_pt, self.weight)
          #self.photon_scaledown_sigmadown.Fill(scaledown_pt, self.weight)
          #self.photon_unscaled.Fill(sel_photon.pt, self.weight)
          
          self.photon_eta.Fill(the_photon.Eta(), self.weight)
          self.photon_phi.Fill(the_photon.Phi(), self.weight)
          self.twoprong_pt.Fill(the_twoprong.Pt(), self.weight)
          self.twoprong_eta.Fill(the_twoprong.Eta(), self.weight)
          self.twoprong_phi.Fill(the_twoprong.Phi(), self.weight)
          self.twoprong_mass.Fill(twoprongs[recophi.twoprongindex].mass, self.weight)
          self.twoprong_masspi0.Fill(twoprongs[recophi.twoprongindex].massPi0, self.weight)
          self.twoprong_masseta.Fill(twoprongs[recophi.twoprongindex].massEta, self.weight)
          self.nphoton.Fill(len(photons), self.weight)
          self.ntwoprong.Fill(ntwoprong, self.weight)
          self.njets.Fill(njets, self.weight)
          self.ht.Fill(ht, self.weight)
          self.met.Fill(event.MET_pt, self.weight)
          self.met_phi.Fill(event.MET_phi, self.weight)
          if event.MET_pt > 30: self.met_phi_cut.Fill(event.MET_phi, self.weight)
          self.npv.Fill(event.PV_npvs, self.weight)
          if self.photon == 'HPID':
            if abs(photons[recophi.photonindex].scEta)<1.4442: self.twoprong_eta_barrel.Fill(the_twoprong.Eta(), self.weight)
            if abs(photons[recophi.photonindex].scEta)>1.566 and abs(photons[recophi.photonindex].scEta)<2.5: self.twoprong_eta_endcap.Fill(twoprong.Eta(), self.weight)
          if self.photon == 'CBL':
            if photons[recophi.photonindex].isScEtaEB: self.twoprong_eta_barrel.Fill(the_twoprong.Eta(), self.weight)
            if photons[recophi.photonindex].isScEtaEE: self.twoprong_eta_endcap.Fill(the_twoprong.Eta(), self.weight)
          self.recomass_uniformbinning.Fill(twoprongs[recophi.twoprongindex].massPi0, recophi.mass, self.weight)
          self.recomass_bare.Fill(twoprongs[recophi.twoprongindex].massPi0, recophi.mass, self.weight)
          if twoprong_subdet == 'barrel': self.recomass_bare_barrel.Fill(twoprongs[recophi.twoprongindex].massPi0, recophi.mass, self.weight)
          if twoprong_subdet == 'endcap': self.recomass_bare_endcap.Fill(twoprongs[recophi.twoprongindex].massPi0, recophi.mass, self.weight)

          self.recomassprime_bare.Fill(twoprongs[recophi.twoprongindex].massEta, recophi.mass, self.weight)
          if twoprong_subdet == 'barrel': self.recomassprime_bare_barrel.Fill(twoprongs[recophi.twoprongindex].massEta, recophi.mass, self.weight)
          if twoprong_subdet == 'endcap': self.recomassprime_bare_endcap.Fill(twoprongs[recophi.twoprongindex].massEta, recophi.mass, self.weight)

          m2p = twoprongs[recophi.twoprongindex].massPi0
          m2peta = twoprongs[recophi.twoprongindex].massEta
          pt2p = twoprongs[recophi.twoprongindex].pt
          m2p_gen = M2P_GEN
          tp_vec = get_vec(twoprongs[recophi.twoprongindex])

          eff_central = 0.9
          eff_up = 0.9 + 0.22
          eff_down = 0.9 - 0.22

          shift_central = -0.01
          shift_up = -0.01 + 0.03
          shift_down = -0.01 - 0.03

          stretch_central = 1.12
          stretch_up = 1.12 + 0.04
          stretch_down = 1.12 - 0.04

          scale_central = 1.07
          scale_up = 1.07 + 0.08
          scale_down = 1.07 - 0.08

          res_central = 1.07
          res_up = 1.07 + 0.07
          res_down = 1.07 - 0.07

          m2p_central = (m2p - m2p_gen + shift_central) * stretch_central + m2p_gen
          m2p_shift_up = (m2p - m2p_gen + shift_up) * stretch_central + m2p_gen
          m2p_shift_down = (m2p - m2p_gen + shift_down) * stretch_central + m2p_gen
          m2p_stretch_up = (m2p - m2p_gen + shift_central) * stretch_up + m2p_gen
          m2p_stretch_down = (m2p - m2p_gen + shift_central) * stretch_down + m2p_gen

          m2peta_central = (m2peta - m2p_gen + shift_central) * stretch_central + m2p_gen
          m2peta_shift_up = (m2peta - m2p_gen + shift_up) * stretch_central + m2p_gen
          m2peta_shift_down = (m2peta - m2p_gen + shift_down) * stretch_central + m2p_gen
          m2peta_stretch_up = (m2peta - m2p_gen + shift_central) * stretch_up + m2p_gen
          m2peta_stretch_down = (m2peta - m2p_gen + shift_central) * stretch_down + m2p_gen

          self.recomass.Fill(m2p_central, recophi.mass, self.weight)
          self.recomass_shiftUp.Fill(m2p_shift_up, recophi.mass, self.weight)
          self.recomass_shiftDown.Fill(m2p_shift_down, recophi.mass, self.weight)
          self.recomass_stretchUp.Fill(m2p_stretch_up, recophi.mass, self.weight)
          self.recomass_stretchDown.Fill(m2p_stretch_down, recophi.mass, self.weight)

          self.recomassprime.Fill(m2peta_central, recophi.mass, self.weight)
          self.recomassprime_shiftUp.Fill(m2peta_shift_up, recophi.mass, self.weight)
          self.recomassprime_shiftDown.Fill(m2peta_shift_down, recophi.mass, self.weight)
          self.recomassprime_stretchUp.Fill(m2peta_stretch_up, recophi.mass, self.weight)
          self.recomassprime_stretchDown.Fill(m2peta_stretch_down, recophi.mass, self.weight)

          if twoprong_subdet == 'barrel': 
              self.recomass_barrel.Fill(m2p_central, recophi.mass, self.weight)
              self.recomass_barrel_shiftUp.Fill(m2p_shift_up, recophi.mass, self.weight)
              self.recomass_barrel_shiftDown.Fill(m2p_shift_down, recophi.mass, self.weight)
              self.recomass_barrel_stretchUp.Fill(m2p_stretch_up, recophi.mass, self.weight)
              self.recomass_barrel_stretchDown.Fill(m2p_stretch_down, recophi.mass, self.weight)

              self.recomassprime_barrel.Fill(m2peta_central, recophi.mass, self.weight)
              self.recomassprime_barrel_shiftUp.Fill(m2peta_shift_up, recophi.mass, self.weight)
              self.recomassprime_barrel_shiftDown.Fill(m2peta_shift_down, recophi.mass, self.weight)
              self.recomassprime_barrel_stretchUp.Fill(m2peta_stretch_up, recophi.mass, self.weight)
              self.recomassprime_barrel_stretchDown.Fill(m2peta_stretch_down, recophi.mass, self.weight)

          if twoprong_subdet == 'endcap':
              self.recomass_endcap.Fill(m2p_central, recophi.mass, self.weight)
              self.recomass_endcap_shiftUp.Fill(m2p_shift_up, recophi.mass, self.weight)
              self.recomass_endcap_shiftDown.Fill(m2p_shift_down, recophi.mass, self.weight)
              self.recomass_endcap_stretchUp.Fill(m2p_stretch_up, recophi.mass, self.weight)
              self.recomass_endcap_stretchDown.Fill(m2p_stretch_down, recophi.mass, self.weight)

              self.recomassprime_endcap.Fill(m2peta_central, recophi.mass, self.weight)
              self.recomassprime_endcap_shiftUp.Fill(m2peta_shift_up, recophi.mass, self.weight)
              self.recomassprime_endcap_shiftDown.Fill(m2peta_shift_down, recophi.mass, self.weight)
              self.recomassprime_endcap_stretchUp.Fill(m2peta_stretch_up, recophi.mass, self.weight)
              self.recomassprime_endcap_stretchDown.Fill(m2peta_stretch_down, recophi.mass, self.weight)

          N_central = self.rand.Gaus(0, res_central - 1)
          N_up = self.rand.Gaus(0, res_up - 1)
          N_down = self.rand.Gaus(0, res_down - 1)
          pt2p_central = pt2p * scale_central + pt2p * N_central
          pt2p_scale_up = pt2p * scale_up + pt2p * N_central
          pt2p_scale_down = pt2p * scale_down + pt2p * N_central
          pt2p_res_up = pt2p * scale_central + pt2p * N_up
          pt2p_res_down = pt2p * scale_central + pt2p * N_down

          tp_vec_res_up = clone_vec(tp_vec)
          tp_vec_res_up.SetPt(pt2p_res_up)
          tp_vec_res_down = clone_vec(tp_vec)
          tp_vec_res_down.SetPt(pt2p_res_down)
          tp_vec_scale_up = clone_vec(tp_vec)
          tp_vec_scale_up.SetPt(pt2p_scale_up)
          tp_vec_scale_down = clone_vec(tp_vec)
          tp_vec_scale_down.SetPt(pt2p_scale_down)

          recophi_vec_res_up = tp_vec_res_up + get_vec(photons[recophi.photonindex])
          recophi_vec_res_down = tp_vec_res_down + get_vec(photons[recophi.photonindex])
          recophi_vec_scale_up = tp_vec_scale_up + get_vec(photons[recophi.photonindex])
          recophi_vec_scale_down = tp_vec_scale_down + get_vec(photons[recophi.photonindex])
          
          self.recomass_resUp.Fill(m2p, recophi_vec_res_up.M(), self.weight)
          self.recomass_resDown.Fill(m2p, recophi_vec_res_down.M(), self.weight)
          self.recomass_scaleUp.Fill(m2p, recophi_vec_scale_up.M(), self.weight)
          self.recomass_scaleDown.Fill(m2p, recophi_vec_scale_down.M(), self.weight)

          self.recomassprime_resUp.Fill(m2peta, recophi_vec_res_up.M(), self.weight)
          self.recomassprime_resDown.Fill(m2peta, recophi_vec_res_down.M(), self.weight)
          self.recomassprime_scaleUp.Fill(m2peta, recophi_vec_scale_up.M(), self.weight)
          self.recomassprime_scaleDown.Fill(m2peta, recophi_vec_scale_down.M(), self.weight)

          if twoprong_subdet == 'barrel': 
              self.recomass_barrel_resUp.Fill(m2p, recophi_vec_res_up.M(), self.weight)
              self.recomass_barrel_resDown.Fill(m2p, recophi_vec_res_down.M(), self.weight)
              self.recomass_barrel_scaleUp.Fill(m2p, recophi_vec_scale_up.M(), self.weight)
              self.recomass_barrel_scaleDown.Fill(m2p, recophi_vec_scale_down.M(), self.weight)

              self.recomassprime_barrel_resUp.Fill(m2peta, recophi_vec_res_up.M(), self.weight)
              self.recomassprime_barrel_resDown.Fill(m2peta, recophi_vec_res_down.M(), self.weight)
              self.recomassprime_barrel_scaleUp.Fill(m2peta, recophi_vec_scale_up.M(), self.weight)
              self.recomassprime_barrel_scaleDown.Fill(m2peta, recophi_vec_scale_down.M(), self.weight)

          if twoprong_subdet == 'endcap': 
              self.recomass_endcap_resUp.Fill(m2p, recophi_vec_res_up.M(), self.weight)
              self.recomass_endcap_resDown.Fill(m2p, recophi_vec_res_down.M(), self.weight)
              self.recomass_endcap_scaleUp.Fill(m2p, recophi_vec_scale_up.M(), self.weight)
              self.recomass_endcap_scaleDown.Fill(m2p, recophi_vec_scale_down.M(), self.weight)

              self.recomassprime_endcap_resUp.Fill(m2peta, recophi_vec_res_up.M(), self.weight)
              self.recomassprime_endcap_resDown.Fill(m2peta, recophi_vec_res_down.M(), self.weight)
              self.recomassprime_endcap_scaleUp.Fill(m2peta, recophi_vec_scale_up.M(), self.weight)
              self.recomassprime_endcap_scaleDown.Fill(m2peta, recophi_vec_scale_down.M(), self.weight)

        tight_photons = []
        loose_photons = []
        for photon in photons:
            if photon.pt > 220 and photon.isScEtaEB and photon.hadTowOverEm < 0.04596: # photon cuts
                if photon.cutBased >= 1:
                    tight_photons.append(photon)
                    break
                else: loose_photons.append(photon)
        if len(loose_photons) == 0 and len(tight_photons) == 0: return False
        if not len(tight_photons) == 0:
          sel_photon = tight_photons[0]
          tight_photon = True
        else: sel_photon = loose_photons[0]
        tight_photon = False
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
        else: return False
        if abs(sel_tp.eta)<1.4442: twoprong_subdet = 'barrel'
        else: twoprong_subdet = 'endcap'
        if len(tight_photons) >= 1 and pass_trigger and len(noniso_asym_tp) != 0 and len(iso_sym_tp) == 0 and len(iso_asym_tp) == 0 and len(noniso_sym_tp) == 0:  # tight photon and loose twoprong (asym and noniso)
          self.cutflow.Fill(5)
          self.recomass_sideband.Fill(sel_tp.massPi0, recophi.mass, self.weight)
          if twoprong_subdet == 'barrel': self.recomass_sideband_barrel.Fill(sel_tp.massPi0, recophi.mass, self.weight)
          if twoprong_subdet == 'endcap': self.recomass_sideband_endcap.Fill(sel_tp.massPi0, recophi.mass, self.weight)

          self.recomassprime_sideband.Fill(sel_tp.massEta, recophi.mass, self.weight)
          if twoprong_subdet == 'barrel': self.recomassprime_sideband_barrel.Fill(sel_tp.massEta, recophi.mass, self.weight)
          if twoprong_subdet == 'endcap': self.recomassprime_sideband_endcap.Fill(sel_tp.massEta, recophi.mass, self.weight)

        if len(loose_photons) >= 1 and len(tight_photons) == 0 and pass_trigger and len(iso_asym_tp) != 0:  # loose photon and loose twoprong
          pass

        # signal
        if self.datamc == 'sigRes':
          genomegas = Collection(event, "GenOmega")
          for i in range(event.nGenOmega):
            self.SIGNAL_decaymode.Fill(event.GenOmega_decaymode[i])
            self.SIGNAL_prongs.Fill(event.GenOmega_prongs[i])
          
          for genomega in genomegas:
            genvec = get_vec(genomega)
            tagged = False
            for twoprong in twoprongs:
              twoprongvec = get_vec(twoprong)
              if ROOT.Math.VectorUtil.DeltaR(genvec, twoprongvec) < 0.1:
                tagged = True
                break
            self.SIGNAL_tag_pt_DENOM.Fill(genvec.Pt())
            self.SIGNAL_tag_eta_DENOM.Fill(genvec.Eta())
            self.SIGNAL_tag_phi_DENOM.Fill(genvec.Phi())
            self.SIGNAL_tag_npv_DENOM.Fill(event.PV_npvs)
            if tagged: self.SIGNAL_tag_pt_NUMER.Fill(genvec.Pt())
            if tagged: self.SIGNAL_tag_eta_NUMER.Fill(genvec.Eta())
            if tagged: self.SIGNAL_tag_phi_NUMER.Fill(genvec.Phi())
            if tagged: self.SIGNAL_tag_npv_NUMER.Fill(event.PV_npvs)
            photon_tagged = False
            if self.photon == 'CBL': 
              photon_tagged = False
              for photon in photons:
                if photon.cutBased < PHOTON_CUTBASED_ID: continue
                if abs(photon.scEta) > PHOTON_BARREL_ETA: continue
                if photon.pt <= self.photon_min_pt: continue
                if photon.hadTowOverEm > PHOTON_HoverE_CUT: continue
                photonvec = get_vec(photon)
                if ROOT.Math.VectorUtil.DeltaR(genvec, photonvec) < 0.1:
                  photon_tagged = True
                  break
            if self.photon == 'HPID': 
              photon_tagged = False
              for photon in photons:
                photonvec = get_vec(photon)
                if ROOT.Math.VectorUtil.DeltaR(genvec, photonvec) < 0.1:
                  photon_tagged = True
                  break
            self.SIGNAL_photontag_pt_DENOM.Fill(genvec.Pt())
            self.SIGNAL_photontag_eta_DENOM.Fill(genvec.Eta())
            self.SIGNAL_photontag_phi_DENOM.Fill(genvec.Phi())
            self.SIGNAL_photontag_npv_DENOM.Fill(event.PV_npvs)
            if photon_tagged: self.SIGNAL_photontag_pt_NUMER.Fill(genvec.Pt())
            if photon_tagged: self.SIGNAL_photontag_eta_NUMER.Fill(genvec.Eta())
            if photon_tagged: self.SIGNAL_photontag_phi_NUMER.Fill(genvec.Phi())
            if photon_tagged: self.SIGNAL_photontag_npv_NUMER.Fill(event.PV_npvs)


        return True
