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

class SanityAnalysis(Module):
    def __init__(self, datamc, lumi=1.0, dict_xs=None, dict_ngen=None, cut='None', photon='HPID'):
        self.writeHistFile = True
        self.lumi = lumi
        self.dict_xs = dict_xs
        self.dict_ngen = dict_ngen
        self.cut = cut

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

        phi_bins = [500, 520, 541, 563, 586, 609, 633, 658, 684, 711, 739, 769, 800, 832, 865, 900, 936, 973, 1012, 1052, 1094, 1138, 1184, 1231, 1280, 1331, 1384, 1439, 1497, 1557, 1619, 1684, 1751, 1821, 1894, 1970, 2049, 2131, 2216, 2305, 2397, 2493, 2593, 2697, 2805, 2917, 3034, 3155, 3281, 3412, 3548, 3690, 3838, 3992]
        omega_bins = [0.40, 0.53, 0.58, 0.64, 0.70, 0.77, 0.85, 0.94, 1.03, 1.13, 1.24, 1.36, 1.50, 1.65, 1.81, 1.99, 2.19, 2.41, 2.65, 2.92, 3.21, 3.5]

        self.recomass_2d = ROOT.TH2D('recomass_2d', 'recomass_2d', 300, 0, 15, 200, 0, 6000)
        self.addObject(self.recomass_2d)
        self.recomass_2d_variable = ROOT.TH2D('recomass_2d_variable', 'recomass_2d_variable', len(omega_bins)-1, array('d', omega_bins), len(phi_bins)-1, array('d', phi_bins))
        self.addObject(self.recomass_2d_variable)
        self.recomass_2d_variable_sb = ROOT.TH2D('recomass_2d_variable_sb', 'recomass_2d_variable_sb', len(omega_bins)-1, array('d', omega_bins), len(phi_bins)-1, array('d', phi_bins))
        self.addObject(self.recomass_2d_variable_sb)
        self.recomass_2d_variable_sb_tightpho = ROOT.TH2D('recomass_2d_variable_sb_tightpho', 'recomass_2d_variable_sb_tightpho', len(omega_bins)-1, array('d', omega_bins), len(phi_bins)-1, array('d', phi_bins))
        self.addObject(self.recomass_2d_variable_sb_tightpho)

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
        
        # weight
        if self.dict_xs and self.dict_ngen:
          dataset_id = event.dataset_id
          xs = self.dict_xs[dataset_id]
          Ngen = self.dict_ngen[dataset_id]
          weight = xs * self.lumi / Ngen
        else:
          weight = 1.0

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
            self.hthat_gjets.Fill(event.htHat_lhe, weight)
            self.hthat_qcd.Fill(event.htHat_lhe, weight)
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
          self.recophi_pt.Fill(recophi.pt, weight)
          self.recophi_eta.Fill(recophi.eta, weight)
          self.recophi_phi.Fill(recophi.phi, weight)
          self.recophi_m.Fill(recophi.mass, weight)
          self.recophi_dphi.Fill(abs(dPhi(the_photon,the_twoprong)), weight)
          self.recophi_deta.Fill(abs(the_photon.Eta() - the_twoprong.Eta()), weight)
          self.recophi_dr.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
          if njets == 0: self.recophi_dr_0jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
          if njets == 1: self.recophi_dr_1jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
          if njets >= 2: self.recophi_dr_2jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
          if photon_subdet == 'barrel' and twoprong_subdet == 'barrel':
            self.recophi_B_B_pt.Fill(recophi.pt, weight)
            self.recophi_B_B_eta.Fill(recophi.eta, weight)
            self.recophi_B_B_phi.Fill(recophi.phi, weight)
            self.recophi_B_B_m.Fill(recophi.mass, weight)
            self.recophi_B_B_dphi.Fill(abs(dPhi(the_photon,the_twoprong)), weight)
            self.recophi_B_B_deta.Fill(abs(the_photon.Eta() - the_twoprong.Eta()), weight)
            self.recophi_B_B_dr.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
            if njets == 0: self.recophi_B_B_dr_0jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
            if njets == 1: self.recophi_B_B_dr_1jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
            if njets >= 2: self.recophi_B_B_dr_2jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
          if photon_subdet == 'barrel' and twoprong_subdet == 'endcap':
            self.recophi_B_E_pt.Fill(recophi.pt, weight)
            self.recophi_B_E_eta.Fill(recophi.eta, weight)
            self.recophi_B_E_phi.Fill(recophi.phi, weight)
            self.recophi_B_E_m.Fill(recophi.mass, weight)
            self.recophi_B_E_dphi.Fill(abs(dPhi(the_photon,the_twoprong)), weight)
            self.recophi_B_E_deta.Fill(abs(the_photon.Eta() - the_twoprong.Eta()), weight)
            self.recophi_B_E_dr.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
            if njets == 0: self.recophi_B_E_dr_0jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
            if njets == 1: self.recophi_B_E_dr_1jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
            if njets >= 2: self.recophi_B_E_dr_2jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
          '''
          if photon_subdet == 'endcap' and twoprong_subdet == 'barrel':
            self.recophi_E_B_pt.Fill(recophi.pt, weight)
            self.recophi_E_B_eta.Fill(recophi.eta, weight)
            self.recophi_E_B_phi.Fill(recophi.phi, weight)
            self.recophi_E_B_m.Fill(recophi.mass, weight)
            self.recophi_E_B_dphi.Fill(abs(ROOT.Math.VectorUtil.DeltaPhi(the_photon,the_twoprong)), weight)
            self.recophi_E_B_deta.Fill(abs(the_photon.Eta() - the_twoprong.Eta()), weight)
            self.recophi_E_B_dr.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
            if event.NJets == 0: self.recophi_E_B_dr_0jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
            if event.NJets == 1: self.recophi_E_B_dr_1jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
            if event.NJets >= 2: self.recophi_E_B_dr_2jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
          if photon_subdet == 'endcap' and twoprong_subdet == 'endcap':
            self.recophi_E_E_pt.Fill(recophi.pt, weight)
            self.recophi_E_E_eta.Fill(recophi.eta, weight)
            self.recophi_E_E_phi.Fill(recophi.phi, weight)
            self.recophi_E_E_m.Fill(recophi.mass, weight)
            self.recophi_E_E_dphi.Fill(abs(ROOT.Math.VectorUtil.DeltaPhi(the_photon,the_twoprong)), weight)
            self.recophi_E_E_deta.Fill(abs(the_photon.Eta() - the_twoprong.Eta()), weight)
            self.recophi_E_E_dr.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
            if event.NJets == 0: self.recophi_E_E_dr_0jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
            if event.NJets == 1: self.recophi_E_E_dr_1jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
            if event.NJets >= 2: self.recophi_E_E_dr_2jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
          '''
          self.photon_pt.Fill(the_photon.Pt(), weight)
            
          # scaling and smearing for photon_pt
          #corr_photon = photons[recophi.photonindex]
        
          #corr_up = corr_photon.dEscaleUp + corr_photon.dEsigmaUp * corr_photon.eCorr
          #photon_pt_up = corr_photon.pt / corr_photon.eCorr * corr_up
          #corr_down = corr_photon.dEscaledown + corr_photon.dEsigmadown * corr_photon.eCorr
          #photon_pt_down = corr_photon.pt / corr_photon.eCorr * corr_down
          
          #self.photon_pt_up.Fill(photon_pt_up, weight)
          #self.photon_pt_down.Fill(photon_pt_down, weight)
          
          self.photon_eta.Fill(the_photon.Eta(), weight)
          self.photon_phi.Fill(the_photon.Phi(), weight)
          self.twoprong_pt.Fill(the_twoprong.Pt(), weight)
          self.twoprong_eta.Fill(the_twoprong.Eta(), weight)
          self.twoprong_phi.Fill(the_twoprong.Phi(), weight)
          self.twoprong_mass.Fill(twoprongs[recophi.twoprongindex].mass, weight)
          self.twoprong_masspi0.Fill(twoprongs[recophi.twoprongindex].massPi0, weight)
          self.twoprong_masseta.Fill(twoprongs[recophi.twoprongindex].massEta, weight)
          self.nphoton.Fill(len(photons), weight)
          self.ntwoprong.Fill(ntwoprong, weight)
          self.njets.Fill(njets, weight)
          self.ht.Fill(ht, weight)
          self.met.Fill(event.MET_pt, weight)
          self.met_phi.Fill(event.MET_phi, weight)
          if event.MET_pt > 30: self.met_phi_cut.Fill(event.MET_phi, weight)
          self.npv.Fill(event.PV_npvs, weight)
          if self.photon == 'HPID':
            if abs(photons[recophi.photonindex].scEta)<1.4442: self.twoprong_eta_barrel.Fill(the_twoprong.Eta(), weight)
            if abs(photons[recophi.photonindex].scEta)>1.566 and abs(photons[recophi.photonindex].scEta)<2.5: self.twoprong_eta_endcap.Fill(twoprong.Eta(), weight)
          if self.photon == 'CBL':
            if photons[recophi.photonindex].isScEtaEB: self.twoprong_eta_barrel.Fill(the_twoprong.Eta(), weight)
            if photons[recophi.photonindex].isScEtaEE: self.twoprong_eta_endcap.Fill(the_twoprong.Eta(), weight)
          self.recomass_2d.Fill(twoprongs[recophi.twoprongindex].massPi0, recophi.mass, weight)
          self.recomass_2d_variable.Fill(twoprongs[recophi.twoprongindex].massPi0, recophi.mass, weight)
          self.recomass_2d_variable.Fill(twoprongs[recophi.twoprongindex].massPi0, recophi.mass, weight)

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
        if len(loose_photons) >= 1 and len(tight_photons) == 0 and pass_trigger and len(iso_asym_tp) != 0:  # loose photon and loose twoprong
          self.cutflow.Fill(4)
          self.recomass_2d_variable_sb.Fill(sel_tp.massPi0, recophi.mass, weight)
        if len(tight_photons) >= 1 and pass_trigger and len(noniso_asym_tp) != 0 and len(iso_sym_tp) == 0 and len(iso_asym_tp) == 0 and len(noniso_sym_tp) == 0:  # tight photon and loose twoprong (asym and noniso)
          self.cutflow.Fill(5)
          self.recomass_2d_variable_sb_tightpho.Fill(sel_tp.massPi0, recophi.mass, weight)

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
