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

class MyAnalysis(Module):
    def __init__(self, datamc, lumi=1.0, dict_xs=None, dict_ngen=None, deta=False, photon='HPID'):
        self.writeHistFile = True
        self.lumi = lumi
        self.dict_xs = dict_xs
        self.dict_ngen = dict_ngen
        self.datamc = datamc
        self.deta = deta
        self.photon = photon

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

        self.recophi_dphi = ROOT.TH1F('recophi_dphi', 'recophi_dphi', 35, 0, 3.5)
        self.addObject(self.recophi_dphi)
        self.recophi_deta = ROOT.TH1F('recophi_deta', 'recophi_deta', 70, 0, 7)
        self.addObject(self.recophi_deta)
        self.recophi_dr = ROOT.TH1F('recophi_dr', 'recophi_dr', 70, 0, 7)
        self.addObject(self.recophi_dr)
        self.recophi_dr_0jet = ROOT.TH1F('recophi_dr_0jet', 'recophi_dr_0jet', 70, 0, 7)
        self.addObject(self.recophi_dr_0jet)
        self.recophi_dr_1jet = ROOT.TH1F('recophi_dr_1jet', 'recophi_dr_1jet', 70, 0, 7)
        self.addObject(self.recophi_dr_1jet)
        self.recophi_dr_2jet = ROOT.TH1F('recophi_dr_2jet', 'recophi_dr_2jet', 70, 0, 7)
        self.addObject(self.recophi_dr_2jet)
        self.recophi_m = ROOT.TH1F('recophi_m', 'recophi_m', 200, 0, 6000)
        self.addObject(self.recophi_m)
        self.recophi_pt = ROOT.TH1F('recophi_pt', 'recophi_pt', 200, 0, 2000)
        self.addObject(self.recophi_pt)
        self.recophi_eta = ROOT.TH1F('recophi_eta', 'recophi_eta', 200, -10, 10)
        self.addObject(self.recophi_eta)
        self.recophi_phi = ROOT.TH1F('recophi_phi', 'recophi_phi', 70, -3.5, 3.5)
        self.addObject(self.recophi_phi)

        self.recophi_b_dphi = ROOT.TH1F('recophi_b_dphi', 'recophi_b_dphi', 35, 0, 3.5)
        self.addObject(self.recophi_b_dphi)
        self.recophi_b_deta = ROOT.TH1F('recophi_b_deta', 'recophi_b_deta', 70, 0, 7)
        self.addObject(self.recophi_b_deta)
        self.recophi_b_dr = ROOT.TH1F('recophi_b_dr', 'recophi_b_dr', 70, 0, 7)
        self.addObject(self.recophi_b_dr)
        self.recophi_b_m = ROOT.TH1F('recophi_b_m', 'recophi_b_m', 200, 0, 6000)
        self.addObject(self.recophi_b_m)
        self.recophi_b_pt = ROOT.TH1F('recophi_b_pt', 'recophi_b_pt', 200, 0, 2000)
        self.addObject(self.recophi_b_pt)
        self.recophi_b_eta = ROOT.TH1F('recophi_b_eta', 'recophi_b_eta', 200, -10, 10)
        self.addObject(self.recophi_b_eta)
        self.recophi_b_phi = ROOT.TH1F('recophi_b_phi', 'recophi_b_phi', 70, -3.5, 3.5)
        self.addObject(self.recophi_b_phi)

        self.recophi_e_dphi = ROOT.TH1F('recophi_e_dphi', 'recophi_e_dphi', 35, 0, 3.5)
        self.addObject(self.recophi_e_dphi)
        self.recophi_e_deta = ROOT.TH1F('recophi_e_deta', 'recophi_e_deta', 70, 0, 7)
        self.addObject(self.recophi_e_deta)
        self.recophi_e_dr = ROOT.TH1F('recophi_e_dr', 'recophi_e_dr', 70, 0, 7)
        self.addObject(self.recophi_e_dr)
        self.recophi_e_m = ROOT.TH1F('recophi_e_m', 'recophi_e_m', 200, 0, 6000)
        self.addObject(self.recophi_e_m)
        self.recophi_e_pt = ROOT.TH1F('recophi_e_pt', 'recophi_e_pt', 200, 0, 2000)
        self.addObject(self.recophi_e_pt)
        self.recophi_e_eta = ROOT.TH1F('recophi_e_eta', 'recophi_e_eta', 200, -10, 10)
        self.addObject(self.recophi_e_eta)
        self.recophi_e_phi = ROOT.TH1F('recophi_e_phi', 'recophi_e_phi', 70, -3.5, 3.5)
        self.addObject(self.recophi_e_phi)

        self.photon_pt = ROOT.TH1F('photon_pt', 'photon_pt', 160, 0, 1600)
        self.addObject(self.photon_pt)
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

        self.hthat_gjets = ROOT.TH1F('GJETS_hthat_lhe', 'hthat', 150, 0, 1500)
        self.addObject(self.hthat_gjets)
        self.hthat_qcd = ROOT.TH1F('QCD_hthat_lhe', 'hthat', 300, 0, 3000)
        self.addObject(self.hthat_qcd)

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
        twoprongs = Collection(event, "TwoProng")
        # FIXME sorting should be done upstream
        twoprongs = sorted(twoprongs, reverse=True, key=lambda obj : obj.pt)
        if self.photon == 'HPID':
          photons = Collection(event, "HighPtIdPhoton")
          recophi = Object(event, "RecoPhi")
          region = event.Region
        elif self.photon == 'CBL':
          photons = Collection(event, "Photon")
          recophi = Object(event, "CutBased_RecoPhi")
          region = event.CutBased_Region
        pass_trigger = event.HLT_Photon200
        ntwoprong = 0
        for twoprong in twoprongs:
          try:
            tight = twoprong.isTight
          except RuntimeError:
            tight = True
          if tight: ntwoprong += 1
        if region == 1:
          photon = get_vec(photons[recophi.photonindex])
          twoprong = get_vec(twoprongs[recophi.twoprongindex])
          deta = abs(photon.Eta() - twoprong.Eta())

        # deta cut
        if region == 1:
          if self.deta and deta < 1.5: pass_deta = True
          elif self.deta and deta > 1.5: pass_deta = False
          else: pass_deta = True
        
        self.cutflow.Fill(0)
        try:
          if len(photons)>=1:
            self.hthat_gjets.Fill(event.htHat_lhe, weight)
            self.hthat_qcd.Fill(event.htHat_lhe, weight)
        except RuntimeError:
          pass

        if region == 1:
          self.cutflow.Fill(1)

        if region == 1 and photon.Pt() > 220 and pass_deta:
          self.cutflow.Fill(2)

        if region == 1 and photon.Pt() > 220 and pass_trigger and pass_deta:
          self.cutflow.Fill(3)
          self.recophi_pt.Fill(recophi.pt, weight)
          self.recophi_eta.Fill(recophi.eta, weight)
          self.recophi_phi.Fill(recophi.phi, weight)
          self.recophi_m.Fill(recophi.mass, weight)
          self.recophi_dphi.Fill(abs(ROOT.Math.VectorUtil.DeltaPhi(photon,twoprong)), weight)
          self.recophi_deta.Fill(abs(photon.Eta() - twoprong.Eta()), weight)
          self.recophi_dr.Fill(ROOT.Math.VectorUtil.DeltaR(photon,twoprong), weight)
          if self.photon == 'HPID':
            if abs(photons[recophi.photonindex].scEta)<1.4442: photon_subdet = 'barrel'
            if abs(photons[recophi.photonindex].scEta)>1.566 and abs(photons[recophi.photonindex].scEta)<2.5: photon_subdet = 'endcap'
          if self.photon == 'CBL':
            if photons[recophi.photonindex].isScEtaEB: photon_subdet = 'barrel'
            if photons[recophi.photonindex].isScEtaEE: photon_subdet = 'endcap'
          if photon_subdet == 'barrel':
            self.recophi_b_pt.Fill(recophi.pt, weight)
            self.recophi_b_eta.Fill(recophi.eta, weight)
            self.recophi_b_phi.Fill(recophi.phi, weight)
            self.recophi_b_m.Fill(recophi.mass, weight)
            self.recophi_b_dphi.Fill(abs(ROOT.Math.VectorUtil.DeltaPhi(photon,twoprong)), weight)
            self.recophi_b_deta.Fill(abs(photon.Eta() - twoprong.Eta()), weight)
            self.recophi_b_dr.Fill(ROOT.Math.VectorUtil.DeltaR(photon,twoprong), weight)
          if photon_subdet == 'endcap':
            self.recophi_e_pt.Fill(recophi.pt, weight)
            self.recophi_e_eta.Fill(recophi.eta, weight)
            self.recophi_e_phi.Fill(recophi.phi, weight)
            self.recophi_e_m.Fill(recophi.mass, weight)
            self.recophi_e_dphi.Fill(abs(ROOT.Math.VectorUtil.DeltaPhi(photon,twoprong)), weight)
            self.recophi_e_deta.Fill(abs(photon.Eta() - twoprong.Eta()), weight)
            self.recophi_e_dr.Fill(ROOT.Math.VectorUtil.DeltaR(photon,twoprong), weight)
          if event.NJets == 0: self.recophi_dr_0jet.Fill(ROOT.Math.VectorUtil.DeltaR(photon,twoprong), weight)
          if event.NJets == 1: self.recophi_dr_1jet.Fill(ROOT.Math.VectorUtil.DeltaR(photon,twoprong), weight)
          if event.NJets >= 2: self.recophi_dr_2jet.Fill(ROOT.Math.VectorUtil.DeltaR(photon,twoprong), weight)
          self.photon_pt.Fill(photon.Pt(), weight)
          self.photon_eta.Fill(photon.Eta(), weight)
          self.photon_phi.Fill(photon.Phi(), weight)
          self.twoprong_pt.Fill(twoprong.Pt(), weight)
          self.twoprong_eta.Fill(twoprong.Eta(), weight)
          self.twoprong_phi.Fill(twoprong.Phi(), weight)
          self.twoprong_mass.Fill(twoprongs[recophi.twoprongindex].mass, weight)
          self.twoprong_masspi0.Fill(twoprongs[recophi.twoprongindex].massPi0, weight)
          self.twoprong_masseta.Fill(twoprongs[recophi.twoprongindex].massEta, weight)
          self.nphoton.Fill(event.nHighPtIdPhoton, weight)
          self.ntwoprong.Fill(ntwoprong, weight)
          self.njets.Fill(event.NJets, weight)
          self.ht.Fill(event.HT, weight)
          self.met.Fill(event.MET_pt, weight)
          self.met_phi.Fill(event.MET_phi, weight)
          if event.MET_pt > 30: self.met_phi_cut.Fill(event.MET_phi, weight)
          self.npv.Fill(event.PV_npvs, weight)
          if self.photon == 'HPID':
            if abs(photons[recophi.photonindex].scEta)<1.4442: self.twoprong_eta_barrel.Fill(twoprong.Eta(), weight)
            if abs(photons[recophi.photonindex].scEta)>1.566 and abs(photons[recophi.photonindex].scEta)<2.5: self.twoprong_eta_endcap.Fill(twoprong.Eta(), weight)
          if self.photon == 'CBL':
            if photons[recophi.photonindex].isScEtaEB: self.twoprong_eta_barrel.Fill(twoprong.Eta(), weight)
            if photons[recophi.photonindex].isScEtaEE: self.twoprong_eta_endcap.Fill(twoprong.Eta(), weight)

        return True
