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
    def __init__(self, datamc, lumi=1.0, dict_xs=None, dict_ngen=None):
        self.writeHistFile = True
        self.lumi = lumi
        self.dict_xs = dict_xs
        self.dict_ngen = dict_ngen
        self.datamc = datamc

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
        self.recophi_m = ROOT.TH1F('reophi_m', 'recophi_m', 200, 0, 6000)
        self.addObject(self.recophi_m)
        self.recophi_pt = ROOT.TH1F('reophi_pt', 'recophi_pt', 200, 0, 2000)
        self.addObject(self.recophi_pt)
        self.recophi_eta = ROOT.TH1F('reophi_eta', 'recophi_eta', 200, -10, 10)
        self.addObject(self.recophi_eta)
        self.recophi_phi = ROOT.TH1F('reophi_phi', 'recophi_phi', 70, -3.5, 3.5)
        self.addObject(self.recophi_phi)

        self.photon_pt = ROOT.TH1F('photon_pt', 'photon_pt', 150, 0, 1500)
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

        self.hthat_gjets = ROOT.TH1F('GJETS_hthat_lhe', 'hthat', 100, 0, 1000)
        self.addObject(self.hthat_gjets)
        self.hthat_qcd = ROOT.TH1F('QCD_hthat_lhe', 'hthat', 300, 0, 3000)
        self.addObject(self.hthat_qcd)

        self.cutflow = ROOT.TH1F('cutflow', 'cutflow', 10, 0, 10)
        self.addObject(self.cutflow)

    def analyze(self, event):
        if self.dict_xs and self.dict_ngen:
          dataset_id = event.dataset_id
          xs = self.dict_xs[dataset_id]
          Ngen = self.dict_ngen[dataset_id]
          weight = xs * self.lumi / Ngen
        else:
          weight = 1.0

        recophi = Object(event, "RecoPhi")
        twoprongs = Collection(event, "TwoProng")
        photons = Collection(event, "HighPtIdPhoton")
        pass_trigger = event.HLT_Photon200
        #if self.datamc == 'mc' or self.datamc == 'sigRes' or self.datamc == 'sigNonRes': pass_trigger = True
        ntwoprong = 0
        for twoprong in twoprongs:
          try:
            tight = twoprong.isTight
          except RuntimeError:
            tight = True
          if tight: ntwoprong += 1
        if event.Region == 1:
          photon = get_vec(photons[recophi.photonindex])
          twoprong = get_vec(twoprongs[recophi.twoprongindex])

        if True:
          self.cutflow.Fill(0)
          try:
            self.hthat_gjets.Fill(event.htHat_lhe, weight)
            self.hthat_qcd.Fill(event.htHat_lhe, weight)
          except RuntimeError:
            pass

        if event.Region == 1:
          self.cutflow.Fill(1)

        if event.Region == 1 and photon.Pt() > 280:
          self.cutflow.Fill(2)

        if event.Region == 1 and photon.Pt() > 280 and pass_trigger:
          self.cutflow.Fill(3)
          self.recophi_pt.Fill(recophi.pt, weight)
          self.recophi_eta.Fill(recophi.eta, weight)
          self.recophi_phi.Fill(recophi.phi, weight)
          self.recophi_m.Fill(recophi.mass, weight)
          self.recophi_dphi.Fill(abs(ROOT.Math.VectorUtil.DeltaPhi(photon,twoprong)), weight)
          self.recophi_deta.Fill(abs(photon.Eta() - twoprong.Eta()), weight)
          self.recophi_dr.Fill(ROOT.Math.VectorUtil.DeltaR(photon,twoprong), weight)
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

        return True
