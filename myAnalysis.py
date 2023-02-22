#!/usr/bin/env python
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

class MyAnalysis(Module):
    def __init__(self, lumi=1.0, dict_xs=None, dict_ngen=None):
        self.writeHistFile = True
        self.lumi = lumi
        self.dict_xs = dict_xs
        self.dict_ngen = dict_ngen

    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)

        self.njets = ROOT.TH1F('njets', 'njets', 20, 0, 20)
        self.addObject(self.njets)
        self.ht = ROOT.TH1F('ht', 'ht', 250, 0, 5000)
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
        self.recophi_m = ROOT.TH1F('reophi_m', 'recophi_m', 300, 0, 6000)
        self.addObject(self.recophi_m)
        self.recophi_pt = ROOT.TH1F('reophi_pt', 'recophi_pt', 200, 0, 2000)
        self.addObject(self.recophi_m)
        self.recophi_eta = ROOT.TH1F('reophi_eta', 'recophi_eta', 100, -5, 5)
        self.addObject(self.recophi_m)
        self.recophi_phi = ROOT.TH1F('reophi_phi', 'recophi_phi', 70, -3.5, 3.5)
        self.addObject(self.recophi_m)

        self.photon_pt = ROOT.TH1F('photon_pt', 'photon_pt', 100, 0, 1000)
        self.addObject(self.photon_pt)
        self.photon_eta = ROOT.TH1F('photon_eta', 'photon_eta', 100, -5, 5)
        self.addObject(self.photon_eta)
        self.photon_phi = ROOT.TH1F('photon_phi', 'photon_phi', 70, -3.5, 3.5)
        self.addObject(self.photon_phi)
        self.nphoton = ROOT.TH1F('nphoton', 'nphoton', 10, 0, 10)
        self.addObject(self.nphoton)

        self.twoprong_pt = ROOT.TH1F('twoprong_pt', 'twoprong_pt', 100, 0, 1000)
        self.addObject(self.twoprong_pt)
        self.twoprong_eta = ROOT.TH1F('twoprong_eta', 'twoprong_eta', 100, -5, 5)
        self.addObject(self.twoprong_eta)
        self.twoprong_phi = ROOT.TH1F('twoprong_phi', 'twoprong_phi', 70, -3.5, 3.5)
        self.addObject(self.twoprong_phi)
        self.twoprong_mass = ROOT.TH1F('twoprong_mass', 'twoprong_mass', 120, 0, 12)
        self.addObject(self.twoprong_mass)
        self.twoprong_masspi0 = ROOT.TH1F('twoprong_masspi0', 'twoprong_masspi0', 120, 0, 12)
        self.addObject(self.twoprong_masspi0)
        self.twoprong_masseta = ROOT.TH1F('twoprong_masseta', 'twoprong_masseta', 120, 0, 12)
        self.addObject(self.twoprong_masseta)
        self.ntwoprong = ROOT.TH1F('ntwoprong', 'ntwoprong', 10, 0, 10)
        self.addObject(self.ntwoprong)

        self.hthat_lhe = ROOT.TH1F('MC_hthat_lhe', 'hthat_lhe', 100, 0, 1000)
        self.addObject(self.hthat_lhe)
        self.hthat_genPart = ROOT.TH1F('MC_hthat_genPart', 'hthat_genPart', 100, 0, 1000)
        self.addObject(self.hthat_genPart)

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

        self.cutflow.Fill(0)
        if event.RecoPhi_pass:
          self.cutflow.Fill(1)
        if event.RecoPhi_pass and recophi.photonLeg_pt>200:
          self.cutflow.Fill(2)
          photon = ROOT.TLorentzVector()
          twoprong = ROOT.TLorentzVector()
          photon.SetPtEtaPhiM(recophi.photonLeg_pt, recophi.photonLeg_eta, recophi.photonLeg_phi, recophi.photonLeg_mass)
          twoprong.SetPtEtaPhiM(recophi.twoprongLeg_pt, recophi.twoprongLeg_eta, recophi.twoprongLeg_phi, recophi.twoprongLeg_mass)
          self.recophi_dr.Fill(photon.DeltaR(twoprong), weight)
          self.recophi_dphi.Fill(photon.DeltaPhi(twoprong), weight)
          self.recophi_deta.Fill(abs(photon.Eta() - twoprong.Eta()), weight)
          if event.NJets == 0: self.recophi_dr_0jet.Fill(photon.DeltaR(twoprong), weight)
          if event.NJets == 1: self.recophi_dr_1jet.Fill(photon.DeltaR(twoprong), weight)
          if event.NJets >= 2: self.recophi_dr_2jet.Fill(photon.DeltaR(twoprong), weight)
          self.recophi_m.Fill(recophi.mass, weight)
          self.recophi_pt.Fill(recophi.pt, weight)
          self.recophi_eta.Fill(recophi.eta, weight)
          self.recophi_phi.Fill(recophi.phi, weight)
          self.photon_pt.Fill(recophi.photonLeg_pt, weight)
          self.photon_eta.Fill(recophi.photonLeg_eta, weight)
          self.photon_phi.Fill(recophi.photonLeg_phi, weight)
          self.twoprong_pt.Fill(recophi.twoprongLeg_pt, weight)
          self.twoprong_eta.Fill(recophi.twoprongLeg_eta, weight)
          self.twoprong_phi.Fill(recophi.twoprongLeg_phi, weight)
          for twoprong in twoprongs:
            if twoprong.isTight:
              self.twoprong_mass.Fill(twoprong.mass, weight)
              self.twoprong_masspi0.Fill(twoprong.massPi0, weight)
              self.twoprong_masseta.Fill(twoprong.massEta, weight)
              break

          self.nphoton.Fill(event.nHighPtIdPhoton, weight)
          ntwoprong = 0
          for twoprong in twoprongs:
            if twoprong.isTight: ntwoprong += 1
          self.ntwoprong.Fill(ntwoprong, weight)

          self.njets.Fill(event.NJets, weight)
          self.ht.Fill(event.HT, weight)
          self.met.Fill(event.MET_pt, weight)
          self.met_phi.Fill(event.MET_phi, weight)
          if event.MET_pt > 30: self.met_phi_cut.Fill(event.MET_phi, weight)
          self.npv.Fill(event.PV_npvs, weight)

        try:
          self.hthat_lhe.Fill(event.htHat_lhe, weight)
          self.hthat_genPart.Fill(event.htHat_genPart, weight)
        except RuntimeError:
          pass

        return True
