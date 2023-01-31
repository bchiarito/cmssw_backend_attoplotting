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

        self.recophi_dr = ROOT.TH1F('recophi_dr', 'recophi_dr', 1000, 0, 10)
        self.addObject(self.recophi_dr)
        self.m_g2 = ROOT.TH1F('m_gr', 'm_gr', 100, 0, 1000)
        self.addObject(self.m_g2)
        self.photon_pt = ROOT.TH1F('photon_pt', 'photon_pt', 200, 0, 2000)
        self.addObject(self.photon_pt)

        self.hthat_lhe = ROOT.TH1F('MC_hthat_lhe', 'hthat_lhe', 100, 0, 1000)
        self.addObject(self.hthat_lhe)
        self.hthat_genPart = ROOT.TH1F('MC_hthat_genPart', 'hthat_genPart', 100, 0, 1000)
        self.addObject(self.hthat_genPart)

        self.cutflow_ = ROOT.TH1F('cutflow', 'cutflow', 100, 0, 100)
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

        self.cutflow.Fill(0)

        if event.RecoPhi_pass and recophi.photonLeg_pt>200:
          self.cutflow.Fill(1)
          photon = ROOT.TLorentzVector()
          twoprong = ROOT.TLorentzVector()
          photon.SetPtEtaPhiM(recophi.photonLeg_pt, recophi.photonLeg_eta, recophi.photonLeg_phi, recophi.photonLeg_mass)
          twoprong.SetPtEtaPhiM(recophi.twoprongLeg_pt, recophi.twoprongLeg_eta, recophi.twoprongLeg_phi, recophi.twoprongLeg_mass)
          self.recophi_dr.Fill(photon.DeltaR(twoprong), weight)
          self.m_g2.Fill(recophi.mass, weight)
          self.photon_pt.Fill(recophi.photonLeg_pt, weight)

        try:
          self.hthat_lhe.Fill(event.htHat_lhe, weight)
          self.hthat_genPart.Fill(event.htHat_genPart, weight)
        except RuntimeError:
          pass

        return True
