#!/usr/bin/env python
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True


class HthatAnalysis(Module):
    def __init__(self, lumi=1.0, dict_xs=None, dict_ngen=None):
        self.writeHistFile = True
        self.lumi = lumi
        self.dict_xs = dict_xs
        self.dict_ngen = dict_ngen

    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)

        self.hthat_lhe = ROOT.TH1F('hthat_lhe', 'hthat_lhe', 100, 0, 1000)
        self.addObject(self.hthat_lhe)
        self.hthat_genPart = ROOT.TH1F('hthat_genPart', 'hthat_genPart', 100, 0, 1000)
        self.addObject(self.hthat_genPart)

    def analyze(self, event):
        if self.dict_xs and self.dict_ngen:
          flag = event.flag
          dataset_id = event.dataset_id
          xs = self.dict_xs[dataset_id]
          Ngen = self.dict_ngen[dataset_id]
          weight = xs * self.lumi / Ngen
        else:
          weight = 1.0

        self.hthat_lhe.Fill(event.htHat_lhe, weight)
        self.hthat_genPart.Fill(event.htHat_genPart, weight)

        return True
