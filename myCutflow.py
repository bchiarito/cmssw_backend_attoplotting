#!/usr/bin/env python
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True


class MyCutflow(Module):
    def __init__(self, filename):
        self.writeHistFile = False
        self.filename = filename
        self.total = 0
        self.passing_og = 0
        self.passing_new = 0
    
    def endJob(self):
        print('hi, this is the cutflow')
        print('total '+str(self.total))
        print('pass og '+str(self.passing_og))
        print('pass new '+str(self.passing_new))
        with open(self.filename, 'w') as f:
          f.write('total '+str(self.total)+'\n')
          f.write('pass og '+str(self.passing_og)+'\n')
          f.write('pass new '+str(self.passing_new)+'\n')

    def beginJob(self, histFile=None, histDirName=None):
        pass

    def analyze(self, event):
        recophi = Object(event, "RecoPhi")

        self.total += 1
        if event.RecoPhi_pass and event.HighPtIdPhoton_pt[0]>200: self.passing_og += 1
        if event.RecoPhi_pass and recophi.photonLeg_pt>200: self.passing_new += 1

        return True
