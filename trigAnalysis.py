#!/usr/bin/env python
#!/usr/bin/env python
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
import array
ROOT.PyConfig.IgnoreCommandLineOptions = True


class TriggerAnalysis(Module):
    def __init__(self, lumi=1.0, dict_xs=None, dict_ngen=None):
        self.writeHistFile = True
        self.year = "2018"
        self.lumi = lumi
        self.dict_xs = dict_xs
        self.dict_ngen = dict_ngen

    def beginJob(self, histFile=None, histDirName=None):
        Module.beginJob(self, histFile, histDirName)
        self.photon_pt_num = ROOT.TH1F('photon_pt_num', '; Photon p_{T}', 160, 0, 1600)
        self.photon_pt_denom = ROOT.TH1F('photon_pt_denom', '; Photon p_{T}', 160, 0, 1600)
        self.photon_eta_num = ROOT.TH1F('photon_eta_num', '; Photon_Eta', 15, -1.5, 1.5)
        self.photon_eta_denom = ROOT.TH1F('photon_eta_denom', '; Photon_Eta', 15, -1.5, 1.5)
        self.photon_phi_num = ROOT.TH1F('photon_phi_num', '; Photon_phi', 15, -3.5, 3.5)
        self.photon_phi_denom = ROOT.TH1F('photon_phi_denom', '; Photon_phi', 15, -3.5, 3.5)

        pt_bin = array.array('f', [220, 260, 300, 350, 400, 450, 500, 600, 800, 2000])
        eta_bin = array.array('f', [0, 0.8, 1.44442])
        self.photon_pt_eta_2d_num = ROOT.TH2F("photon_pt_eta_2d_num", "; Photon p_{T}; Photon Eta", len(pt_bin)-1, pt_bin, len(eta_bin)-1, eta_bin)
        self.photon_pt_eta_2d_denom = ROOT.TH2F("photon_pt_eta_2d_denom", "; Photon p_{T}; Photon Eta", len(pt_bin)-1, pt_bin, len(eta_bin)-1, eta_bin)
         
        self.addObject(self.photon_pt_num)
        self.addObject(self.photon_pt_denom)
        self.addObject(self.photon_eta_num)
        self.addObject(self.photon_eta_denom)
        self.addObject(self.photon_phi_num)
        self.addObject(self.photon_phi_denom)
        self.addObject(self.photon_pt_eta_2d_num)
        self.addObject(self.photon_pt_eta_2d_denom)

    
    def mygetattr(self, my_obj, my_branch, default_bool):
        try:
            return getattr(my_obj, my_branch)
        except RuntimeError:
            return default_bool


    def analyze(self, event):
        twoprongs = Collection(event, "TwoProng")
        flags     = Object(event, "Flag")
        photons   = Collection(event, "Photon")
        passtrig = event.HLT_Photon200
        
        # Filters
        pass_filter = (
            self.mygetattr(flags, 'goodVertices', True)
            and self.mygetattr(flags, 'HBHENoiseFilter', True)
            and self.mygetattr(flags, 'HBHENoiseIsoFilter', True)
            and self.mygetattr(flags, 'EcalDeadCellTriggerPrimitiveFilter', True)
            and self.mygetattr(flags, 'BadPFMuonFilter', True)
            and self.mygetattr(flags, 'BadChargedCandidateFilter', True)
            and self.mygetattr(flags, 'ecalBadCalibFilter', True)
            and self.mygetattr(flags, 'globalSuperTightHalo2016Filter', True)
            and self.mygetattr(flags, 'eeBadScFilter', True)
            )
        if not pass_filter: return True 

        ps = []
        for p in photons:
            #if p.cutBased >= 1 and p.electronVeto == 1 and abs(p.scEta) < 1.4442:
            if p.cutBased >= 1 and p.electronVeto == 1 and p.isScEtaEB:
                ps.append(p)
                break
        
        if len(ps) == 0: return True
        
        # Get MC weights
        if self.dict_xs and self.dict_ngen:
            dataset_id = event.dataset_id
            xs = self.dict_xs[dataset_id]
            Ngen = self.dict_ngen[dataset_id]
            weight = xs * self.lumi / Ngen
        else:
            weight = 1.0

        self.photon_pt_denom.Fill(ps[0].pt, weight)
        if ps[0].pt > 220:
            self.photon_eta_denom.Fill(ps[0].eta, weight)
            self.photon_phi_denom.Fill(ps[0].phi, weight)
            self.photon_pt_eta_2d_denom.Fill(ps[0].pt, ps[0].eta, weight)
        if passtrig:
            self.photon_pt_num.Fill(ps[0].pt, weight)
            if ps[0].pt > 220:
                self.photon_eta_num.Fill(ps[0].eta, weight)
                self.photon_phi_num.Fill(ps[0].phi, weight)
                self.photon_pt_eta_2d_num.Fill(ps[0].pt, ps[0].eta, weight)

        return True
