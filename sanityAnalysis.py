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

class SanityAnalysis(Module):
    def __init__(self, datamc, lumi=1.0, dict_xs=None, dict_ngen=None, deta=False, photon='HPID'):
        self.writeHistFile = True
        self.lumi = lumi
        self.dict_xs = dict_xs
        self.dict_ngen = dict_ngen
        self.datamc = datamc
        self.deta = deta
        self.photon = photon

    def book_histo(self, name, bins, low, high, title=None):
        if not title: title = name
        exec("self."+name+" = ROOT.TH1F(name, title, bins, low, high)")
        exec("self.addObject(self."+name+")")

    def book_manyEBEB(self, name, var, bins, low, high):
        self.book_histo(name+'_'+var, bins, low, high)
        self.book_histo(name+'_B_B_'+var, bins, low, high)
        self.book_histo(name+'_B_E_'+var, bins, low, high)
        self.book_histo(name+'_E_B_'+var, bins, low, high)
        self.book_histo(name+'_E_E_'+var, bins, low, high)

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

        self.book_histo('SIGNAL_decaymode', 21, -1, 20, title='decaymode')
        self.book_histo('SIGNAL_prongs', 5, -1, 4, title='prongs')
        self.book_histo('SIGNAL_tag_num_pt', 150, 0, 1500, title='tag_num_pt')
        self.book_histo('SIGNAL_tag_den_pt', 150, 0, 1500, title='tag_den_pt')
        self.book_histo('SIGNAL_tag_eff_pt', 150, 0, 1500, title='tag_eff_pt')

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
        genomegas = Collection(event, "GenOmega")
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
          the_photon = get_vec(photons[recophi.photonindex])
          the_twoprong = get_vec(twoprongs[recophi.twoprongindex])
          deta = abs(the_photon.Eta() - the_twoprong.Eta())

        # deta cut
        if region == 1:
          if self.deta and deta < 1.5: pass_deta = True
          elif self.deta and deta > 1.5: pass_deta = False
          else: pass_deta = True
        
        # mc hthat
        self.cutflow.Fill(0)
        try:
          #if len(photons)>=1:
          if True:
            self.hthat_gjets.Fill(event.htHat_lhe, weight)
            self.hthat_qcd.Fill(event.htHat_lhe, weight)
        except RuntimeError:
          pass

        if region == 1:
          self.cutflow.Fill(1)

        if region == 1 and the_photon.Pt() > 220 and abs(photons[recophi.photonindex].scEta) < 1.4442 and pass_deta:
          self.cutflow.Fill(2)

        if region == 1 and the_photon.Pt() > 220 and abs(photons[recophi.photonindex].scEta) < 1.4442 and pass_trigger and pass_deta:
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
          #elif abs(twoprongs[recophi.twoprongindex].eta)>1.566 and abs(twoprongs[recophi.twoprongindex].eta)<2.5: twoprong_subdet = 'endcap'
          else: twoprong_subdet = 'endcap'
          self.recophi_pt.Fill(recophi.pt, weight)
          self.recophi_eta.Fill(recophi.eta, weight)
          self.recophi_phi.Fill(recophi.phi, weight)
          self.recophi_m.Fill(recophi.mass, weight)
          self.recophi_dphi.Fill(abs(ROOT.Math.VectorUtil.DeltaPhi(the_photon,the_twoprong)), weight)
          self.recophi_deta.Fill(abs(the_photon.Eta() - the_twoprong.Eta()), weight)
          self.recophi_dr.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
          if event.NJets == 0: self.recophi_dr_0jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
          if event.NJets == 1: self.recophi_dr_1jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
          if event.NJets >= 2: self.recophi_dr_2jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
          if photon_subdet == 'barrel' and twoprong_subdet == 'barrel':
            self.recophi_B_B_pt.Fill(recophi.pt, weight)
            self.recophi_B_B_eta.Fill(recophi.eta, weight)
            self.recophi_B_B_phi.Fill(recophi.phi, weight)
            self.recophi_B_B_m.Fill(recophi.mass, weight)
            self.recophi_B_B_dphi.Fill(abs(ROOT.Math.VectorUtil.DeltaPhi(the_photon,the_twoprong)), weight)
            self.recophi_B_B_deta.Fill(abs(the_photon.Eta() - the_twoprong.Eta()), weight)
            self.recophi_B_B_dr.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
            if event.NJets == 0: self.recophi_B_B_dr_0jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
            if event.NJets == 1: self.recophi_B_B_dr_1jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
            if event.NJets >= 2: self.recophi_B_B_dr_2jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
          if photon_subdet == 'barrel' and twoprong_subdet == 'endcap':
            self.recophi_B_E_pt.Fill(recophi.pt, weight)
            self.recophi_B_E_eta.Fill(recophi.eta, weight)
            self.recophi_B_E_phi.Fill(recophi.phi, weight)
            self.recophi_B_E_m.Fill(recophi.mass, weight)
            self.recophi_B_E_dphi.Fill(abs(ROOT.Math.VectorUtil.DeltaPhi(the_photon,the_twoprong)), weight)
            self.recophi_B_E_deta.Fill(abs(the_photon.Eta() - the_twoprong.Eta()), weight)
            self.recophi_B_E_dr.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
            if event.NJets == 0: self.recophi_B_E_dr_0jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
            if event.NJets == 1: self.recophi_B_E_dr_1jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
            if event.NJets >= 2: self.recophi_B_E_dr_2jet.Fill(ROOT.Math.VectorUtil.DeltaR(the_photon,the_twoprong), weight)
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
          self.photon_pt.Fill(the_photon.Pt(), weight)
          self.photon_eta.Fill(the_photon.Eta(), weight)
          self.photon_phi.Fill(the_photon.Phi(), weight)
          self.twoprong_pt.Fill(the_twoprong.Pt(), weight)
          self.twoprong_eta.Fill(the_twoprong.Eta(), weight)
          self.twoprong_phi.Fill(the_twoprong.Phi(), weight)
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
            if abs(photons[recophi.photonindex].scEta)<1.4442: self.twoprong_eta_barrel.Fill(the_twoprong.Eta(), weight)
            if abs(photons[recophi.photonindex].scEta)>1.566 and abs(photons[recophi.photonindex].scEta)<2.5: self.twoprong_eta_endcap.Fill(twoprong.Eta(), weight)
          if self.photon == 'CBL':
            if photons[recophi.photonindex].isScEtaEB: self.twoprong_eta_barrel.Fill(the_twoprong.Eta(), weight)
            if photons[recophi.photonindex].isScEtaEE: self.twoprong_eta_endcap.Fill(the_twoprong.Eta(), weight)

        # signal
        try:
          for i in range(event.nGenOmega):
            self.SIGNAL_decaymode.Fill(event.GenOmega_decaymode[i])
            self.SIGNAL_prongs.Fill(event.GenOmega_prongs[i])
          
          for genomega in genomegas:
            genvec = get_vec(genomega)
            tagged = False
            for twoprong in twoprongs:
              twoprongvec = get_vec(twoprong)
              ROOT.Math.VectorUtil.DeltaPhi(genvec, twoprongvec) # line does nothing but must be run
              if ROOT.Math.VectorUtil.DeltaR(genvec, twoprongvec) < 0.1 and ROOT.Math.VectorUtil.DeltaPhi(genvec, twoprongvec) < 0.1:
                tagged = True
                break
            self.SIGNAL_tag_den_pt.Fill(genvec.Pt())
            if tagged: self.SIGNAL_tag_num_pt.Fill(genvec.Pt())
            if tagged: self.SIGNAL_tag_eff_pt.Fill(genvec.Pt())
          #self.SIGNAL_tag_eff_pt = (self.SIGNAL_tag_num_pt.Clone(self.SIGNAL_tag_eff_pt.GetName()))
          self.SIGNAL_tag_eff_pt.Divide(self.SIGNAL_tag_den_pt)
          #self.SIGNAL_tag_eff_pt.Print('All')
          
        except RuntimeError:
          pass

        return True