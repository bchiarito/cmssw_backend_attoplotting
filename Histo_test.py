from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
import numpy as np
import math
import ctypes
import json
ROOT.PyConfig.IgnoreCommandLineOptions = True


class HistProd(Module):
    def __init__(self,datamc):
        self.writeHistFile = True
        self.kSpreadMC = None
    	self.kSmearMC = None
    	self.kScaleDT = None
    	self.roccin = None
	self.datamc = datamc
        self.histograms = {}
        self.string_list = [
	('Dataset_Scale_Factor','GSF_dataset','sf','NA'),
        ('Central_Scale_Dist','GSF_central','sf','NA'),
        ('Upper_Scale_Dist','GSF_plus','sf','NA'),
        ('Lower_Scale_Dist','GSF_minus','sf','NA'),
        ('Rochester_Scale_Dist','GSF_rochester','sfroc','NA'),
        ('ID_ISO_Trig_Scale_Dist','GSF_id_iso_trig','sfid','NA'),
        ('HEM_Scale_Dist','GSF_HEM','sfhem','NA'),
        ('Pileup_Central_Scale_Dist','GSF_pileup_central','sfpile','NA'),
        ('Pileup_Upper_Scale_Dist','GSF_pileup_plus','sfpile','NA'),
        ('Pileup_Lower_Scale_Dist','GSF_pileup_minus','sfpile','NA'),
        ('GnJets_corrected','GnJets_corrected','num','Jet'),
#	('nJets','nJets','num','Jet'),
	('GJet_mass','GJet_mass','mass','Jet'),
        ('GJet_pt','GJet_pt','pt','Jet'),
        ('GJet_eta','GJet_eta','etaphi','Jet'),
        ('GJet_phi','GJet_phi','etaphi','Jet'),
        ('GMet','GMET','met','met'),
        ('GMuon_pt','GMuon_pt','pt','Muon'),
        ('GnMuon','GnMuon','num','Muon'),
        ('GMuon_eta','GMuon_eta','etaphi','Muon'),
        ('GMuon_phi','GMuon_phi','etaphi','Muon'),
        ('GMuon_mass','GMuon_mass','massM','Muon'),
	('GMuon_genpt','GMuon_genpt','pt','Muon'),
        ('GMuon_charge','GMuon_charge','charge','Muon'),
	('GMuon_ntrackerlayers','GMuon_ntrackerlayers','num','Muon'),
        ('GTau_pt','GTau_pt','pt','Tau'),
        ('GnTau','GnTau','num','Tau'),
        ('GTau_eta','GTau_eta','etaphi','Tau'),
        ('GTau_phi','GTau_phi','etaphi','Tau'),
        ('GTau_mass','GTau_mass','massT','Tau'),
        ('GTau_charge','GTau_charge','charge','Tau'),
        ('GTwoProng_pt','GTwoProng_pt','pt','TwoProng'),
        ('GTPnTracks','GTPnTracks','num','TwoProng'),
        ('GTwoProng_eta','GTwoProng_eta','etaphi','TwoProng'),
        ('GTwoProng_phi','GTwoProng_phi','etaphi','TwoProng'),
        ('GTwoProng_massl','GTwoProng_massl','massTP','TwoProng'),
        ('GTwoProng_massEta','GTwoProng_massEta','massTP','TwoProng'),
        ('GTwoProng_massPi0','GTwoProng_massPi0','massTP','TwoProng'),
        ('GChextra_charge','GChextra_charge','charge','TwoProng'),
        ('GChpospt','GChpospt','pt','TwoProng'),
        ('GChnegpt','GChnegpt','pt','TwoProng'),
        ('Z_massTau','Z_massTau','massZ','Tau'),
        ('Z_ptTau','Z_ptTau','pt','Tau'),
        ('Z_etaTau','Z_etaTau','etaphi','Tau'),
        ('Z_phiTau','Z_phiTau','etaphi','Tau'),
        ('Z_massTP','Z_massTP','massZ','TP'),
        ('Z_ptTP','Z_ptTP','pt','TP'),
        ('Z_etaTP','Z_etaTP','etaphi','TP'),
        ('Z_phiTP','Z_phiTP','etaphi','TP'),
	('TransMass','GTransMass','mass','Muon'),
	('GNPV','GNPV','num','Tau')
        ]

        self.syst_list = [
        ('_c',0),('_p',1),('_m',2)]

        self.types = [
        {'type': 'sf','bins':1000 ,'min': 0, 'max':100},
        {'type': 'met','bins':100 ,'min': 0, 'max':800},
        {'type': 'sfpile','bins':1000 ,'min': 0, 'max':3},
        {'type': 'sfroc','bins':1000 ,'min': 0.98, 'max':1.02},
        {'type': 'sfid','bins':1000 ,'min': 0.94, 'max':1.05},
        {'type': 'sfhem','bins':50 ,'min': 0, 'max':1.2},
        {'type': 'pt','bins':100 ,'min': 0, 'max':400},
        {'type': 'massZ','bins':100 ,'min': 0, 'max':400},
        {'type': 'massM','bins':100 ,'min': 0, 'max':0.5},\
        {'type': 'massT','bins':100 ,'min': 0, 'max':2.5},
        {'type': 'massTP','bins':100 ,'min': 0, 'max':20},
        {'type': 'etaphi','bins':100 ,'min': -4, 'max':4},
        {'type': 'met','bins':100 ,'min': 0, 'max':800},
        {'type': 'tracks','bins':100 ,'min': 0, 'max':10},
        {'type': 'delr','bins':100 ,'min': 0, 'max':0.6},
        {'type': 'charge','bins':100 ,'min': -1.2, 'max':1.2},
        {'type': 'num','bins':10 ,'min': 0, 'max':10}
        ]
    	rc = ctypes.CDLL('/cms/smd376/CMSSW_11_1_0/src/RocWrapper.so')
	print(rc)
	roccorfile = "/cms/smd376/CMSSW_11_1_0/src/PhysicsTools/NanoAODTools/python/postprocessing/data/roccor.Run2.v3/RoccoR2018.txt"
	RoccoR_instance = rc.RoccoRinit
	RoccoR_instance.restype = ctypes.c_void_p
	Roccin = RoccoR_instance(roccorfile)
	self.roccin = Roccin
	print("Instance of RoccoR created")
	rc.kScaleDT.argtypes = [ctypes.c_void_p,ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_int,ctypes.c_int]

	rc.kScaleDT.restype = ctypes.c_double

	rc.kSpreadMC.argtypes = [ctypes.c_void_p,ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_int,ctypes.c_int]

	rc.kSpreadMC.restype = ctypes.c_double

	rc.kSmearMC.argtypes = [ctypes.c_void_p,ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double,ctypes.c_int,ctypes.c_double,ctypes.c_int,ctypes.c_int]

	rc.kSmearMC.restype = ctypes.c_double
	self.kSpreadMC = rc.kSpreadMC
	self.kSmearMC = rc.kSmearMC
	self.kScaleDT = rc.kScaleDT
    def set_type_vals(self,type):
	bins, min,max = 100,0,400
	for item in self.types:
		if item['type'] == type:
       			bins = item['bins']
			min = item['min']
			max = item['max']
        		break
	return bins,min,max


    def beginJob(self, histFile=None, histDirName=None):
	Module.beginJob(self, histFile, histDirName)



        for iter in self.string_list:
            name=iter[0]
            branch=iter[1]

            type=iter[2]
            cat=iter[3]
            bins,min,max = self.set_type_vals(type)
            if cat == "NA":
                self.histograms[name] = ROOT.TH1D(name, name, bins, min, max)
                self.addObject(self.histograms[name])
            if cat != "NA" :
                name_ss_tp = name +"_ss_tp"
                name_ss_tau = name +"_ss_tau"
                name_os_tp = name +"_os_tp"
                name_os_tau = name +"_os_tau"
                self.histograms[name_ss_tp] = ROOT.TH1D(name_ss_tp, name_ss_tp, bins, min, max)
                self.addObject(self.histograms[name_ss_tp])
                self.histograms[name_ss_tau] = ROOT.TH1D(name_ss_tau, name_ss_tau, bins, min, max)
                self.addObject(self.histograms[name_ss_tau])
                self.histograms[name_os_tp] = ROOT.TH1D(name_os_tp, name_os_tp, bins, min, max)
                self.addObject(self.histograms[name_os_tp])
                self.histograms[name_os_tau] = ROOT.TH1D(name_os_tau, name_os_tau, bins, min, max)
            	self.addObject(self.histograms[name_os_tau])
            for iter2 in self.syst_list:
                systname = iter2[0]
                num = iter2[1]
                if cat != "NA":
                	name_ss_tp_syst = name +"_ss_tp" + systname
                        name_ss_tau_syst = name +"_ss_tau" + systname
                        name_os_tp_syst = name +"_os_tp" + systname
                        name_os_tau_syst = name +"_os_tau" + systname
                        self.histograms[name_ss_tp_syst] = ROOT.TH1D(name_ss_tp_syst, name_ss_tp_syst, bins, min, max)
                        self.addObject(self.histograms[name_ss_tp_syst])
                        self.histograms[name_ss_tau_syst] = ROOT.TH1D(name_ss_tau_syst, name_ss_tau_syst, bins, min, max)
                        self.addObject(self.histograms[name_ss_tau_syst])
                        self.histograms[name_os_tp_syst] = ROOT.TH1D(name_os_tp_syst, name_os_tp_syst, bins, min, max)
                        self.addObject(self.histograms[name_os_tp_syst])
                        self.histograms[name_os_tau_syst] = ROOT.TH1D(name_os_tau_syst, name_os_tau_syst, bins, min, max)
                        self.addObject(self.histograms[name_os_tau_syst])

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        dataset=1
        weights = [dataset*event.GSF_pileup_central*event.GSF_id_iso_trig,
        dataset*event.GSF_pileup_plus*event.GSF_id_iso_trig,
        dataset*event.GSF_pileup_minus*event.GSF_id_iso_trig]
        TwoProng_charge=0
        if event.GTPnTracks != 3:
        	if event.GChnegpt >event.GChpospt:
                    TwoProng_charge= -1
        	else:
           	        TwoProng_charge=1
        if event.GTPnTracks == 3:
            TwoProng_charge= np.signbit(event.GChextra_charge)
        for iter in self.string_list:
            name=iter[0]
            branch=iter[1]
            type=iter[2]
            cat=iter[3]
            name_ss_tp = name +"_ss_tp"
            name_ss_tau = name +"_ss_tau"
            name_os_tp = name +"_os_tp"
            name_os_tau = name +"_os_tau"
            name_ss = name +"_ss"
            if cat == "NA": self.histograms[name].Fill(getattr(event,branch))
            if cat != "NA":

                if  np.signbit(event.GMuon_charge) == np.signbit(event.GTau_charge) and event.GTau_pt >0 and event.GPZeta_Tau: self.histograms[name_ss_tau].Fill(getattr(event,branch))
                if  np.signbit(event.GMuon_charge) == np.signbit(TwoProng_charge) and event.GMuon_charge != 0 and TwoProng_charge !=0 and event.GTwoProng_pt >0 and event.GTPnTracks == 3 and event.GPZeta_TP: self.histograms[name_ss_tp].Fill(getattr(event,branch))
                if  np.signbit(event.GMuon_charge) != np.signbit(event.GTau_charge) and event.GTau_pt >0 and event.GPZeta_Tau: self.histograms[name_os_tau].Fill(getattr(event,branch))
                if  np.signbit(event.GMuon_charge) != np.signbit(TwoProng_charge) and event.GMuon_charge != 0 and TwoProng_charge !=0 and event.GTwoProng_pt >0 and event.GTPnTracks == 3 and event.GPZeta_TP: self.histograms[name_os_tp].Fill(getattr(event,branch))

            for iter2 in self.syst_list:
                systname = iter2[0]
                num = iter2[1]
		name_ss = name +"_ss"
                name_syst = name + systname
                name_syst_ss = name_ss + systname
                name_ss_tp_syst = name +"_ss_tp" + systname
                name_ss_tau_syst = name +"_ss_tau" + systname
                name_os_tp_syst = name +"_os_tp" + systname
                name_os_tau_syst = name +"_os_tau" + systname
                if cat != "NA":
                    if  np.signbit(event.GMuon_charge) == np.signbit(event.GTau_charge) and event.GTau_pt >0 and event.GPZeta_Tau: self.histograms[name_ss_tau_syst].Fill(getattr(event,branch),weights[num])
                    if  np.signbit(event.GMuon_charge) == np.signbit(TwoProng_charge) and event.GMuon_charge != 0 and TwoProng_charge !=0 and event.GTwoProng_pt >0 and event.GTPnTracks == 3 and event.GPZeta_TP: self.histograms[name_ss_tp_syst].Fill(getattr(event,branch),weights[num])
                    if  np.signbit(event.GMuon_charge) != np.signbit(event.GTau_charge) and event.GTau_pt >0 and event.GPZeta_Tau: self.histograms[name_os_tau_syst].Fill(getattr(event,branch),weights[num])
                    if  np.signbit(event.GMuon_charge) != np.signbit(TwoProng_charge) and event.GMuon_charge != 0 and TwoProng_charge !=0 and event.GTwoProng_pt >0 and event.GTPnTracks == 3 and event.GPZeta_TP: self.histograms[name_os_tp_syst].Fill(getattr(event,branch),weights[num])













        return True



        #This is how the rochester corrections are called in the Atto step. One supplies the charge, pt, eta, phi, and then either gen pt, tracker layer number, or nothing else for data.
        #if DMC==1:
            #if GoodMuon_genpt != -1 : sf_roc= self.kSpreadMC(self.roccin, int(GoodMuon_charge), GoodMuon_pt, GoodMuon_eta ,GoodMuon_phi,GoodMuon_genpt, 0,0 )
            #else: sf_roc= self.kSmearMC(self.roccin, int(GoodMuon_charge), GoodMuon_pt, GoodMuon_eta ,GoodMuon_phi ,int(GoodMuon_ntrackerlayer), nrandom, 0,0 )
        #else: sf_roc= self.kScaleDT(self.roccin, int(GoodMuon_charge), GoodMuon_pt, GoodMuon_eta, GoodMuon_phi, 0, 0)


        def PZeta_Cut(self, muon_pt,muon_eta,muon_phi,tau_pt,tau_eta,tau_phi, met,mcorrections):

            muvec= ROOT.TVector3()
            muvec.SetPtEtaPhi(muon_pt*mcorrections,muon_eta,0)
            tauvec= ROOT.TVector3()
            tauvec.SetPtEtaPhi(tau_pt,tau_eta,0)
            etvec= ROOT.TVector3()
            etvec.SetPtEtaPhi(met_pt,0,met_phi)
            zetaU=ROOT.TVector3()
            zetaU= muvec*(1/muvec.Mag()) +tauvec*(1/tauvec.Mag())
            pall=(muvec + tauvec + etvec).Dot(zetaU)
            pvis=(muvec+ tauvec).Dot(zetaU)
            pzeta=pall-0.85*(pvis)
            if pzeta > -25: return True
            else: return False
