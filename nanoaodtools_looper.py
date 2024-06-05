#!/usr/bin/env python
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from importlib import import_module
import os
import sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import argparse
import glob
import subprocess
import plotting_util

def filename_num(s):
  s = s[:-5]
  i = s.rfind('_')
  num = s[i+1:]
  return int(num)

# command line options
parser = argparse.ArgumentParser(description="", usage="./%(prog)s INPUT")

# input/output
io_args = parser.add_argument_group('input/output options')
io_args.add_argument("input", metavar='INPUT',help="")
datamc_options = parser.add_mutually_exclusive_group()
datamc_options.add_argument("--data", action="store_true", default=False, help="running on data")
datamc_options.add_argument("--mc", action="store_true", default=False, help="running on bkg mc")
datamc_options.add_argument("--sigRes", action="store_true", default=False, help="running on resonant signal mc")
datamc_options.add_argument("--sigNonRes", action="store_true", default=False, help="running on nonresonant signal mc")
parser.add_argument("--out", default='analysis_out.root', help='name of root output file with histos (include .root)')
parser.add_argument("--loc", help='original location of atto for metadata')
parser.add_argument("--fast", action='store_true', default=False, help='test run')
parser.add_argument("--lumi", default=1.0, help='')
parser.add_argument("--cut", type=str, default='None', help='')
parser.add_argument("--photonchoice", default="HPID", help='')
parser.add_argument("--plotter", default="", help='')
parser.add_argument("--year", default="UL18", choices=["UL17", "UL18"], help='')
parser.add_argument("--phislice", type=bool, default=False, help='')
args = parser.parse_args()

# import modules
from PhysicsTools.NanoAODTools.fmk_plotting.sanityAnalysis import SanityAnalysis
from PhysicsTools.NanoAODTools.fmk_plotting.myAnalysis import MyAnalysis
from PhysicsTools.NanoAODTools.fmk_plotting.sigAnalysis import SigAnalysis
from PhysicsTools.NanoAODTools.fmk_plotting.trigAnalysis import TriggerAnalysis
#from PhysicsTools.NanoAODTools.fmk_plotting.HistogramProducerModVer import HistProd

if args.data: datamc = 'data'
elif args.mc: datamc = 'mc'
elif args.sigRes: datamc = 'sigRes'
elif args.sigNonRes: datamc = 'sigNonRes'
else: raise SystemExit('ERROR: Must specify one of --data / --mc / --sigRes / --sigNonRes !')

# append the year to the datamc string 
datamc += args.year[2:]

print("DATAMC: " + datamc)

files = []
if not args.loc and not args.plotter=="zttplot": metadata_chain = ROOT.TChain('Metadata')
if args.input == 'local':
    for fi in os.listdir("."):
        if fi == 'scalefactor.root': continue
        if fi.endswith(".root"):
            files.append(fi)
            if not args.loc and not args.plotter=="zttplot": metadata_chain.Add(fi)
else:
  if '.root' in args.input:
    files = [args.input]
    if not args.loc and not args.plotter=="zttplot": metadata_chain.Add(args.input)
  elif '/store/' in args.input:
    list_of_files = (subprocess.check_output("xrdfs root://cmseos.fnal.gov ls " + args.input, shell=True)).split('\n')
    list_of_files = [x for x in list_of_files if x]
    list_of_files.sort(key=filename_num)
    for line in list_of_files:
        files.append('root://cmseos.fnal.gov/'+line)
        if not args.loc and args.plotter=="zttplot": metadata_chain.Add('root://cmseos.fnal.gov/'+line)
  else:
    with open(args.input) as fi:
      for line in fi:
        # process line if .dat
        if args.input[-4:] == '.dat':
          newline = line.strip()
          i = newline.rfind('/')
          newline = newline[i+1:len(line)]
          files.append(newline)
        # dont process for .txt
        if args.input[-4:] == '.txt':
          files.append(line.strip())
        # set of eos locations for .eos
        if args.input[-4:] == '.eos':
          list_of_files = (subprocess.check_output("xrdfs root://cmseos.fnal.gov ls " + line.strip(), shell=True)).split('\n')
          list_of_files = [x for x in list_of_files if x]
          list_of_files.sort(key=filename_num)
          for fi in list_of_files:
            print fi
            files.append('root://cmseos.fnal.gov/'+fi)
            if not args.loc and not args.plotter=="zttplot": metadata_chain.Add('root://cmseos.fnal.gov/'+fi)

if args.loc and not args.plotter=="zttplot" and not args.data:
  d = plotting_util.get_meta(args.loc, jobdir=False)
  print(d)
  lookup_xs = {}; lookup_xs[d['dataset_id']] = d['xs']
  lookup_ngen = {}; lookup_ngen[d['dataset_id']] = d['evtProcessed']
  for key in lookup_xs:
    print(key, '->', lookup_xs[key])
  for key in lookup_ngen:
    print(key, '->', lookup_ngen[key])
  
if not args.loc and not args.plotter=="zttplot" and not args.data:
  lookup_xs = {}
  lookup_ngen = {}
  for event in metadata_chain:
    print(event.dataset)
    print(event.dataset_id)
    print(event.flag)
    print(event.xs)
    print(event.evtProcessed)
    print('')
    if int(event.dataset_id) not in lookup_xs:
      lookup_xs[int(event.dataset_id)] = float(event.xs)
    if int(event.dataset_id) not in lookup_ngen:
      lookup_ngen[int(event.dataset_id)] = int(event.evtProcessed)
    else:
      lookup_ngen[int(event.dataset_id)] += int(event.evtProcessed)
  for key in lookup_xs:
    print(key, '->', lookup_xs[key])
  for key in lookup_ngen:
    print(key, '->', lookup_ngen[key])

if args.data:
  lookup_xs = None
  lookup_ngen = None

modules = []
if args.plotter == 'sanity': modules += [SanityAnalysis(datamc, float(args.lumi), lookup_xs, lookup_ngen, args.cut, args.photonchoice)]
if args.plotter == 'bkg': modules += [MyAnalysis(datamc, float(args.lumi), lookup_xs, lookup_ngen, args.phislice)]
if args.plotter == 'sigeff': modules += [SigAnalysis(datamc, float(args.lumi), lookup_xs, lookup_ngen)]
if args.plotter == 'trig': modules += [TriggerAnalysis(float(args.lumi), lookup_xs, lookup_ngen)]
#if args.plotter == 'zttplot': modules += [HistProd(datamc)]

if args.fast: n = 1
else: n = None
  
p = PostProcessor('.',
                  files,
                  noOut=True,
                  histFileName=args.out,
                  histDirName='plots',
                  modules=modules,
                  maxEntries=n,
                  #totalEntries=args.numEvents,
                  )
p.run()
