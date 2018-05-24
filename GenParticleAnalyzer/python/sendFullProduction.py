#! /usr/bin/env python
import os
import sys


if (len(sys.argv) != 2):
    print "usage sendFullProduction.py prodName"
    sys.exit(1)

prodName = sys.argv[1]


os.system('source /afs/cern.ch/work/p/pandolf/CMSSW_9_4_1_CMG/src/Daphne/GenParticleAnalyzer/python/setAAA.sh')
os.system('eos mkdir /eos/cms/store/user/pandolf/Dafne/GenAnalysis/'+prodName)
os.system('python sendOnBatch.py '+prodName+' QCD_Pt_15to30 2')
