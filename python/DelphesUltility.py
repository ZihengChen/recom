import uproot
import numpy as np
import pandas as pd
from collections import defaultdict
from skhep.math.vectors import LorentzVector



BASEDIR = "/Users/zihengchen/Documents/Graviton"


def fill_kinematics(out, lep1, lep2, gamma, jet1, jet2):
    out['leptonOne_pt' ].append(lep1.pt)
    out['leptonOne_eta'].append(lep1.eta)
    out['leptonTwo_pt' ].append(lep2.pt)
    out['leptonTwo_eta'].append(lep2.eta)
    out['gamma_pt'     ].append(gamma.pt)
    out['gamma_eta'    ].append(gamma.eta)
    out['jetOne_pt'    ].append(jet1.pt)
    out['jetOne_eta'   ].append(jet1.eta)
    out['jetTwo_pt'    ].append(jet2.pt)
    out['jetTwo_eta'   ].append(jet2.eta)


    dilep = lep1 + lep2
    out['dilepton_pt' ].append(dilep.pt)
    out['dilepton_eta'].append(dilep.eta)
    out['dilepton_m'].append(dilep.m)
    out['dilepton_deltar'].append(lep1.deltar(lep2))
    out['dilepton_deltaphi'].append(lep1.deltaphi(lep2))


    dilepg = dilep + gamma
    out['dileptonG_pt' ].append(dilepg.pt)
    out['dileptonG_eta'].append(dilepg.eta)
    out['dileptonG_m'].append(dilepg.m)
    out['dileptonG_deltar'].append(dilep.deltar(gamma))
    out['dileptonG_deltaphi'].append(dilep.deltaphi(gamma))   


    dijet = jet1 + jet2
    out['dijet_pt' ].append(dijet.pt)
    out['dijet_eta'].append(dijet.eta)
    out['dijet_m'].append(dijet.m)
    out['dijet_deltar'].append(jet1.deltar(jet2))
    out['dijet_deltaphi'].append(jet1.deltaphi(jet2))


    dilepgDijet = dijet + dilepg
    out['dileptonGDijet_pt' ].append(dilepgDijet.pt)
    out['dileptonGDijet_eta'].append(dilepgDijet.eta)
    out['dileptonGDijet_m'].append(dilepgDijet.m)
    out['dileptonGDijet_deltar'].append(dilepg.deltar(dijet))
    out['dileptonGDijet_deltaphi'].append(dilepg.deltaphi(dijet))