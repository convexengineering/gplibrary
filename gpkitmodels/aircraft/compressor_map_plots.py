from engine_components import FanMap, LPCMap
from gpkit import Model, units
import numpy as np
import matplotlib.pyplot as plt
from engine import EngineOnDesign

pivec = []
nvec = []

engineOnD = EngineOnDesign()
    
sol = engineOnD.localsolve(verbosity = 1, kktsolver="ldl")

mhtD, mltD, NlpcD, NhpcD, A5, A7 = engineOnD.sizing(sol)

mvec = np.linspace(20,50,2)
#mvec = 20

nvec = np.linspace(1,20,1)
##nvec = 1

lpcmap = LPCMap()
lpcmap.substitutions.update({
    'T_{ref}': 1000,
    'T_{t_2}': sol('T_{t_2}'),
    'm_{core}': ('sweep', mvec),
    'm_{core_D}': sol('m_{core}'),
    'P_{ref}': 22,
    'P_{t_2}': sol('P_{t_2}'),
    '\pi_{lc_D}': sol('\pi_{lc}'),
    'N_{{bar}_D_lc}': 1,
    'N_1': ('sweep',nvec)
    })
sol = lpcmap.localsolve(kktsolver="ldl", verbosity = 4, skipsweepfailures=True)


