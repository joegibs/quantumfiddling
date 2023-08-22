# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 16:58:42 2023

@author: jogib
"""

import quimb as qu
import numpy as np
#%%
state=qu.ghz_state(2,qtype='dop')
qu.partial_trace(state,[2]*2,1)
evv,evec=np.linalg.eigh(state)

print(evv)
print(sum(evv))

#%%
state_p=qu.computational_state('11',qtype='dop')+qu.computational_state('00',qtype='dop')
state_p=state_p/norm(state_p)
evv,evec=np.linalg.eigh(state_p)
print(sum(evv))