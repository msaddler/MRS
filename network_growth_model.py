import numpy as np
import matplotlib.pyplot as mp

dT = 0.5
T = np.arange(0, 1, dT)

N = 10
inh_fract = 0.20
inh_flag = np.zeros(N)
rand_idx = np.arange(N)
np.random.shuffle(rand_idx)
inh_idx = rand_idx[0:round(inh_fract * N)]
exc_idx = [i for i in rand_idx if i not in inh_idx]
inh_flag[inh_idx] = 1

pos = np.random.uniform(size = (N, 2))
rad = np.zeros(shape = (N, len(T)))

print(rad)