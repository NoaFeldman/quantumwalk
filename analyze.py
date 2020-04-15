from matplotlib import pyplot as plt
import pickle
import numpy as np

d = 9
res = []
lambdas = [0.01 * i for i in range(21)]
Es = [0.1 * e for e in range(-20, 21)]
for disorderParam in lambdas:
    dRes = []
    for E in Es:
        strE = str(E)
        if strE[0] == '-':
            strE = strE[:4]
        else:
            strE = strE[:3]
        with open('vtestSameV/randomTrees_' + str(d) + '_' + str(float(0)) +
          '_' + str(round(disorderParam, 2)) + '_' + strE + '.data', 'rb') as filehandle:
            amp = pickle.load(filehandle)
            dRes.append(amp)
    res.append(dRes)
plt.pcolormesh(Es, lambdas, res)
plt.colorbar()
plt.show()
