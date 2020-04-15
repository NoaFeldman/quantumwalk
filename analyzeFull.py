from matplotlib import pyplot as plt
import pickle
import numpy as np
import math
import cmath
from matplotlib.animation import FuncAnimation

d = 9
dataX = [i for i in range(10 * 2*d )]
dataY = [i for i in range(int(-10 * (2**d / 2)), int(10 * (2**d / 2)))]

data = np.zeros((len(dataX), len(dataY)))

def getRowNum(index):
    return int(math.log(index + 1, 2))

def getSiteDataIndex(region: str, index: int):
    if region == 'tree1':
        rowNum = getRowNum(index)
        return [10 * rowNum, 10 * (index + 1 - 2**rowNum) + int(len(dataY) / 2) - int(10 * 2 ** (rowNum - 1))]
    if region == 'tree2':
        rowNum = getRowNum(index)
        return  [10 *  (2 * d - rowNum), 10 * (index + 1 - 2**rowNum) + int(len(dataY) / 2) - int(10 * 2 ** (rowNum - 1))]

from state import Graph, State
filename = 'l_0.1_21700.data'
with open(filename, 'rb') as filehandle:
    state = pickle.load(filehandle)

# import math
# from state import Graph, State
# for l in [0, 0.1, 0.2]:
#     amps = []
#     if l == 0:
#         iVals = [100 * i for i in range(int(1 * d * math.pi/1e-2 / 100))]
#     else:
#         iVals = [100 * i for i in range(int(1 * d * math.pi/(l * 1e-2) / 100))]
#     for i in iVals:
#         with open('vtestSameV/fullData/l_' + str(float(l)) + '_' + str(i) + '.data', 'rb') as filehandle:
#             state = pickle.load(filehandle)
#             amps.append(abs(state.tree2Coeffs[0])**2)
#     plt.plot([i / len(iVals) for i in iVals], amps)
# for l in [0, 0.1, 0.2]:
#     amps = []
#     if l == 0:
#         iVals = [100 * i for i in range(int(1 * d * math.pi/1e-2 / 100))]
#     else:
#         iVals = [100 * i for i in range(int(1 * d * math.pi/(l * 1e-2) / 100))]
#     for i in iVals:
#         with open('vtestSameV/fullData/l_' + str(float(l)) + '_' + str(i) + '.data', 'rb') as filehandle:
#             state = pickle.load(filehandle)
#             amps.append(abs(state.tree2Coeffs[2**d - 6])**2 * 2**d)
#     plt.plot([i / len(iVals) for i in iVals], amps)
# plt.show()

for i in range(len(state.tree1Coeffs)):
    x, y = getSiteDataIndex('tree1', i)
    for j in range(x, x + 10):
        for k in range(y, y + 10):
            data[j, k] = abs(state.tree1Coeffs[i])**2 * 2**getRowNum(i)
for i in range(len(state.tree2Coeffs)):
    getSiteDataIndex('tree2', i)
    for j in range(x, x + 10):
        for k in range(y, y + 10):
            data[j, k] = abs(state.tree2Coeffs[i])**2 * 2**getRowNum(i)

plt.pcolormesh([i for i in range(len(dataY))], [i for i in range(len(dataX))], data)
plt.colorbar()
plt.show()
