#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

nis_las = np.loadtxt('build/nis_las.txt', delimiter='\n')
nis_rad = np.loadtxt('build/nis_rad.txt', delimiter='\n')
f, (ax1, ax2) = plt.subplots(2, 1)
n = len(nis_las)
ax1.plot(range(n), nis_las)
ax1.plot(range(n), [7.815] * n)
n = len(nis_rad)
ax2.plot(range(n), nis_rad)
ax2.plot(range(n), [7.815] * n)
plt.show()

# df = pd.read_csv('build/result.txt')
# plt.plot(df['gt_px'], df['gt_py'], 'b-')
# plt.plot(df['est_px'], df['est_py'], 'r-')
# plt.show()
