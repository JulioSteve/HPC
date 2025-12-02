import numpy as np
import matplotlib.pyplot as plt

Times = np.loadtxt("pi.dat", usecols=(1), delimiter=",")*1e-3
Times_2 = np.loadtxt("pi2.dat", usecols=(1), delimiter=",")*1e-3
N = np.arange(1,17)

plt.xlabel('Number of threads')
plt.ylabel('Time (s)')
plt.xlim(0,17)
plt.hlines(0.106, -1, 18, alpha=0.6, color="grey", ls="--", lw=4, label="No //")
plt.plot(N,Times, marker=".", ms=15, lw=4, color="black", label="Raw //")
plt.plot(N,Times, marker="x", ms=15, lw=4, ls=":",color="red", label="Reduction //")

plt.legend(loc="center right")
plt.show()