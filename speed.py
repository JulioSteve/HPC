import matplotlib.pyplot as plt

T = [0.675, 0.353, 0.246, 0.187, 0.159, 0.135, 0.122, 0.114, 0.131, 0.127, 0.120, 0.122, 0.127, .109, 0.122, 0.114]
C = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]

plt.plot(C,T, color="black", lw=3, marker=".", ls="solid", ms=10)
plt.xlabel('Threads')
plt.ylabel('Time (s)')
plt.show()