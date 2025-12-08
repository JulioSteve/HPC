import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.animation as anim

# Load data
Tsteps, N, mod = np.loadtxt("settings.dat", unpack=True, skiprows=1, usecols=(0,1,3), dtype=int)
dt = np.loadtxt("settings.dat", unpack=True, skiprows=1, usecols=(2), dtype=np.float32)

t = np.array([x*dt for x in range(0,Tsteps+1, mod)]) #We tested the reconstruction to ensure it is the right formula

POS = np.fromfile("pos.bin", dtype=np.float32).reshape((len(t),N,3), order="F")
VEL = np.fromfile("vel.bin", dtype=np.float32).reshape((len(t),N,3), order="F")


# # Setup figure
# plt.rc('figure', facecolor="black")
# plt.rc('axes', facecolor='black', edgecolor='red', labelcolor='red')
# plt.rc('grid', color='grey')
# plt.rc('text', color='green')

# f = plt.figure()
# ax = f.add_subplot(projection='3d')

# # Initial scatter
# scatp = dict(s=100, c='gold', marker="*")
# scat = ax.scatter(x[:N], y[:N], z[:N], **scatp)

# # Animation update function
# def update(frame):
#     start = frame * N
#     end = start + N
#     X = x[start:end]
#     Y = y[start:end]
#     Z = z[start:end]
#     # print(f"Frame {frame}: X = {X}, Y = {Y}, Z = {Z}")  # Debugging
#     progress = (frame+1)/timesteps*100
#     print(f"\rGenerating animation... {progress:5.1f}%", end='')
#     scat._offsets3d = (X, Y, Z)
#     return scat,

# # Aesthetics
# ax.xaxis.set_pane_color((0, 0, 0, 1))
# ax.yaxis.set_pane_color((0, 0, 0, 1))
# ax.zaxis.set_pane_color((0, 0, 0, 1))

# # Working frame
# boxlim = (-1,1)
# ax.set(xlim=boxlim,ylim=boxlim,zlim=boxlim)

# ax.tick_params(colors='white')

# # Create animation (store reference in a variable to avoid garbage collection)
# freq = 60 #FPS (Hz) 
# tau = 1000/freq#interval in ms between images of animation.
# # ani = anim.FuncAnimation(fig=f, func=update, frames=timesteps, interval=tau, blit=False)

# # Generating a mp4 video of the simulation | or | showing the animation in live to user.
# animtitle = f"N{N}_T5e3_dt5e-4_NOeps_60fps"
# # ani.save(animtitle+".gif", writer="pillow",fps=freq)
# # ani.save(animtitle+".avi", writer="ffmpeg", bitrate=5000,fps=freq)
# # plt.show()

# # Generating the plot of the trajectories, fixed on the center of mass:
# # t2, x2, y2, z2 = np.loadtxt("2body.dat", skiprows=1,unpack=True, usecols=[0,1,2,3])
# # def traj2(T,X,Y,Z):
# #     pos1 = []
# #     pos2 = []
# #     for i in range(len(T)):
# #         if i%2==0:
# #             pos1.append([X[i],Y[i],Z[i]])
# #             print(pos1)
# #         else:
# #             pos2.append([X[i],Y[i],Z[i]])
    
# #     TRAJ = (pos1[0],pos1[1],pos2[0],pos2[1])
# #     CM = np.mean(TRAJ)
# #     return CM

# # CM = traj2(t2, x2, y2, z2)

# # print(CM)





# print("") #Just so the percentage in animation doesn't ruin the terminal layout.