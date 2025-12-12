import numpy as np 
import matplotlib

matplotlib.use("Qt5Agg")

import matplotlib.pyplot as plt
import matplotlib.animation as anim

# Load data
Tsteps, N, mod = np.loadtxt("settings.dat", unpack=True, skiprows=1, usecols=(0,1,3), dtype=int)
dt = np.loadtxt("settings.dat", unpack=True, skiprows=1, usecols=(2), dtype=np.float32)

t = np.array([x*dt for x in range(0,Tsteps+1, mod)]) #We tested the reconstruction to ensure it is the right formula

POS = np.fromfile("pos.bin", dtype=np.float32).reshape((len(t),N,3))
VEL = np.fromfile("vel.bin", dtype=np.float32).reshape((len(t),N,3))
ENERG = np.fromfile("energies.bin", dtype=np.float32).reshape((len(t),2))
LTOT = np.fromfile("angularmomentum.bin", dtype=np.float32)

def animation(flag):
    if flag:
        # Setup figure
        plt.rc('figure', facecolor="black")
        plt.rc('axes', facecolor='black', edgecolor='red', labelcolor='red')
        plt.rc('grid', color='grey')
        plt.rc('text', color='green')

        f = plt.figure(figsize=(2560/300,1440/300), dpi=300)
        ax = f.add_subplot(projection='3d')
        ax.set_xlabel(r"$x$ (reduced units)")
        ax.set_ylabel(r"$y$ (reduced units)")
        ax.set_zlabel(r"$z$ (reduced units)")

        # Initial scatter
        scatp = dict(s=20, c='gold', marker="*")
        x = POS[0,:,0]
        y = POS[0,:,1]
        z = POS[0,:,2]
        scat = ax.scatter(x, y, z, **scatp)

        # Animation update function
        def update(frame):
            start = frame * N
            end = start + N
            X = POS[frame,:,0]
            Y = POS[frame,:,1]
            Z = POS[frame,:,2]
            # print(f"Frame {frame}: X = {X}, Y = {Y}, Z = {Z}")  # Debugging
            progress = (frame+1)/len(t)*100
            print(f"\rGenerating animation... {progress:5.1f}%", end='')
            scat._offsets3d = (X, Y, Z)
            f.suptitle(f"T = {t[frame]:.2f} s\n N={N:,} bodies")
            return scat,

        # Aesthetics
        ax.xaxis.set_pane_color((0, 0, 0, 1))
        ax.yaxis.set_pane_color((0, 0, 0, 1))
        ax.zaxis.set_pane_color((0, 0, 0, 1))

        # Working frame
        lim = 2
        boxlim = (-lim,lim)
        ax.set(xlim=boxlim,ylim=boxlim,zlim=boxlim)
        ticks_uniform = np.linspace(-lim,lim,5)
        tu_labels = [str(t) for t in ticks_uniform]
        ax.set_xticks(ticks=ticks_uniform, labels=tu_labels)
        ax.set_yticks(ticks=ticks_uniform, labels=tu_labels)
        ax.set_zticks(ticks=ticks_uniform, labels=tu_labels)

        ax.tick_params(colors='white')

        # Create animation (store reference in a variable to avoid garbage collection)
        
        ani = anim.FuncAnimation(fig=f, func=update, frames=len(t))

        # Generating a mp4 video of the simulation or showing the animation
        animtitle = f"N{N}"
        
        ani.save(animtitle+".mp4", writer="ffmpeg",fps=60, dpi=300)
        # plt.show()

        # Generating the plot of the trajectories, fixed on the center of mass:
        # t2, x2, y2, z2 = np.loadtxt("2body.dat", skiprows=1,unpack=True, usecols=[0,1,2,3])
        # def traj2(T,X,Y,Z):
        #     pos1 = []
        #     pos2 = []
        #     for i in range(len(T)):
        #         if i%2==0:
        #             pos1.append([X[i],Y[i],Z[i]])
        #             print(pos1)
        #         else:
        #             pos2.append([X[i],Y[i],Z[i]])
            
        #     TRAJ = (pos1[0],pos1[1],pos2[0],pos2[1])
        #     CM = np.mean(TRAJ)
        #     return CM

        # CM = traj2(t2, x2, y2, z2)

        # print(CM)

        print("") #Just so the percentage in animation doesn't ruin the terminal layout.


def plots(plot, save):
    if (plot or save):
        fig, ([energ,angmom],[Zenerg,Zangmom]) = plt.subplots(
            2,2, layout="constrained", sharex=True, figsize=(2560/300, 1440/300), dpi=300
            )
        fig.supxlabel("Time (s)", color="white")
        fig.patch.set_facecolor("black")
        for ax in [energ,angmom,Zenerg,Zangmom]:
            ax.set_facecolor('dimgray')
            ax.tick_params(axis="both", colors="white")

        energ.plot(t,ENERG[:,0], label=r"$E_k$", lw=4, color="indianred")
        energ.plot(t,ENERG[:,1], label=r"$E_p$", lw=4, color="goldenrod")
        energ.plot(t,ENERG[:,0]+ENERG[:,1], label=r"$E_{tot}$", lw=4, color="mediumpurple")
        energ.set_title("Conservation of total energy:", color="white")
        energ.legend(loc="upper right")
        energ.set_ylabel("Energy (reduced unit)", color="white")

        ymin = 0.99*LTOT.min()
        ymax = 1.01*LTOT.max()
        angmom.set_ylim(ymin,ymax)
        angmom.plot(t,LTOT, label=r"$L_0$", lw=4, color="paleturquoise")
        angmom.set_title(r"Conservation of total angular momentum $L_0$", color="white")
        angmom.legend(loc="upper right")
        angmom.set_ylabel(r"$||\vec{L_0}||$ (reduced unit)", color="white")

        Zenerg.plot(t,ENERG[:,0]+ENERG[:,1], label=r"$E_{tot}$", lw=1, color="mediumpurple")
        Zenerg.grid(color="darkgrey")
        Zenerg.set_ylabel("Energy (reduced unit)", color="white")

        Zangmom.plot(t,LTOT, label=r"$L_0$", lw=1, color="paleturquoise")
        Zangmom.grid(color="darkgray")
        Zangmom.set_ylabel(r"$||\vec{L_0}||$ (reduced unit)", color="white")

        if plot:
            plt.show()
        if save:
            # plt.savefig("Conservation.png", format='png', dpi=300)
            plt.savefig("Conservation.pdf", format="pdf", bbox_inches="tight")

animation(True)
plots(plot=False, save=False)
