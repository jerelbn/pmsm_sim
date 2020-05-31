import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(linewidth=150)
plt.rc("text", usetex=True)
plt.rc("font", family="serif")

x = np.fromfile("/tmp/pmsm.log", dtype=np.float64).reshape(-1, 5) # [t,id,iq,theta,omega]

fig, axes = plt.subplots(2, 2)
plt.suptitle("PMSM States")

axes[0,0].grid(True)
axes[0,0].set_ylabel("Direct Current (A)")
axes[0,0].set_xlabel("Time (s)")
axes[0,0].plot(x[:,0], x[:,1], color="b", linewidth=2.0, label="id")
axes[0,0].legend()

axes[0,1].grid(True)
axes[0,1].set_ylabel("Quadrature Current (A)")
axes[0,1].set_xlabel("Time (s)")
axes[0,1].plot(x[:,0], x[:,2], color="b", linewidth=2.0, label="iq")
axes[0,1].legend()

axes[1,0].grid(True)
axes[1,0].set_ylabel("Angle (rad)")
axes[1,0].set_xlabel("Time (s)")
axes[1,0].plot(x[:,0], x[:,3], color="b", linewidth=2.0, label="angle")
axes[1,0].legend()

axes[1,1].grid(True)
axes[1,1].set_ylabel("Rate (rad/s)")
axes[1,1].set_xlabel("Time (s)")
axes[1,1].plot(x[:,0], x[:,4], color="b", linewidth=2.0, label="rate")
axes[1,1].legend()

plt.show()
