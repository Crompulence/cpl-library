import numpy as np 
import matplotlib.pyplot as plt

dragModels = ['Stokes']
dp_values = [0.01, 0.02, 0.03, 0.04]
fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)
plt.ion()
plt.show()
for dp in dp_values:
    for dragModel in dragModels:
        data = np.load('./results/fig_dp_{}_{}_displacement.npz'.format(dp, dragModel))
        t = data['t']
        xy = data['xy']
        xySol = data['xySol']

        data = np.load('./results/fig_dp_{}_{}_velocity.npz'.format(dp, dragModel))
        vy = data['vy']
        vySol = data['vySol']
        vyTer = data['vyTer']

        ax1.plot(t, xy, linewidth=3.0, label='dp = {}'.format(dp))
        ax1.plot(t, xySol, 'k--')
        plt.pause(0.1)

        ax2.plot(t, vy, linewidth=3.0, label='dp = {}'.format(dp))
        ax2.plot(t, vySol, 'k--')
        ax2.plot(t, vyTer, 'k:')
        plt.pause(0.1)

ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Position (cm)')
ax1.legend()
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Velocity (cm/s)')
plt.tight_layout()
plt.savefig('fig_displacement_velocity.png')
plt.close()

