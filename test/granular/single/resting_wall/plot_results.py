import numpy as np 
import matplotlib.pyplot as plt

dragModel = 'DiFelice'
y0 = 1.0

dp_values = [0.05, 0.06, 0.07, 0.08, 0.09, 0.10]
Uf_values = [0., 40.]
fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)
plt.ion()
plt.show()
for dp in dp_values:
    
    data = np.load('./results/fig_{}_dp_{}_Uf_{}_y0_{}.npz'.format(dragModel, dp, Uf_values[0], y0))
    t = data['t']
    xy = data['xy']
    xySol = data['xySol']

    ax1.plot(t, xy-xy[0], linewidth=3.0, label='Uf = {:.0f}, dp = {:.2f}'.format(Uf_values[0], dp))
    ax1.plot(t, (xySol-xy[0])*np.ones_like(t), 'k:')
    plt.pause(0.1)

    data = np.load('./results/fig_{}_dp_{}_Uf_{}_y0_{}.npz'.format(dragModel, dp, Uf_values[1], y0))
    t = data['t']
    xy = data['xy']
    xySol = data['xySol']

    ax2.plot(t, xy-xy[0], linewidth=3.0, label='Uf = {:.0f}, dp = {:.2f}'.format(Uf_values[1], dp))
    ax2.plot(t, (xySol-xy[0])*np.ones_like(t), 'k--')
    plt.pause(0.1)

ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Position (cm)')
ax1.set_ylim((-0.0009, 0.0001))
ax1.legend(fontsize='xx-small', loc='upper right', ncol=2)
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Position (cm)')
ax2.set_ylim((-0.0009, 0.0001))
# ax2.legend(fontsize='xx-small')
plt.tight_layout()
plt.savefig('fig_displacement.png')
plt.close()

