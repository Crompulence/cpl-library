import numpy as np 
import matplotlib.pyplot as plt

dragModels = ['DiFelice', 'Ergun']
Uf_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
fig = plt.figure(figsize=(12.8, 9.6))
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)
plt.ion()
plt.show()
for Uf in Uf_values:
    data = np.load('./results/fig_Uf_{}_DiFelice_displacement.npz'.format(Uf))
    t = data['t']
    xy = data['xy']
    xySol = data['xySol']

    ax1.plot(t, xy, linewidth=3.0, label='Uf = {}'.format(Uf))
    ax1.plot(t, np.ones_like(t)*xySol, 'k--')
    plt.pause(0.1)

    data = np.load('./results/fig_Uf_{}_Ergun_displacement.npz'.format(Uf))
    t = data['t']
    xy = data['xy']
    xySol = data['xySol']

    ax2.plot(t, xy, linewidth=3.0)
    ax2.plot(t, np.ones_like(t)*xySol, 'k--')
    plt.pause(0.1)

ax1.set_title('Di Felice')
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Position (cm)')
ax1.set_ylim((9.905,9.915))
ax1.legend()
ax2.set_title('Ergun')
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Position (cm)')
ax2.set_ylim((9.905,9.915))
plt.tight_layout()
plt.savefig('fig_displacement.png')
plt.close()



dragModels = ['DiFelice', 'Ergun']
Uf_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
fig = plt.figure(figsize=(12.8, 9.6))
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)
plt.ion()
plt.show()
for Uf in Uf_values:
    data = np.load('./results/fig_Uf_{}_DiFelice_pressure.npz'.format(Uf))
    h = data['h']
    p = data['p']
    pSol = data['pSol']
    ylo = data['ylo']
    yhi = data['yhi']

    ax1.plot(h, p, linewidth=3.0, label='Uf = {}'.format(Uf))
    ax1.plot([ylo,yhi], [pSol,0.], 'k--')
    plt.pause(0.1)

    data = np.load('./results/fig_Uf_{}_Ergun_pressure.npz'.format(Uf))
    h = data['h']
    p = data['p']
    pSol = data['pSol']
    ylo = data['ylo']
    yhi = data['yhi']

    ax2.plot(h, p, linewidth=3.0, label='Uf = {}'.format(Uf))
    ax2.plot([ylo,yhi], [pSol,0.], 'k--')
    plt.pause(0.1)

ax1.set_title('Di Felice')
ax1.set_xlabel('Height (cm)')
ax1.set_ylabel('Pressure (0.1Pa)')
ax1.legend()
ax2.set_title('Ergun')
ax2.set_xlabel('Height (cm)')
ax2.set_ylabel('Pressure (0.1Pa)')
plt.tight_layout()
plt.savefig('fig_pressure.png')
plt.close()