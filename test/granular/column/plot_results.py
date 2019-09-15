import numpy as np 
import matplotlib.pyplot as plt

dragModels = ['DiFelice', 'Ergun']
Uf_values = [0.1, 0.2, 0.3, 0.4, 0.5]
fig = plt.figure()
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)
plt.ion()
plt.show()
for Uf in Uf_values:
    data = np.load('./results/fig_pressure_Uf_{}_DiFelice_moving.npz'.format(Uf))
    h = data['h']
    p = data['p']
    pSol = data['pSol']
    xyz_orig = data['xyz_orig']
    xyzL = data['xyzL']

    ax1.plot(h, p, linewidth=3.0, label='Uf = {}'.format(Uf))
    ax1.plot([xyz_orig[1],xyzL[1]], [pSol,0.], 'k--')
    plt.pause(0.1)

    data = np.load('./results/fig_pressure_Uf_{}_Ergun_moving.npz'.format(Uf))
    h = data['h']
    p = data['p']
    pSol = data['pSol']
    xyz_orig = data['xyz_orig']
    xyzL = data['xyzL']

    ax2.plot(h, p, linewidth=3.0, label='Uf = {}'.format(Uf))
    ax2.plot([xyz_orig[1],xyzL[1]], [pSol,0.], 'k--')
    plt.pause(0.1)

ax1.set_title('Di Felice')
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Position (cm)')
ax1.legend()
ax2.set_title('Ergun')
ax2.set_xlabel('Height (cm)')
ax2.set_ylabel('Pressure (0.1Pa)')
plt.tight_layout()
plt.savefig('fig_pressure.png')
plt.close()

