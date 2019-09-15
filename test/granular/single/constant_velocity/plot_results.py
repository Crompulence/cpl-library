import numpy as np 
import matplotlib.pyplot as plt

dragModel = 'DiFelice'
dp_values = [0.01, 0.02, 0.03]
vy0_values = [0.5, 1.0, 1.5, 2.0, 2.5]

plt.figure()
plt.ion()
plt.show()
for vy0 in vy0_values:
    plotData = np.zeros((len(dp_values),3))
    counter = 0
    for dp in dp_values:
        file_name = './results/data_{}_dp_{}_vy0_{}.npz'.format(dragModel, dp, vy0)
        data = np.load(file_name)
        plotData[counter,:] = [dp, data['fy'], data['fySol']]
        counter += 1

    plt.plot(plotData[:,0], plotData[:,1], linewidth=3.0, label='vy0 = {}'.format(vy0))
    plt.plot(plotData[:,0], plotData[:,2], 'k--')
    plt.pause(0.1)

plt.xlabel('Particle diameter (cm)')
plt.ylabel('Force')
plt.legend()
plt.tight_layout()
plt.savefig('fig_force.png')
plt.close()

