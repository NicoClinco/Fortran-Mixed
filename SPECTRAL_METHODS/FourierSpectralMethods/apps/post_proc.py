
"""
 File used for doing post-processing of the field
"""

if __name__=='__main__':
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    #sol = pd.read_csv('Solution.csv',header=None)
    sol = pd.read_csv('sol_euler.dat',header=None,delim_whitespace=True,dtype=float)
    
    sol_numpy = sol.to_numpy()
    
    fig,axs = plt.subplots(2,1,figsize=(7,7))
    num_points = sol_numpy.shape[1]
    x = np.linspace(0,2*np.pi,num_points)
    axs[0].plot(x,sol_numpy[0,:],label='Solution t=0s')
    print(sol_numpy[0,:])
    #axs[0].plot(x,np.exp(-100.0*(x-1.0)**2),label='Initial exact solution')
    axs[1].plot(x,sol_numpy[-1,:],label='Solution t=2s')
    [axs[i].set_xlabel('x') for i in range(2)]
    [axs[i].set_ylabel('u(t,x)') for i in range(2)]
    [axs[i].legend(loc='best') for i in range(2)]
    [axs[i].grid(True) for i in range(2)]
    plt.show()
    
