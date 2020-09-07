def setAxProp(recip_vec, real_vec):
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    
    #Set up plots
    fig1 = plt.figure(figsize=(8,4))
    ax1 = fig1.add_subplot(1, 2, 1, projection='3d')
    ax2 = fig1.add_subplot(1,2,2, projection='3d')
    
    ps = np.linspace(-5, 5, 6)
    qs = np.linspace(-5, 5, 6)
    rs = np.linspace(-5,5,6)

    pqr = np.array(np.meshgrid(ps, qs, rs)).reshape(3,-1)
    xyz = real_vec.dot(pqr)
    xyz2 = recip_vec.dot(pqr)

    #zs = xy[0,:]
    ax1.scatter(xyz[0], xyz[1], xyz[2])
    ax2.scatter(xyz2[0], xyz2[1], xyz2[2])
    
    ax1.set(xlim=(-15, 15), ylim=(-15, 15), zlim=(-15, 15),
            xlabel='x', ylabel='y', zlabel='z', 
            xticklabels=[], yticklabels=[], zticklabels=[],
           title = 'Real Lattice')
    ax2.set(xlim=(-15, 15), ylim=(-15, 15), zlim=(-15, 15),
            xlabel='x\'', ylabel='y\'', zlabel='z\'', 
            xticklabels=[], yticklabels=[], zticklabels=[],
           title = 'Reciprocal Lattice')
    plt.show()