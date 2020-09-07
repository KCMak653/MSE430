def recipLattice(recip_lat, a1):
    from pymatgen import Structure, Lattice
    import numpy as np
    import nglview as ngl
    import MSE430Funcs.CrysStrucFuncs as csfunc
    
    unit_cell = Structure(recip_lat, ['Cl'], [[0,0,0]], to_unit_cell=True, coords_are_cartesian=False)
    unit_cell_conv = csfunc.cubicCell(unit_cell, a1)
   
    selec=[]
    for ind, site in enumerate(unit_cell_conv.sites):
        if all(site.coords>=-a1*.15) & all(site.coords<=a1*1.15):
            selec=selec+[ind]
    view6 = ngl.show_pymatgen(unit_cell_conv)
    view6.clear_representations()
    view6.add_representation('point', aspectRatio=10, selection=selec)
    #[view6.add_representation('ball+stick', aspectRatio=5, selection=[i]) for i in selec]
 
    return(view6)

def brilZone(bz, recip_lat):
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    import matplotlib.pyplot as plt
    
    fig = plt.figure(figsize=(4,4))

    ax = fig.add_subplot(1, 1, 1, projection='3d')
    
    
    recip_vec=recip_lat.matrix

    for facet in bz:

        ab = [[tuple(it) for it in facet]]
        ax.add_collection3d(Poly3DCollection(ab, alpha=0.3, facecolor = 'red', edgecolor='k'))
    #ax.set_edgecolor('k')
    ax.set(zlim=[-4,4], xlim=[-4, 4], ylim=[-4,4], xlabel='x\'', ylabel='y\'', zlabel='z\'',
           xticklabels=[], yticklabels=[], zticklabels=[])

    ps = np.linspace(-1, 1, 3)
    qs = np.linspace(-1, 1, 3)
    rs = np.linspace(-1,1,3)
    pqr = np.array(np.meshgrid(ps, qs, rs)).reshape(3,-1)

    xyz2 = recip_vec.dot(pqr)

    ax.scatter(xyz2[0], xyz2[1], xyz2[2])
    
    plt.show()
    
