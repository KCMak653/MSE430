"""Contains functions used for DOS module"""

def DOS(dim):
    import numpy as np
    Es = np.linspace(-2, 10, 600)
    m_e = 0.8 * 9.11e-31 #kg
    hbar = 1.054e-34 # m^2kg/s
    if dim == 3:
        gE = 1/(2*np.pi**2)*(2*m_e/hbar**2)**(3/2)*np.sqrt(abs(Es))       
    if dim == 2:
        gE= m_e/(np.pi*hbar**2)*np.ones(len(abs(Es)))
    if dim == 1:
        gE = m_e/(np.pi*hbar)*np.sqrt(m_e/(2*abs(Es)))
    gE[Es<0]=np.nan   
    return(gE, Es)    

def center_axis(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])

    #Redraw axis
    labls = [r'$k_x$', r'$k_y$', r'$k_z$']
    val = [1,0, 0]
    a1 =[-12, 12]
    a2 = [0, 0]


    ax.plot(a1, a2, a2, color = 'black', lw = 3, alpha=0.8)
    ax.plot(a2, a1, a2, color = 'black', lw = 3, alpha =0.8)
    ax.plot(a2, a2, a1, color = 'black', lw = 3, alpha = 0.8)
    ax.text(-14, 0, 0, labls[0])
    ax.text(0, -14, 0, labls[1])
    ax.text(0, 0, -14, labls[2])
    
def k_diag(dim):
    import numpy as np
    pts = np.linspace(-10, 10, 11)
    if dim == 3:
        #plot grid points
        (kx, ky, kz) = np.meshgrid(pts, pts, pts)
        #ax1.scatter(kx, ky, kz, s=0.7)

        #make surface plot
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)

        x = 10 * np.outer(np.cos(u), np.sin(v))
        y = 10 * np.outer(np.sin(u), np.sin(v))
        z = 10 * np.outer(np.ones(np.size(u)), np.cos(v))

        #vector
        k_vec=np.sqrt(100/dim)
        (k_vec_x, k_vec_y, k_vec_z)=[[0, k_vec], [0,k_vec], [0,k_vec]]
        (kx_t, ky_t, kz_t) = [k_vec, k_vec/2 +2, k_vec/2]

    if dim ==2:
        #plot grid points
        (kx, ky, kz) = np.meshgrid(pts, pts, np.zeros(len(pts)))
        #ax1.scatter(kx,ky, kz, s=0.7)

        #make surface plot
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)

        x = 10 * np.outer(np.cos(u), np.sin(v))
        y = 10 * np.outer(np.sin(u), np.sin(v))
        z = 0 * np.outer(np.ones(np.size(u)), np.cos(v))

        #ax1.plot_surface(x, y, z,  rstride=4, cstride=4, color='b', linewidth=0, alpha=0.3)

        #add vector
        k_vec=np.sqrt(100/dim)
        print(k_vec)
        (k_vec_x, k_vec_y, k_vec_z)=[[0, k_vec], [0,k_vec], [0,0]]
        (kx_t, ky_t, kz_t) = [k_vec/2, k_vec/2+2, 0]  


    if dim==1:
        #plot grid points
        kx = pts
        ky = np.ones(len(pts))
        kz = np.zeros(len(pts))
        #ax1.scatter(kx,ky, s=0.7)
        x=np.zeros([2,2])
        y=x
        z=x
        #add vector
        (k_vec_x, k_vec_y, k_vec_z)=[[0, 10], [1,1], [0,0]]
        (kx_t, ky_t, kz_t) = [5, 3, 0]   
    return(kx, ky, kz, k_vec_x, k_vec_y, k_vec_z, kx_t, ky_t, kz_t, x, y, z)
def makeVector(k_vec_x, k_vec_y, k_vec_z,ax):
    from matplotlib.patches import FancyArrowPatch
    class Arrow3D(FancyArrowPatch):
        from matplotlib.patches import FancyArrowPatch
        def __init__(self, xs, ys, zs, *args, **kwargs):
            FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
            self._verts3d = xs, ys, zs

        def draw(self, renderer):
            xs3d, ys3d, zs3d = self._verts3d
            xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
            self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
            FancyArrowPatch.draw(self, renderer)
    a = Arrow3D(k_vec_x, k_vec_y, k_vec_z, mutation_scale=20,
             lw=1, arrowstyle="-|>", color="r")
    ax.add_artist(a)

def graphProp(ax1, ax2, dim):
    ymx = ax2.get_ylim()
    multy =[1,2,1]
    ax2.set(title = 'Density of States in {}D'.format(dim), xlabel = 'E', ylabel = r'$\rho (E)$', 
        yticks =[], xticks=[0], xticklabels=[r'$E_C$'], xlim=[-2,10], ylim=[0, ymx[1]*multy[dim-1]])

    dep = ['E^{-1/2}', 'E^0', 'E^{1/2}']
    ax2.text(6, ymx[1]/2, r'$\propto~~{}$'.format(dep[dim-1]))
