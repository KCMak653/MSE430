def potProf(b, V):
    import numpy as np
    
    d=1e-9 #m
    a=d-b

    #Draw the potential barriers

    x = np.linspace(0, a+b, 200)
    Vs = np.zeros(len(x))
    Vs[x>a]=V
    x = np.linspace(0, 5*(a+b), 1000)
    Vs = np.tile(Vs, 5)
    return(x, Vs, a, d)

def graphSetUp(ax1, a, b, d, V0, Etot):
    import matplotlib.pyplot as plt
    ax1.set_yticks([0, V0, 30*Etot])
    ax1.set_yticklabels(['V=0', r'$V_0$', r'$30*E_{tot}$'])

    ax1.set_xticks([0, a, a+b, a+d])
    ax1.set_xticklabels(['0', 'a', 'a+b', 'a+d'])
    ax1.set_title('Figure A: Kronig-Penney Model - Potential Profile')
    ax1.set(xlabel='x', ylabel='V(x)')

def plotRHS(b, V0):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle
    
    hbar=1.05e-34
    d=1e-9
    m_e=9.11e-31
    a = d-b

    fig2 = plt.figure()
    ax2 = plt.axes()


    rep = 200
    npts = rep * 24+25
    #Independent variable: alpha * a 
    alpha = np.linspace(-6*np.pi/a, 6*np.pi/a, npts)
    alpha_a = alpha * a

    E = hbar**2 * alpha**2/(2*m_e)

    B = np.sqrt(2*m_e/hbar**2 * (V0 - E))

    RHS = np.cos(alpha_a)*np.cosh(b*B)-(alpha**2 - B**2)/(2*alpha*B)*np.sin(alpha_a)*np.sinh(b*B)

    #Loop through data, record where RHS = 1, -1
    #a2 = a.copy()
    RHS2 = RHS.copy()
    #a2 = a2[a>=0]
    RHS2=RHS2[alpha_a>=0]

    intr= abs(RHS2) - 1

    intr =intr<0
    intr =[int(inr) for inr in intr]
    intr = np.diff(intr)

    inds = intr.nonzero()[0]


    #bx = 1 if intr[inds[0]]>0 else -1
    #up_down = 1 if RHS[inds[0]]>0 else -1
    n = int(len(alpha_a)/2)
    for ind, n_ind in zip(inds[:-1:2], inds[1::2]):
        w = alpha_a[n_ind+n] - alpha_a[ind+n]
        h = RHS[n_ind+n] - RHS[ind+n]
        rect = plt.Rectangle((alpha_a[ind+n],RHS[ind+n]), w, h, facecolor='gray', alpha=0.3)
        rect2 = plt.Rectangle((-alpha_a[ind+n], RHS[ind+n]), -w, h, facecolor='gray', alpha = 0.3)
        ax2.add_patch(rect)
        ax2.add_patch(rect2)

    if len(inds) % 2 != 0:
        w = alpha_a[-1] - alpha_a[inds[-1]+n]
        h = 2 if RHS[inds[-1]+n] < 0 else -2
        rect = plt.Rectangle((alpha_a[inds[-1]+n],RHS[inds[-1]+n]), w, h, facecolor='gray', alpha=0.3)
        rect2 = plt.Rectangle((-alpha_a[inds[-1]+n], RHS[inds[-1]+n]), -w, h, facecolor='gray', alpha = 0.3)
        ax2.add_patch(rect)
        ax2.add_patch(rect2)

    ax2.plot(alpha_a, RHS)
    ax2.set_ylim([-2, 2])
    ax2.set_xlim([-10, 10])

    ax2.set_yticks(np.arange(-2,3))
    ax2.hlines([-1, 1], -10, 10, colors = 'gray', linestyle=':')
    ax2.set_xlabel(r'$\alpha a$')
    ax2.set_ylabel('RHS')
    ax2.set_title('Figure B')
    plt.show()
    return(RHS, alpha, inds+n)

def Egraphs(RHS, alpha, b, inds):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle
    
    hbar=1.05e-34
    d=1e-9
    m_e=9.11e-31
    a = d-b
    
    fig3 = plt.figure(figsize =(8,6))
    fig3.suptitle('Figure C')
    ax3=fig3.add_subplot(1,2,1)
    ax4=fig3.add_subplot(1,2,2)

    #Convert alpha * a to Energy (E) for free electrons
    #This forms the base for the energy spectrum
    Es = hbar**2 * alpha[alpha>=0]**2/(2*m_e)

    #Truncate the RHS values to those that correspond to positive alpha * a 
    RHS_t = RHS[alpha>=0]

    #Using the free electron E, work backwards to find the k values
    ks = np.sqrt(Es * 2 * m_e / (hbar**2))

    #Flip and append the energy spectrum to obtain parabola
    Es_k = Es
    #Es_k = np.append(np.flip(Es_k), Es_k)

    #Repeat with the k-values
    #ks = np.append(-np.flip(ks), ks)

    #Since the RHS data is truncated we need to subtract the length of 
        #the truncated value from the intersection indices
    
    inds_Ek = inds - len(Es)
    
    
    if len(inds_Ek) % 2 !=0:
        inds_Ek=inds_Ek[:-1]
    inds_Ek = inds_Ek[inds_Ek>0][1:]

    #Calculate the bandgaps
    Egs = [Es[EgU] - Es[EgL] for EgU, EgL in zip(inds_Ek[1::2], inds_Ek[::2])]


    E_min = 0
    E_shift = 0 
    #E_max = (hbar**2 * (np.pi/a)**2 / (2*m_e)) 
    E_max = Es[inds_Ek[0]]
    ks_m = np.linspace(0, np.pi/a, 1000)
    for n, ind in zip(range(len(Egs)), inds_Ek[::2]):
        E_max = Es[ind]
        k_plot = ks_m + n*np.pi/a

        E_mult =E_max - E_min

        E_base = (np.cos(ks_m*a+np.pi) + 1)/2
        E_plot = (E_base ) * E_mult + E_shift
        ax3.plot(k_plot,E_plot, color='blue')
        ax3.plot(-k_plot, E_plot, color='blue')
        E_shift = E_max + Egs[n] 
        E_min = E_max + Egs[n]
        #E_max = (hbar**2 * (np.pi/a * (n+2))**2 / (2*m_e)) 


    ymx = E_plot.max()*1.15
    


    ax3.plot(ks, Es, ':', color = 'red', alpha =0.3)
    ax3.plot(-ks, Es, ':', color = 'red', alpha =0.3)
    ax4.plot(RHS_t, Es)
    ax4.set_xlim([-2, 2])
    ax4.vlines([-1, 1], 0, ymx, colors='gray', linestyle = ':')
    hlin = [Es[ind] for ind in inds_Ek[inds_Ek<len(Es)]]
    
    ax4.hlines([hlin], -2, 2, colors='blue', linestyle='-.', alpha=0.6)
    ax3.hlines([hlin], -6*np.pi/a, 6*np.pi/a, colors = 'blue', linestyle ='-.', alpha=0.6)

    #Set block shading
    for h1, h2 in zip(hlin[1::2], hlin[2::2]):
        w = 2
        w2 = 6*np.pi/a*2
        h = h2-h1
        rect = plt.Rectangle((-1, h1), w, h, facecolor='gray', alpha=0.3)
        rect2 = plt.Rectangle((-6*np.pi/a, h1), w2, h, facecolor='gray', alpha = 0.3)
        ax4.add_patch(rect)
        ax3.add_patch(rect2)
    
    ax3.set_xlim([-6*np.pi/a, 6*np.pi/a])
    ax3.set_ylim([0, ymx])
    ax4.set_ylim([0, ymx])
    #ax3.vlines([np.pi/a, np.pi*2/a], 0, ymx)
    ax3.set_yticks([])
    ax4.set_yticks([])
    ax3.set_title('E vs. k (Extended Scheme)')
    ax4.set_title('E vs. RHS')
    ax3.set_xticks(np.arange(-5*np.pi/a, 6*np.pi/a, np.pi/a))
    ax3.set(xlabel = 'k', ylabel='E')
    ax4.set(xlabel = 'RHS', ylabel = 'E')
    xtck_lab = [r'$\frac{-5\pi}{a}$',r'$\frac{-4\pi}{a}$', 
                r'$\frac{-3\pi}{a}$', r'$\frac{-2\pi}{a}$', r'$\frac{-\pi}{a}$', '0',
                r'$\frac{\pi}{a}$', r'$\frac{2\pi}{a}$', r'$\frac{3\pi}{a}$', 
                r'$\frac{4\pi}{a}$', r'$\frac{5\pi}{a}$']
    ax3.set_xticklabels(xtck_lab)

    plt.show()