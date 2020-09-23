def potGraph(a, Etot, V0, ax1):
    
    import numpy as np
    import matplotlib.pyplot as plt
    xtk = np.linspace(0,a*8,9)
    ax1.set_xticks(xtk)
    xlab = [0]+['{}a'.format(n) for n in range(1,9)]
    ax1.set_xticklabels(xlab)

    ytk = [-Etot, -V0, 0, V0, Etot]
    ylab = [r'$-E_{tot}$', r'$-V_0$', '0', r'$V_0$', r'$E_{tot}$']
    ax1.hlines([-Etot, Etot], 0, 8*a, linestyles ="dashed", colors ="gray", alpha=0.5)
    ax1.set_yticks(ytk)
    ax1.set_yticklabels(ylab)
    ax1.set_title('Figure A: Potential Profile - NFE')
    ax1.set_xlabel('x')
    ax1.set_ylabel('V(x)')
    ax1.set_xlim([0, 8*a])
    

def exScheme(ax3, ax4, VG):
    
    def quicksmooth(data):
        n=5 #degree of smoothing
        smoothed=[]

        for i in range(2, len(data) - 2):
            smoothed = np.append(smoothed, sum(data[i-2:i+3])/n)

        smoothed = np.append(data[0:2], smoothed)
        smoothed = np.append(smoothed, data[-2:])
        return smoothed

    import numpy as np
    import matplotlib.pyplot as plt
    hbar=1.05e-34
    a=1e-10
    m_e=9.11e-31
    rep=100
    k_range = np.linspace(0,4*np.pi/a, 9+8*rep)
    ind = (rep+1)*3
    Num = 3

    Eks = hbar**2*k_range**2/(2*m_e)
    colors = ['blue', 'red', 'green', 'magenta']
    marks = ['-', '-.', ':', '--']
    t = 0
    k=k_range.copy()
    for n in range(Num+1):
        
        lwr = (rep+1)*4*n
        lwr = 0 
        
        mid = rep+1
        
        uppr = (rep+1)*2

        nxt_uppr = (rep+1)*(4*n+4)
        nxt_uppr = (rep+1)*4
        #Reduced scheme
        for m in range(Num+1):
            k_red = k[lwr:uppr+1] + n*2*np.pi/a
            k_r = -k_red + n*4*np.pi/a
            Ek = Eks[lwr:uppr+1]
            Ek_m = Eks[lwr:uppr+1]
            Ek_G = np.flip(Eks[uppr:nxt_uppr+1])
            E_new = Ek.copy()

            E_new=1/2*(Ek_m+Ek_G)-np.sqrt(1/4*(Ek_m - Ek_G)**2 + VG**2)
            #
            if m > 0:
                Ek_G = np.flip(Eks[prev_lwr:lwr+1])
                E_new2 =1/2*(Ek_m+Ek_G)+np.sqrt(1/4*(Ek_m - Ek_G)**2 + VG**2)
                E_new = np.append(E_new2[:rep+1], E_new[rep+1:])
                E_new[30:-30] = quicksmooth(E_new[30:-30])


            ax3.plot(k_red, E_new, color = colors[n], linestyle = marks[m])
            ax3.plot(-k_red, E_new, color = colors[n], linestyle = marks[m])
            ax3.plot(k_r, E_new, color = colors[n], linestyle = marks[m])
            ax3.plot(-k_r, E_new, color = colors[n], linestyle = marks[m])
            if n==0:
                ax4.plot(k_red - t*2*np.pi/a, E_new,color = colors[n], linestyle = marks[m])
                ax4.plot(-k_red + t*2*np.pi/a, E_new, color = colors[n], linestyle = marks[m])
                t=t+1 if (m % 2 == 0) else t


            prev_lwr = lwr
            prev_uppr = uppr
            lwr = uppr
            uppr = uppr + (rep+1)*2
            mid = mid + (rep+1)*2
            nxt_uppr = nxt_uppr + (rep+1)*2


def graphSetUp(ax3, ax4):
    import numpy as np
    import matplotlib.pyplot as plt
    a=1e-10
    Num=3
    xtck = [-np.pi/a, np.pi/a]
    ax3.set_xticks(xtck)       
    xtck_lab = [r'$\frac{-\pi}{a}$', r'$\frac{\pi}{a}$']

    ax3.set_xticklabels(xtck_lab)

    ax3.set(title='Extended Zone Scheme')
    xtck = np.append(np.linspace(-(Num)*2*np.pi/a, (Num)*2*np.pi/a, Num*2+1), np.array([-np.pi/a, np.pi/a]))

    ax3.set_xticks(xtck)
    xtck_lab = ['G{}'.format(n) for n in range(-(Num), Num+1)] + [r'$\frac{-\pi}{a}$', r'$\frac{\pi}{a}$']

    ax3.set_yticks([])
    ax3.set_yticklabels([])
    ax3.set_xticklabels(xtck_lab)
    ax4.set_yticks([])
    ax4.set_yticklabels([])
    ax4.set_xticks([-np.pi/a, np.pi/a])
    ax4.set_xticklabels([r'$\frac{-\pi}{a}$', r'$\frac{\pi}{a}$'])
    ax4.set(title='Reduced Zone Scheme')
    ax3.set_xlabel('k')
    ax3.set_ylabel('E(k)')
    ax4.set_xlabel('k')
    ax4.set_ylabel('E(k)')

    
    
    ylims = [-4.85e-18, 9.82e-17]
    ax3.set_ylim(ylims)
    ax4.set_ylim(ylims)
    xlims = [-12*np.pi/a, 12*np.pi/a]
    ax4.set_xlim(xlims)
    ax3.set_xlim(xlims)
    ax3.vlines(np.pi/a, ylims[0],ylims[1], linestyles ="dashed", colors ="gray", alpha=0.5)
    ax3.vlines(-np.pi/a, ylims[0],ylims[1], linestyles ="dashed", colors ="gray", alpha=0.5)
    ax4.vlines(np.pi/a, ylims[0],ylims[1], linestyles ="dashed", colors ="gray", alpha=0.5)
    ax4.vlines(-np.pi/a, ylims[0],ylims[1], linestyles ="dashed", colors ="gray", alpha=0.5)
    plt.show()    


        
