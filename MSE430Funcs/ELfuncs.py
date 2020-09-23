def ExScheme(Num, k_range, E, ax1):
    import numpy as np
    import matplotlib.pyplot as plt


    a=1e-10
    k = k_range.copy()
    colors = ['blue', 'red', 'green', 'magenta']
    marks = ['-', '-.', ':', '--']

    
    
    for n in range(Num+1):
        k_r = -k + n*4*np.pi/a
        lwr = np.pi/a * n * 2
        uppr = lwr + np.pi/a
        for m in range(Num+1):
            ax1.plot(k[(k<=uppr) & (k>=lwr)],E[(k<=uppr) & (k>=lwr)], color = colors[n], linestyle = marks[m])
            ax1.plot(k_r[(k<=uppr) & (k>=lwr)],E[(k<=uppr) & (k>=lwr)], color = colors[n], linestyle = marks[m])
            ax1.plot(-k[(k<=uppr) & (k>=lwr)],E[(k<=uppr) & (k>=lwr)], color = colors[n], linestyle = marks[m])
            ax1.plot(-k_r[(k<=uppr) & (k>=lwr)],E[(k<=uppr) & (k>=lwr)], color = colors[n], linestyle = marks[m])
            lwr = uppr
            uppr = uppr+np.pi/a
        k = k + 2*np.pi/a
        
    ax1.set(title='Extended Zone Scheme')
    xtck = np.append(np.linspace(-(Num)*2*np.pi/a, (Num)*2*np.pi/a, Num*2+1), np.array([-np.pi/a, np.pi/a]))

    ax1.set_xticks(xtck)
    xtck_lab = ['G{}'.format(n) for n in range(-(Num), Num+1)] + [r'$\frac{-\pi}{a}$', r'$\frac{\pi}{a}$']

    ax1.set_yticks([])
    ax1.set_yticklabels([])
    ax1.set_xticklabels(xtck_lab)
    
    
def ReScheme(Num, k_range, E, ax2):
    import numpy as np
    import matplotlib.pyplot as plt
    a=1e-10
    t=0
    k= k_range.copy()
    colors = ['blue', 'red', 'green', 'magenta']
    marks = ['-', '-.', ':', '--']
    for m in range(Num+1):

        lwr = np.pi/a * m
        uppr = lwr + np.pi/a

        ax2.plot(k[(k<=uppr) & (k>=lwr)]-2*t*np.pi/a,E[(k<=uppr) & (k>=lwr)], color = colors[0], linestyle = marks[m])
        ax2.plot(-k[(k<=uppr) & (k>=lwr)]+2*t*np.pi/a,E[(k<=uppr) & (k>=lwr)], color = colors[0], linestyle = marks[m])
        t=t+1 if (m % 2 == 0) else t
        
def graphSetUp(ax1, ax2):
    import matplotlib.pyplot as plt
    import numpy as np
    a=1e-10
    ylims = ax1.get_ylim()
    ax2.set_ylim([0, ylims[1]])
    ax1.set_ylim([0, ylims[1]])

    ax2.set_xlim(ax1.get_xlim())

    ax2.vlines(np.pi/a, ylims[0],ylims[1], linestyles ="dashed", colors ="gray", alpha=0.5)
    ax2.vlines(-np.pi/a, ylims[0],ylims[1], linestyles ="dashed", colors ="gray", alpha=0.5)
    ax1.vlines(np.pi/a, ylims[0],ylims[1], linestyles ="dashed", colors ="gray", alpha=0.5)
    ax1.vlines(-np.pi/a, ylims[0],ylims[1], linestyles ="dashed", colors ="gray", alpha=0.5)

    ax2.set(title='Reduced Zone Scheme')
    ax2.set_yticks([])
    ax2.set_yticklabels([])
    xtck = [-np.pi/a, np.pi/a]
    ax2.set_xticks(xtck)
    xtck_lab = [r'$\frac{-\pi}{a}$', r'$\frac{\pi}{a}$']

    ax2.set_xticklabels(xtck_lab)

    plt.show()