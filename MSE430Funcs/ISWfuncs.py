def potentialProf(xs, L):
    import numpy as np
    Vs = np.zeros(len(xs))
    Vs[(xs<0) | (xs>L)] = 1e6
    return(Vs)

def psi_isw(xs, n, L):
    import numpy as np
    #Calculate Psi at each x
    psi_x = np.sqrt(2/L)*np.sin(n*np.pi*xs/L)
    #Set Psi to zero outside of the well
    psi_x[(xs<0) | (xs>L)] = 0
    return(psi_x)

def graphProp(ax2, ax3, ax4, L, Ltot):

    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider, Button, RadioButtons
    
    ax2.set(title = 'Infinite Square Well', ylabel = 'V', yticks = [0, 1e6], yticklabels=['0', r'$\infty$'],
           xticks = [], xticklabels=[])
    ymin, ymx= ax2.get_ylim()
    ax2.text(-Ltot*.9, ymx/2, 'Potential Profile')

    ax3.set(ylabel = r'$\Psi$', yticks=[0], xticks = [],
            xticklabels = [])
    ax4.set(xlabel = 'x', ylabel = r'$\Psi*\Psi$', xticks=[0,L], xticklabels=['0', 'L'])
    



