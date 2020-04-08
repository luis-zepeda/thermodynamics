from matplotlib.pyplot import scatter, axhline, xlabel, ylabel, legend
from numpy import linspace, log10
from .equilibrium import solve_eos

def pressure_volume(T,tc,pc,acentric):
    volumes = []
    pres = linspace(1,pc)
    R=83.14
    pressures=[]
    for p in pres:
        x=solve_eos(T,p,tc,pc,acentric,diagram=True)
        if(isinstance(x,tuple)):
            volumes.append(x[0]*R*T/p)
            volumes.append(x[1]*R*T/p)
            volumes.append(x[2]*R*T/p)
            pressures.append(p)
            pressures.append(p)
            pressures.append(p)
        else:
            volumes.append(x*R*T/p)
            pressures.append(p)

    scatter(log10(volumes),pressures)
    axhline(pc, color='k', linestyle='--')
    xlabel('log Volume [cm3]')
    ylabel('Pressure [bar]')
    legend(['Critical pressure','PV'])
    
    return 0 