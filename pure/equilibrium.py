from ..helpers.equationsOfState import peng_robinson, redlich_kwong_soave
from ..helpers.alfaFunctions import pr78
from ..helpers.stateFunctionsHelpers import A_fun, B_fun, getCubicCoefficients
from ..solvers.cubicSolver import cubic_solver
from numpy import log, exp


def solve_eos(t,p,tc,pc,acentrico,method='pr',alfa='pr78',diagram=False):
    R=83.14
    
    # Method selection
    if(method == 'pr'):
        u,w,omega_a,omega_b = peng_robinson()
    elif(method=='rks'):
        u,w,omega_a,omega_b = redlich_kwong_soave()
    else:
        return 'Method: '+ method+ ' does not exist, define an allowed method'
    #print(u,w,omega_a,omega_b)
    
    # Alpha funciton selection
    if(alfa == 'pr78'):
        alfa = pr78(t,tc,acentrico)
    else:
        return 'Alpha: '+ alfa+ ' does not exist, define an allowed alfa'
    #print(alfa)
    B = B_fun(t,p,tc,pc,omega_b)
    #print(B)
    A = A_fun(t,p,tc,pc,acentrico,omega_a,alfa)
    #print(A)
    coefficients = getCubicCoefficients(A,B,u,w)
    #print(coefficients)
    x= cubic_solver(coefficients,diagram,B)
   
    if(diagram):
        return x
    
    if(isinstance(x,tuple)):
        z_liq = x[0]
        z_vap = x[1]
        
    else:
        z_liq=x
        z_vap=x
        
    
    ln_liq_fugacity = (z_liq-1)-log(z_liq-B) + A/B * 1/(u-w)*log((z_liq+w*B)/(z_liq+u*B))
    ln_vap_fugacity = (z_vap-1)-log(z_vap-B) + A/B * 1/(u-w)*log((z_vap+w*B)/(z_vap+u*B))
    
    liq_fugacity = exp(ln_liq_fugacity)
    vap_fugacity = exp(ln_vap_fugacity)
    
    return (liq_fugacity, vap_fugacity)