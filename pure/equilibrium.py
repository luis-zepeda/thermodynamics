from ..helpers.equationsOfState import peng_robinson, redlich_kwong_soave
from ..helpers.alfaFunctions import alfa_peng_robinson
from ..helpers.stateFunctionsHelpers import A_fun, B_fun, getCubicCoefficients
from ..solvers.cubicSolver import cubic_solver
from numpy import log, exp, sqrt
from scipy.optimize import fsolve



def solve_eos(t,p,tc,pc,acentric,method='pr',alfa='alfa_peng_robinson',diagram=False,properties=False):
    R=83.14
    
    # Method selection
    if(method == 'pr'):
        u,w,omega_a,omega_b,L = peng_robinson()
    elif(method=='rks'):
        u,w,omega_a,omega_b,L = redlich_kwong_soave()
    else:
        return 'Method: '+ method+ ' does not exist, define an allowed method'
    #print(u,w,omega_a,omega_b)
    
    # Alpha funciton selection
    if(alfa == 'alfa_peng_robinson'):
        alfa = alfa_peng_robinson(t,tc,acentric)
    else:
        return 'Alpha: '+ alfa+ ' does not exist, define an allowed alfa'
    
    B = B_fun(t,p,tc,pc,omega_b)
    
    A = A_fun(t,p,tc,pc,acentric,omega_a,alfa)
    
    coefficients = getCubicCoefficients(A,B,u,w)
  
    x= cubic_solver(coefficients,diagram,B)
   
    if(diagram):
        return x
    
    if(isinstance(x,tuple)):
        z_liq = x[0]
        z_vap = x[1]
        
    else:
        z_liq=x
        z_vap=x
        
    
        
    ln_liq_fugacity_coef = -log(z_liq-B) + (z_liq-1) + A/B *L(z_liq,B)
    ln_vap_fugacity_coef = -log(z_vap-B)+(z_vap-1) + A/B * L(z_vap,B)
   
    liq_fugacity = exp(ln_liq_fugacity_coef)*p
    vap_fugacity = exp(ln_vap_fugacity_coef)*p
    
    
    if(properties):
        
 
    return (liq_fugacity, vap_fugacity)



def solve_VLE(t,p,tc,pc,acentric,solving_for='pressure',method='pr',alfa='alfa_peng_robinson'):
    if(solving_for=='pressure'):
        print('Solving for pressure, given temperature')
        print(fsolve(vle_pressure_objective_function,p,args=(t,tc,pc,acentric,method,alfa),full_output=1))
    elif(solving_for=='temperature'):
        print('Solving for temperature, given pressure')
        print(fsolve(vle_temperature_objective_function,t,args=(p,tc,pc,acentric,method,alfa),full_output=1))
        
def vle_pressure_objective_function(p,t,tc,pc,acentric,method='pr',alfa='alfa_peng_robinson'):
    liq_fugacity, vap_fugacity = solve_eos(t,p,tc,pc,acentric,method,diagram=False)
    return liq_fugacity-vap_fugacity

def vle_temperature_objective_function(t,p,tc,pc,acentric,method='pr',alfa='alfa_peng_robinson'):
    liq_fugacity, vap_fugacity = solve_eos(t,p,tc,pc,acentric,method,diagram=False)
    return liq_fugacity-vap_fugacity