from ...helpers import eos
from ...helpers import alfaFunctions
from ...helpers.eosHelpers import A_fun, B_fun, getCubicCoefficients, getMixFugacity
from ...solvers.cubicSolver import cubic_solver
from ...helpers import temperatureCorrelations as tempCorr
from ...helpers import mixing_rules

from numpy import log, exp, sqrt,absolute, array,sum
from scipy.optimize import fsolve, newton, root
from scipy.integrate import quad


def solve_eos(t,p,tc,pc,acentric,liq_compositions,vap_compositions,kij,method='pr',alfa_function='alfa_peng_robinson',mixing_rule='van_der_waals',diagram=False,properties=False,heat_capacity=None):
    # Vectorization
    tc = array(tc)
    pc= array(pc)
    acentric = array(acentric)
    liq_compositions=array(liq_compositions)
    vap_compositions = array(vap_compositions)
    kij = array(kij)
    
    # Method selection
    eos_fun = eos.selector(method)
    u,w,omega_a,omega_b,L = eos_fun()

    # Alfa function selection    
    alfa_fun = alfaFunctions.selector(alfa_function)
    alfa= alfa_fun(t,tc,acentric)
    
    Ai = A_fun(t,p,tc,pc,acentric,omega_a,alfa)
    Bi = B_fun(t,p,tc,pc,omega_b)

    # Mixing rules
    mixing_rule_used = mixing_rules.selector(mixing_rule)
    A_liq,B_liq,A_i_liq,Aij_liq,dAdT_liq = mixing_rule_used(liq_compositions,tc,acentric,kij,Ai,Bi,alfa,alfa_fun,t)
    A_vap,B_vap,A_i_vap,Aij_vap,dAdT_vap = mixing_rule_used(vap_compositions,tc,acentric,kij,Ai,Bi,alfa,alfa_fun,t)


    coefficients_liq = getCubicCoefficients(A_liq,B_liq,u,w)
    coefficients_vap = getCubicCoefficients(A_vap,B_vap,u,w)
   
    z_liq= cubic_solver(coefficients_liq,diagram,B_liq)
    z_vap = cubic_solver(coefficients_vap,diagram,B_vap)

    #z,A,B,A_i,Bi,L,compositions
    liq_fugacity = getMixFugacity(z_liq,A_liq,B_liq,A_i_liq,Bi,L,liq_compositions,p)
    vap_fugacity = getMixFugacity(z_vap,A_vap,B_vap,A_i_vap,Bi,L,vap_compositions,p)
    
    return (liq_fugacity,vap_fugacity)
