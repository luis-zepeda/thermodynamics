from ..helpers import eos
from ..helpers import alfaFunctions
from ..helpers.eosHelpers import A_fun, B_fun, getCubicCoefficients, dAdT_fun
from ..solvers.cubicSolver import cubic_solver
from ..helpers import temperatureCorrelations as tempCorr

from numpy import log, exp, sqrt,absolute, array,sum
from scipy.optimize import fsolve, newton, root
from scipy.integrate import quad


def solve_eos(t,p,tc,pc,acentric,compositions,kij,method='pr',alfa_function='alfa_peng_robinson',mixing_rule='van_der_waals',diagram=False,properties=False,heat_capacity=None):
    # Vectorization
    tc = array(tc)
    pc= array(pc)
    acentric = array(acentric)
    compositions=array(compositions)
    kij = array(kij)
    
    # Method selection
    eos_fun() = eos.selector(method)
    u,w,omega_a,omega_b,L = eos_fun()

    # Alfa function selection    
    alfa_fun = alfaFunctions.selector(alfa_function)
    alfa= alfa_fun(t,tc,acentric)

    Ai = A_fun(t,p,tc,pc,acentric,omega_a,alfa)
    Bi = B_fun(t,p,tc,pc,omega_b)

    # Mixing rules
    

