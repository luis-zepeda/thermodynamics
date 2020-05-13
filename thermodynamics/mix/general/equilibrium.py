from ..helpers import eos
from ..helpers import alfaFunctions
from ..helpers.eosHelpers import A_fun, B_fun, getCubicCoefficients, dAdT_fun
from ..solvers.cubicSolver import cubic_solver
from ..helpers import temperatureCorrelations as tempCorr

from numpy import log, exp, sqrt,absolute, array
from scipy.optimize import fsolve, newton, root
from scipy.integrate import quad


def solve_eos(t,p,tc,pc,acentric,compositions,method='pr',alfa_function='alfa_peng_robinson',diagram=False,properties=False,heat_capacity=None):
    # Vectorization
    tc = array(tc)
    pc= array(pc)
    compositions=array(compositions)
    
    
    # Method selection
    u,w,omega_a,omega_b,L = eos.selector(method)

    # Alfa function selection    
    alfa_fun = alfaFunctions.selector(alfa_function)
    alfa= alfa_fun(t,tc,acentric)
    