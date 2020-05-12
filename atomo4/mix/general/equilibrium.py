from ..helpers import eos
from ..helpers import alfaFunctions
from ..helpers.eosHelpers import A_fun, B_fun, getCubicCoefficients, dAdT_fun
from ..solvers.cubicSolver import cubic_solver
from ..helpers import temperatureCorrelations as tempCorr

from numpy import log, exp, sqrt,absolute
from scipy.optimize import fsolve, newton, root
from scipy.integrate import quad


def solve_eos():
    pass