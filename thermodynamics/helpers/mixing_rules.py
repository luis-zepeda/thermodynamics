from np import array,sum, sqrt

def selector(mixing_rule):
    if(mixing_rule == 'van_der_waals'):
        return van_der_waals
    else:
        return 'Mixing rule does not exist'

def van_der_waals(compositions,kij,Ai,Bi):
    B = sum(compositions*Bi)
    Aij=[sqrt(A[i]*A[j])*(1-kij[i][j]) for i in range(0,len(A)) for j in range(0,len(A))]
    A = [sum(compositions[i]*compositions[j]*A[i+j])]
