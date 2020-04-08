from numpy import sqrt
        
    
def peng_robinson():
    u=2
    w=-1
    omega_a = 0.45723553
    omega_b = 0.077796074
    return (u,w,omega_a,omega_b)

def redlich_kwong_soave():
    u=1
    w=0
    omega_a = 0.42748023
    omega_b = 0.08664035
    return (u,w,omega_a,omega_b)
