from numpy import sqrt, arccos, absolute, cos, pi, log, exp

def cubic_solver_not_used(coefficients,diagram,B):
    alfa,beta,gamma = coefficients
    alfa= -1*alfa
    p=(3*beta-alfa**2)/3
    q=(27*gamma-9*alfa*beta+2*alfa**3)/27
    R=(p/3)**3+(q/2)**2

    if(R>0):
        first = -q*0.5+sqrt(R)
        second = -q*0.5-sqrt(R)
        
        if(second >0 and first >0 ):
            x = first**(1/3)+second**(1/3)-alfa/3
        elif(second>0 and first<0):
            x = -(-1*first)**(1/3) + second**(1/3)-alfa/3
        elif(second<0 and first>0):
            x = first**(1/3) -(-1*second)**(1/3)-alfa/3
        else:
            x = -(-1*first)**(1/3)-(-1*second)**(1/3)-alfa/3
        return x
        
    elif(R<0):
        theta = arccos(((-27/4)*(q**2/p**3))**(0.5))
        x1=-2*(absolute(q)/q)*(-p/3)**0.5 * cos(theta/3)-alfa/3
        x2=-2*(absolute(q)/q)*(-p/3)**0.5 * cos(theta/3+(2*pi/3))-alfa/3
        x3=-2*(absolute(q)/q)*(-p/3)**0.5 * cos(theta/3+(4*pi/3))-alfa/3
        x = [x1,x2,x3]
        return (min(x),max(x))
    else:
        raise Exception("cubic solver error")
    

def cubic_solver(coefficients,diagram,B):
    alfa,beta,gamma = coefficients
    C = 3*beta - alfa**2
    D = -alfa**3 + 4.5*alfa*beta-13.5*gamma
    Q = C **3 + D **2

    if(Q<= 0):
        theta = arccos(-D/(sqrt(-C**3)))
        
        z_liq = (alfa + 2* sqrt(-C)*cos(theta/3+(2*pi/3)))/3
        z_vap = (alfa + 2* sqrt(-C)*cos(theta/3))/3
        
        
        if (diagram):
            z_extra =  (alfa + 2* sqrt(-C)*cos(theta/3+(4*pi/3)))/3
            return(z_vap,z_liq,z_extra)
        
        if(z_liq < B):
            z_liq = (alfa + 2* sqrt(-C)*cos(theta/3))/3
        
        return (z_liq,z_vap)
        
    elif(Q>0):
        first = -D+sqrt(Q)
        second = -D-sqrt(Q)
        aux1=False
        aux2=False
        if(first<0):
            first *= -1
            aux1 = True
        if(second<0):
            second *= -1
            aux2 = True
            
            
        z = (alfa+(-1 if(aux1)else 1)*(first)**(1/3)+(-1 if(aux2) else 1)*(second)**(1/3))/3
        
        return z
    
    else:
        raise Exception("cubic solver error")
    

def validate_vap_solution(A, B, A_0, rho_0, z, u, w):
    rho = B/z
    F_vap = F(rho, A, B, u, w)
    first_vap_condition = (A/B) < A_0
    second_vap_condition = rho < rho_0 and F_vap > 0.1
    is_z_vap_valid = first_vap_condition or second_vap_condition

    if(not is_z_vap_valid):
        rho_1 = 0.1
        F_1 = (rho_1/(1-rho_1)) - (((A/B)*rho_1**2) / (1+u*rho_1+w*rho_1**2))
        F_2 = (F_vap*((rho_0-rho_1)/2)-F_1) / (F_1**2*((rho_0-rho_1)/2)**2)
        F_3 = 2*F_2*rho_1+F_vap/(F_1**2)
        F_0 = (1/F_1)+F_3*rho_1-F_2*rho_1**2
        print(F_3**2-4*F_2*(F_0-(1/B)))
        new_rho = (F_3-sqrt(F_3**2-4*F_2*(F_0-(1/B))))/(2*F_2)
        z = B/new_rho
    return z

def validate_liq_solution(A, B, A_0, rho_0, z, u, w):
    rho = B/z
    F_liq = F(rho, A, B, u, w)
    first_liq_condition = rho > rho_0
    second_liq_condition = F_liq > 0.1
    is_z_liq_valid = first_liq_condition and second_liq_condition
    
    if(not is_z_liq_valid):
        rho_1 = 0.8
        F_1 = (rho_1/(1-rho_1)) - (((A/B)*rho_1**2) / (1+u*rho_1+w*rho_1**2))
        F_2 = (rho_1-0.7*rho_0)*F_liq
        F_0 = F_1-F_2*log(rho_1-0.7*rho_0)
        new_rho = exp((B-F_0)/F_2)+0.7*rho_0
        z = B/new_rho
        B_0 = (new_rho/(1-new_rho)) - ((A*new_rho**2)/(B*(1+u*new_rho+w*new_rho**2)))
        return (z, B_0)
    
    return z

def F(rho, A, B, u, w):
    return (1/(1-rho)**2)-(((A/B)*rho*(2+u*rho))/(1+u*rho+w*rho**2)**2)
