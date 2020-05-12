from np import array,sum, sqrt

def van_der_waals(compositions,binary_interactions,a,b):
    if(len(compositions) != len(binary_interactions) or len(a) != len(b)):
        return 'Compositions and binary interactions must be same dimensional'
    
    compositions= array(compositions)
    binary_interactions = array(binary_interactions)
    a=array(a)
    b=array(b)
    #Yt = (1.0/(M*N)) * sum([Y[i][j] for i in range(M) for j in range(N)])
    b = sum(compositions*b)
    a = sum(compositions[i]*compositions[j]*sqrt(a[i]*a[j])*(1-binary_interactions[i]) for i in range(len(compositions)) for j in range(len(compositions)))
    return (a,b)
