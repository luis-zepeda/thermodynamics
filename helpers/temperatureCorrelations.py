from numpy import exp

class TemperatureCorrelations:
    
    def equation_selector(self,number):
        if(number == '16'):
            return self.eq_16

        else:
            return None



    def eq_16(self,t,constants):
        a,b,c,d,e = constants
        return (a+exp(b/t)+c+d*t+e*t**2)
