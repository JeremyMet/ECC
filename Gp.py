import random ;

class Gp(object):

    p = 23 ; # Field Characteristic

    @classmethod
    def gcd(cls, a, b):
        rev = False ;
        if a < b:
            a, b = b, a ; rev = True ;
        r0 = a ; u0 = 1 ; v0 = 0 ;
        r1 = b ; u1 = 0 ; v1 = 1 ;
        while(r0!=0):
            c = r1//r0 ;
            r1 = r1 - r0*c ;
            u1 = u1 - u0*c ;
            v1 = v1 - v0*c ;
            r0, r1 = r1, r0;
            u0, u1 = u1, u0;
            v0, v1 = v1, v0;
        if rev:
            u1, v1 = v1, u1 ;
        return (r1, u1, v1) ;


    def __init__(self, val = 0):
        self.x = val ;

    def __add__(self, other):
        x = (self.x+other.x)%Gp.p ;
        return Gp(x) ;

    def __neg__(self):
        x = (-self.x)%Gp.p ;
        return Gp(x) ;

    def __sub__(self, other):
        return self+(-other) ;

    def __mul__(self, other):
        x = (self.x*other.x)%Gp.p ;
        return Gp(x) ;

    def __invert__(self):
        x = (Gp.gcd(self.x, Gp.p)[1])%Gp.p ;
        return Gp(x) ;

    def __truediv__(self, other):
        return self*(~other) ;

    def __str__(self):
        return str(self.x)

    def __pow__(self, power):
        x = self.x
        if power < 0:
            power = -power ;
            x = (Gp.gcd(x, Gp.p)[1])%Gp.p ;
        return Gp(pow(x, power, Gp.p)) ;

    def __eq__(self, other):
        return self.x == other.x ;

    def __ne__(self, other):
        return not(self == other) ;

    @classmethod
    def random(cls):
        return Gp(random.randint(0, Gp.p)) ;

## Class Constants
Gp.zero = Gp(0) ;
Gp.one = Gp(1) ;

if __name__ == "__main__":
    A = Gp(5) ;
    print(A**-3) ; 
