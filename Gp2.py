from Gp import Gp ;
import copy ;
import random ;
from Tonelli_Shanks import * ;

class Gp2(object):

    def conjugate(self):
        return Gp2(self.x, -self.y)

    def norm(self):
        return self.x ** 2 + self.y ** 2 ;

    def __init__(self, x=Gp(0), y=Gp(0)):
        self.x = x;
        self.y = y;

    def __add__(self, other):
        x = self.x + other.x ;
        y = self.y + other.y ;
        return Gp2(x, y) ;

    def __neg__(self):
        x = (-self.x) ;
        y = (-self.y)  ;
        return Gp2(x, y);

    def __sub__(self, other):
        return self + (-other);

    def __mul__(self, other):
        # Karatsuba
        z2 = self.x*other.x ;
        z0 = self.y*other.y ;
        z1 = ((self.x+self.y)*(other.x+other.y) - z0 - z2) ;
        x = (z2-z0) ;
        y = z1 ;
        return Gp2(x, y) ;

    def __eq__(self, other):
        return (self.x == other.x) and (self.y == other.y) ;

    def __invert__(self):
        conj = self.conjugate()
        norm = self.norm() ;
        norm = ~norm ;
        norm = Gp2(norm, Gp(0)) ;
        return norm*conj ;

    def __str__(self):
        ret = ""
        if self.x != Gp.zero:
            ret+=str(self.x) ;
        if self.y != Gp.zero:
            if ret:
                ret+="+" ;
            ret+=str(self.y)+"*i"
        if not(ret):
            ret = "0" ;
        return ret ;

    def __truediv__(self, other):
        inv = other.__invert__() ;
        return self*inv ;

    def __pow__(self, power):
        ret = Gp2(Gp(1)) ;
        if power < 0:
            power = -power ;
            buff = self.__invert__() ;
        else:
            buff = copy.deepcopy(self) ;
        while(power>0):
            if power&1:
                ret = ret*buff ;
            buff = buff*buff ;
            power = power >> 1 ;
        return ret ;


    @classmethod
    def random(cls):
        x = Gp(random.randint(0, Gp.p-1)) ;
        y = Gp(random.randint(0, Gp.p-1)) ;
        return Gp2(x, y) ;

    def sqrt(self):
        sqrt = lambda v : Gp(Tonelli_Shanks(v.x, Gp.p)) ;
        a = self.x ;
        b = self.y ;
        val = (a + sqrt(a ** 2 + b ** 2)) / Gp(2) ;
        if not(IsQuadraticResidue(val.x, Gp.p)):
            val = (a - sqrt(a ** 2 + b ** 2)) / Gp(2);
        x = sqrt(val) ;
        val = x**2-a ;
        y = sqrt(val) ;
        if Gp(2)*x*y != b:
            y = -y ;
        ret = Gp2(x, y) ;
        if ret**2 != self:
            return Gp2(Gp(-1), Gp(-1)) ;
        return ret ;


Gp2.zero = Gp2(Gp(0), Gp(0)) ;
Gp2.one = Gp2(Gp(1), Gp(0)) ;
Gp2.i = Gp2(Gp(0), Gp(1)) ;
Gp2.is_quadratic = Gp2(Gp(-1), Gp(-1))
Gp2.infinite = Gp2(Gp(-1), Gp(0)) ;

if __name__ == "__main__":
    # for i in range(100):
    # a = Gp2.random() ;
    a = Gp2(Gp(9), Gp(5)) ;
    b = a**2 ;
    c = b.sqrt() ;
    print(str(a)+" || "+str(b)+" || "+str(c));
    d = Gp2.zero ;

