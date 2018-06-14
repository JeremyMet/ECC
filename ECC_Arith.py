import copy ;
import random ;
from Gp import Gp ;
from Gp2 import Gp2 ;


class ECC_Arith:


    field = Gp ; # Define field as F_p
    a = Gp(3) ;  # Curve in the form y^2 = x^3+a*x+b where a and b belong to "field"
    b = Gp(0) ;

    @staticmethod
    def Inf():
        tmp = -ECC_Arith.field.infinite ;
        return ECC_Arith(tmp, tmp, False) ; # Every value should be positive so negative is a way to distinguish Infinity point.
                                            # "Simulate Constant, not perfect, I don't like this solution that much".


    @classmethod
    def set_curve(cls, field, a, b):
        ECC_Arith.field = field;
        ECC_Arith.a = a ;
        ECC_Arith.b = b ;

    @classmethod
    def is_in_E(cls, x, y):
        return ((y**2 - x**3-ECC_Arith.a*x-ECC_Arith.b) == ECC_Arith.field.zero) or (x == ECC_Arith.field.infinite and y == ECC_Arith.field.infinite) ;

    def __init__(self, x, y, check = True):
        if check and not(ECC_Arith.is_in_E(x, y)):
            raise("Point is not of the curve") ;
        self.x = x ;
        self.y = y ;


    def __add__(self, other):
        if self == ECC_Arith.Inf():
            return copy.deepcopy(other) ;
        if other == ECC_Arith.Inf():
            return copy.deepcopy(self) ;
        if self == - other:
            return ECC_Arith.Inf() ;
        if (self == other and self.y == ECC_Arith.field.zero):
            return ECC_Arith.Inf() ;
        diff = (other.x-self.x) ;
        if (diff == ECC_Arith.field.zero):
            inv = (self.y+self.y).__invert__();
            tmp = self.x**2 ;
            l = (tmp+tmp+tmp+ECC_Arith.a)*inv ;
        else:
            inv = diff.__invert__() ;
            l = ((other.y-self.y)*inv) ;
        x_r = (l**2-other.x-self.x) ;
        y_r = (l*(self.x-x_r)-self.y) ;
        return ECC_Arith(x_r, y_r) ;

    def __eq__(self, other):
        return (self.x == other.x) and (self.y == other.y) ;

    def __neg__(self):
        return ECC_Arith(self.x, -self.y);

    def __sub__(self, other):
        return self+(-other) ;

    def __rmul__(self, n):
        ret = ECC_Arith.Inf() ;
        buff = copy.deepcopy(self) ;
        if (n<0):
            n = -n ;
            buff = -buff ;
        while(n > 0):
            if (n&1):
                if ret:
                    ret = ret+buff ;
                else:
                    ret = copy.deepcopy(buff) ;
            n = n >> 1 ;
            buff = buff+buff ;
        return ret ;

    def __str__(self):
        return str((str(self.x), str(self.y))) ;


    @staticmethod
    # For now, only work for Gp for now (no extension).
    def random_point():
        has_square_root = False ;
        while(not(has_square_root)):
            x = ECC_Arith.field.random() ;
            tmp = (x**3+ECC_Arith.a*x+ECC_Arith.b) ;
            y = tmp.sqrt() ;
            has_square_root = (y != ECC_Arith.field.is_quadratic)
        return ECC_Arith(x, y) ;

    # Naive method for now :-(
    def get_order(self):
        tmp = copy.deepcopy(self) ;
        cpt = 1 ;
        while(tmp != ECC_Arith.Inf()):
            tmp = tmp+self ;
            cpt+=1 ;
        return cpt ;

if __name__ == '__main__':
    Gp.set_p(23) ;
    ECC_Arith.set_curve(Gp2, Gp2(Gp(22),Gp(0)), Gp2(Gp(0),Gp(0))) ;
    p_1 = ECC_Arith(Gp2(Gp(21), Gp(0)), Gp2(Gp(0), Gp(12)));
    p_2 = ECC_Arith.random_point() ;
    print(p_1.get_order()) ;
    for i in range(p_1.get_order()):
        print(i*p_1) ;
    print(p_2.get_order()) ;
    for i in range(p_2.get_order()):
        print(i*p_2) ;
