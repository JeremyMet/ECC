import math ;
from ECC_Arith import ECC_Arith ;
from Gp import Gp ;
from p_2 import p_2 ;
from Gp2 import Gp2 ;
import copy ;


## Divisors are in the form [(P_0, n_0), ..., (P_{m-1}, n_{m-1}]

class ECC_Pairing(object):

    @staticmethod
    def line(P, Q, Div):
        # Special cases
        if P == -Q:
            return ECC_Pairing.vertical(P, Div) ;
        if P == ECC_Arith.Inf():
            return ECC_Pairing.vertical(Q, Div) ;
        if Q == ECC_Arith.Inf():
            return ECC_Pairing.vertical(P, Div) ;
        ## curve is the Vanishing Set of y^2 = (x^3+a*x+b)
        a = P.a ;
        b = P.b ;
        x_diff = (Q.x - P.x) ;
        y_diff = (Q.y - P.y) ;
        if x_diff == P.field.zero:
            tmp = P.x**2
            slope = ((tmp+tmp+tmp+a) * ~(P.y+P.y)) ;
        else:
            slope = (y_diff * ~x_diff) ;
        ## A line equation is given by y - A*x - B
        ## where A is the slope and B the intercept
        ## Or l(x,y) = y - A*(x-P.x)+P.y
        ret = P.field.one ;
        for d in Div:
            pt = d[0] ;
            n = d[1] ;
            tmp = (pt.y - slope*(pt.x - P.x) - P.y) ;
            if n < 0:
                n = -n ;
                tmp = tmp.__invert__() ;
            tmp = tmp**n ;
            ret *= tmp ;
        return ret ;

    @staticmethod
    def vertical(P, Div):
        # v(x,y) = x-P.x
        if P == ECC_Arith.Inf():
            return P.field.one ;
        ret = P.field.one ;
        for d in Div:
            pt = d[0] ;
            n = d[1] ;
            tmp = (pt.x - P.x) ;
            if n < 0:
                n = -n ;
                tmp = tmp.__invert__() ;
            tmp = tmp**n ;
            ret *= tmp ;
        return ret ;

    @staticmethod
        # Let D be divisor m(P)-m(0) ;
        # where P belongs to E[m] (m-torsion point)
        # Let f be a rational function such D = (f) ;
        # Then, this function returns f(Div) ;
    def Miller(P, Div, m):
        str_scalar = bin(m)[3:] ; # Do not consider the most significant bit
        R = copy.deepcopy(P) ;
        f = P.field.one ;
        for s in str_scalar:
            l_RR = ECC_Pairing.line(R, R, Div) ;
            R = R+R ;
            v_R = ECC_Pairing.vertical(R, Div) ;
            f = ((f ** 2) * l_RR * ~v_R) ;
            if s == '1':
                l_RP = ECC_Pairing.line(R, P, Div);
                R = R+P ;
                v_R = ECC_Pairing.vertical(R, Div);
                f = (f * l_RP * ~v_R)
        return f ;

    @staticmethod
    def Pairing(P, Q):
        num, den = P.field.zero, P.field.zero ;
        while(num == P.field.zero or den == P.field.zero):
            R = ECC_Arith.random_point() ;
            S = ECC_Arith.random_point() ;
            P2 = P + R;
            Q2 = Q + S;
            ordP = P.get_order();
            ordQ = Q.get_order();
            DivP = [(P2, 1), (R, -1)];
            DivQ = [(Q2, 1), (S, -1)];
            num = ECC_Pairing.Miller(P, DivQ, ordP) * (ECC_Pairing.vertical(P2, DivQ) / ECC_Pairing.line(P, R, DivQ)) ** ordP;
            den = ECC_Pairing.Miller(Q, DivP, ordQ) * (ECC_Pairing.vertical(Q2, DivP) / ECC_Pairing.line(Q, S, DivP)) ** ordQ;
        w = num / den;
        return w ;




if __name__ == "__main__":

    Gp.set_p(23) ;
    ECC_Arith.set_curve(Gp2, Gp2(Gp(22),Gp(0)), Gp2(Gp(0), Gp(0))) ;
    P = ECC_Arith(Gp2(Gp(2), Gp(0)), Gp2(Gp(11), Gp(0))) ;
    Q = ECC_Arith(Gp2(Gp(21), Gp(0)), Gp2(Gp(0), Gp(12))) ;

    while(True):
        print(ECC_Pairing.Pairing(P, Q)) ;



