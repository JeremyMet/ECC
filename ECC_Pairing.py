import math ;
from ECC_Arith import ECC_Arith ;
from Gp import Gp ;
from p_2 import p_2 ;
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
            print("R => "+str(R)+str(" vs ")+str(P)) ;
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
    def Pairing(P, ordP, Q, ordQ):
        P2 = P+P ;
        Q2 = Q+Q ;
        DivP = [(P2,1), (P,-1)] ;
        DivQ = [(Q2, 1), (Q, -1)];
        num = ECC_Pairing.Miller(P, DivQ, ordP);
        den = ECC_Pairing.Miller(Q, DivP, ordQ);
        den = ECC_Arith.gcd(den, ECC_Arith.p)[1];
        return (num * den % ECC_Arith.p);



if __name__ == "__main__":

    ECC_Arith.set_curve(p_2, p_2(22), p_2(0)) ;

    P = ECC_Arith(p_2(2), p_2(11)) ;
    Q = ECC_Arith(p_2(21, 0), p_2(0, 12))

    R = ECC_Arith(p_2(0, 17), p_2(21, 2))
    S = ECC_Arith(p_2(18, 10), p_2(13, 13))

    P2 = P+R ;
    Q2 = Q+S ;

    ordP = P.get_order();
    ordQ = Q.get_order();

    DivP = [(P2, 1), (R, -1)] ;
    DivQ = [(Q2, 1), (S, -1)];

    num = ECC_Pairing.Miller(P, DivQ, ordP);
    den = ECC_Pairing.Miller(Q, DivP, ordQ);
    w = num*~den ;

    print(w) ;
    print(w**3) ;

    # ECC_Arith.set_curve(Gp, Gp(17), Gp(6)) ;
    # P = ECC_Arith(Gp(10), Gp(7)) ;
    # Q = ECC_Arith(Gp(16), Gp(2)) ;
    # print("<P>") ;
    # for i in range(P.get_order()):
    #     print(i*P) ;
    # print("<Q>");
    # for i in range(Q.get_order()):
    #     print(i*Q) ;
    #
    # DivQ = [(Q, 1)];
    # Eval = ECC_Pairing.Miller(P, DivQ, 5) ;
    #
    # x = Q.x ;
    # y = Q.y ;
    # Cp = (x+Gp(22))*y+Gp(5)*x**2+Gp(3)*x+Gp(5) ;
    # # Cp = (y+Gp(2)*x+Gp(19))*~(x+Gp(16)) ;
    # # Cp = (Gp(3)*y+x**2+Gp(9)*x+Gp(19))*~(x+Gp(16)) ;
    # # Cp = ((x+Gp(22))*y+Gp(5)*x**2+Gp(3)*x+Gp(5))*~(x+Gp(13)) ;
    #
    # print(Eval) ;
    # print(Cp) ;