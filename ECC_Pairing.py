import math ;
from ECC_Arith import ECC_Arith ;
import copy ;


## Divisors are in the form [(P_0, n_0), ..., (P_{m-1}, n_{m-1}]

class ECC_Pairing(object):

    @staticmethod
    def line(P, Q, Div):
        ## curve is the Vanishing Set of y^2 = (x^3+a*x+b)
        a = P.a ;
        b = P.b ;
        p = P.p ;
        x_diff = (Q.x - P.x)%p ;
        y_diff = (Q.y - P.y)%p ;
        # TODO IF P == -Q
        if x_diff == 0:
            slope = ((3*P.x**2+a)*ECC_Arith.gcd(2*P.y, p)[1])%p ;
        else:
            slope = (y_diff*ECC_Arith.gcd(x_diff, p)[1])%p ;
        ## A line equation is given by y - A*x - B
        ## where A is the slope and B the intercept
        ## Or l(x,y) = y - A*(x-P.x)+P.y
        ret = 0 ;
        for d in Div:
            pt = d[0] ;
            n = d[1] ;
            tmp = (pt.y - slope*(pt.x - P.x) + P.y)%p ;
            if n < 0:
                n = -n ;
                tmp = ECC_Arith.gcd(tmp, p)[1] ;
            tmp = pow(tmp, n, p) ;
            ret += tmp ;
        return ret ;

    @staticmethod
    def vertical(P, Div):
        # v(x,y) = x-P.x
        p = P.p;
        ret = 0 ;
        for d in Div:
            pt = d[0] ;
            n = d[1] ;
            tmp = (pt.x - P.x)%p ;
            if n < 0:
                n = -n ;
                tmp = ECC_Arith.gcd(tmp, p)[1] ;
            tmp = pow(tmp, n, p) ;
            ret += tmp ;
        return ret ;

    @staticmethod
        # Let D be divisor m(P)-m(0) ;
        # where P belongs to E[m] (m-torsion point)
        # Let f be a rational function such D = (f) ;
        # Then, this function returns f(Div) ;
    def Miller(P, Div, m):
        str_scalar = bin(m)[3:] ; # do not consider the most significant bit
        R = copy.deepcopy(P) ;
        f = 1 ;
        p = P.p ;
        for s in str_scalar:
            l_RR = ECC_Pairing.line(R, R, Div) ;
            R = R+R ;
            v_R = ECC_Pairing.vertical(R, Div) ;
            f = (f**2*l_RR*ECC_Arith.gcd(v_R, p)[1])%p ;
            if s == '1':
                l_RP = ECC_Pairing.line(R, P, Div);
                R = R+P ;
                v_R = ECC_Pairing.vertical(R, Div);
                f = (f * l_RP * ECC_Arith.gcd(v_R, p)[1]) % p;
        return f ;

if __name__ == "__main__":
    # ECC_Arith.set_curve(1,1,23) ;
    gen = ECC_Arith(2,5) ; # Order 12
    points = [i*gen for i in range(12)] ;

    print("All the points on E") ;
    for p in points:
        print(p) ;
    p_1 = points[4] ; order_p1 = 3 ;
    p_2 = points[9] ; order_p2 = 4 ;
    p_3 = points[3] ; order_p3 = 4 ;

    print("Chosen Points for Testing Pairing") ;
    print(p_1) ;
    print(p_2) ;
    print(p_3) ;

    p_1_3 = p_1+p_3 ; order_p_1_3 = 12 ;
    print(p_1_3) ;


    print("Values of Pairing")
    Div_1 = [(p_1+p_1, 1), (p_1, -1)] ;
    Div_2 = [(p_2 + p_2, 1), (p_2, -1)];
    Div_1_3 = [(p_1_3 + p_1_3, 1), (p_1_3, -1)];
    Div_3 = [(p_3 + p_3, 1), (p_3, -1)];

    num = ECC_Pairing.Miller(p_2, Div_1, order_p2) ;
    den = ECC_Arith.gcd(ECC_Pairing.Miller(p_1, Div_2, order_p1), ECC_Arith.p)[1] ;
    weil = (num*den)%ECC_Arith.p ;
    print(weil) ;

    num = ECC_Pairing.Miller(p_2, Div_1_3, order_p2) ;
    den = ECC_Arith.gcd(ECC_Pairing.Miller(p_1_3, Div_2, order_p_1_3), ECC_Arith.p)[1] ;
    weil = (num*den)%ECC_Arith.p ;
    print(weil);

    num = ECC_Pairing.Miller(p_2, Div_3, order_p2);
    den = ECC_Arith.gcd(ECC_Pairing.Miller(p_3, Div_2, order_p3), ECC_Arith.p)[1];
    weil = (num * den) % ECC_Arith.p;
    print(weil);



