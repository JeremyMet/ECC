import copy ;

class ECC_Arith:


    p = 11 ; # Define field as F_p
    a = 3 ;  # Curve in the form y^2 = x^3+a*x+b mod p
    b = 0 ;

    @staticmethod
    def Inf():
        return ECC_Arith(-1, -1, False) ; # Every value should be positive so negative is a way to distinguish Infinity point.
                                          # "Simulate Constant, not perfect, I don't like this solution that much".

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


    @classmethod
    def set_curve(cls, a, b, p):
        ECC_Arith.a = a ;
        ECC_Arith.b = b ;
        ECC_Arith.p = p ;

    @classmethod
    def is_in_E(cls, x, y):
        return ((y**2 - x**3-ECC_Arith.a*x-ECC_Arith.b)%ECC_Arith.p == 0) or (x == -1 and y == -1) ;

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
        if (self == other and self.y ==0):
            return ECC_Arith.Inf() ;
        diff = (other.x-self.x)%ECC_Arith.p ;
        if (diff == 0):
            inv = ECC_Arith.gcd(2*self.y%ECC_Arith.p, ECC_Arith.p)[1];
            l = ((3*self.x**2+ECC_Arith.a)*inv)%ECC_Arith.p ;
        else:
            inv = ECC_Arith.gcd(diff, ECC_Arith.p)[1] ;
            l = ((other.y-self.y)*inv%ECC_Arith.p) ;
        x_r = (l**2-other.x-self.x)%ECC_Arith.p ;
        y_r = (l*(self.x-x_r)-self.y)%ECC_Arith.p ;
        return ECC_Arith(x_r, y_r) ;

    def __eq__(self, other):
        return (self.x == other.x) and (self.y == other.y) ;

    def __neg__(self):
        return ECC_Arith(self.x, ECC_Arith.p - self.y);

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
        return str((self.x, self.y)) ;






if __name__ == '__main__':




    print("---------------")
    p_0 = ECC_Arith(0, 0) ;
    p_1 = ECC_Arith(3, 5) ;
    p_2 = ECC_Arith(1, 2) ;
    print(p_0+p_0) ;
    print(p_0+p_1) ;
    print(p_0-p_1);
    print(p_1-p_1);
    print(2*p_1) ;
    print(p_2) ;
    print("----------") ;
    for i in range(16):
        print(i*p_2) ;


    R = p_1 ;
    S = p_1 ;
    T = p_2 ;
    U = p_2-p_1+p_2+p_0 ;
    R2 = R+U ;
    S2 = S+U ;
    T2 = T+U ;

    print(R2, S2, T2, U)
    res_U = (R2.y-R2.x)*(S2.x-S2.y)*(T2.x-T2.y)/((U.x-U.y)**3) ;
    print(res_U) ;





