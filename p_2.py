
import copy ;
import random ;

class p_2(): # mod (i^2+1)

    p = 23 ;

    @classmethod
    def gcd(cls, a, b):
        rev = False;
        if a < b:
            a, b = b, a;
            rev = True;
        r0 = a;
        u0 = 1;
        v0 = 0;
        r1 = b;
        u1 = 0;
        v1 = 1;
        while (r0 != 0):
            c = r1 // r0;
            r1 = r1 - r0 * c;
            u1 = u1 - u0 * c;
            v1 = v1 - v0 * c;
            r0, r1 = r1, r0;
            u0, u1 = u1, u0;
            v0, v1 = v1, v0;
        if rev:
            u1, v1 = v1, u1;
        return (r1, u1, v1);


    def conjugate(self):
        return p_2(self.x, -self.y%p_2.p)

    def norm(self):
        return (self.x**2+self.y**2)%p_2.p  ;

    def __init__(self, x=0, y=0):
        self.x = x ;
        self.y = y ;

    def __add__(self, other):
        if (type(other) == int):
            x = (self.x+other)%p_2.p ;
            y = self.y ;
        else:
            x = (self.x + other.x) % (p_2.p) ;
            y = (self.y + other.y) % (p_2.p);
        return p_2(x, y) ;

    def __neg__(self):
        x = (-self.x)%p_2.p ;
        y = (-self.y)%p_2.p ;
        return p_2(x, y) ;

    def __sub__(self, other):
        return self+(-other) ;

    def __mul__(self, other):
        if type(other) == int:
            x = self.x * other % p_2.p ;
            y = self.y * other % p_2.p ;
        else:
            # Karatsuba
            p = p_2.p ;
            z2 = self.x*other.x%p ;
            z0 = self.y*other.y%p ;
            z1 = ((self.x+self.y)*(other.x+other.y) - z0 - z2)%p ;
            x = (z2-z0)%p ;
            y = z1 ;
        return p_2(x, y) ;

    def __eq__(self, other):
        if (type(other) == int):
            return (self.y == 0 and self.x == other)
        return ((self.x == other.x) and (self.y == other.y)) ;

    def __invert__(self):
        conj = self.conjugate()
        norm = self.norm() ;
        norm = (p_2.gcd(norm, p_2.p)[1])%p_2.p ;
        norm = p_2(norm, 0) ;
        return norm*conj ;

    def __str__(self):
        ret = ""
        if self.x != 0:
            ret+=str(self.x) ;
        if self.y != 0:
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
        ret = p_2(1, 0) ;
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

    def __mod__(self, modulo):
        x = self.x%modulo ;
        y = self.y%modulo ;
        return p_2(x, y) ;

    def __rsub__(self, other):
        return (-self)+other ;

    def __rmul__(self, other):
        return self*other ;

    @classmethod
    def random(cls):
        x = random.randint(0, p_2.p-1) ;
        y = random.randint(0, p_2.p - 1);
        return p_2(x, y) ;


p_2.zero = p_2(0, 0) ;
p_2.one = p_2(1, 0) ;

if __name__ == "__main__":
    a = p_2(21,8) ;
    print(a*a) ;
    print(a*a*a*a*a) ;
    print(a**5)
    print(type(a)==int) ;

