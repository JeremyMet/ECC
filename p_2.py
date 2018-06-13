

class p_2(): # mod (u^2+1)

    p = 23 ;

    def __init__(self, x, y):
        self.x = x ;
        self.y = y ;

    def __add__(self, other):
        x = (self.x + other.x) % (p_2.p) ;
        y = (self.y + other.y) % (p_2.p);
        return p_2(x, y) ;


    def __sub__(self, other):
        x = (self.x - other.x) % (p_2.p) ;
        y = (self.y - other.y) % (p_2.p);
        return p_2(x, y) ;

    def __mul__(self, other):
        # Karatsuba
        p = p_2.p ;
        z2 = self.x*other.x%p ;
        z0 = self.y*other.y%p ;
        z1 = ((self.x+other.x)*(self.y+other.y) - z0 - z2)%p ;
        x = (z0-z2)%p ;
        y = z1 ;
        return p_2(x, y) ;