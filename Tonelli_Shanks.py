


def Tonelli_Shanks(n, p):

    # Find is n is a quadratic residue.
    if pow(n, (p-1) >> 1, p) == p-1:
        return 0 ;

    # Find Q and S such p-1 = Q*2^S.
    S = 0 ;
    Q = p-1 ;
    while(Q & 1 == 0):
        S = S+1 ;
        Q = Q >> 1 ;

    # Find z, a non-quadratic residue of Z/pZ.
    z = 1 ;
    while (pow(z, (p-1) >> 1, p) != p-1):
        z+=1 ;

    # Loop Initialization.
    M = S ;
    c = pow(z, Q, p) ;
    t = pow(n, Q, p) ;
    R = pow(n, (Q+1) >> 1, p) ;

    # Main Loop (...)
    while(t != 0 and t != 1):
        # Find the smallest integer i such that t^(2^i) == 1
        t_tmp = t ;
        for i in range(1, M):
            t_tmp = (t_tmp*t_tmp)%p ;
            if t_tmp == 1:
                break ;
        b = pow(c, 1 >> (M-i-1), p) ;
        M = i ;
        c = (b*b)%p ;
        t = (t*c)%p ;
        R = (R*b)%p ;

    # Post treatment.
    if t == 0:
        ret = 0 ;
    else:
        ret = R ;
    return ret ;

if __name__ == "__main__":
    n = 8 ; p = 17 ;
    r  = Tonelli_Shanks(n, p) ;
    print((n, r)) ;
    print(r**2%p)