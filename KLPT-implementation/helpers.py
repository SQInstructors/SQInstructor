from sage.all import vector, IntegralLattice, matrix,ZZ, QQ, matrix, IntegralLattice
from sage.all import GF, sqrt, ZZ, gcd, i, Integer
from sage.all import log
from sage.all import is_prime, vector, polygen, NumberField
from sage.all import floor
import random
# For tests
from sage.all import QuaternionAlgebra, Matrix

###Other helpers
def size(x):
    return log(x,2).numerical_approx(16)

def print_dbg(string, debug = False):
    if debug:
        print(string)

def check_power_of_l(x,l=2):
    pass


#Cornacchia
def Cornacchia_prime(M,d=1):
    if(M==1):
        return(1,0)
    if(M==2):
        if(d==1):
            return(1,1)
        elif(d==2):
            return(0,1)
        else:
            return False
    if(ZZ(M).is_prime()):
        one = GF(M)(-d)
        u = one.sqrt()
        if(u not in GF(M)):
            return False
        a = M
        b = ZZ(u)
        r = M
        while(r**2>=M):
            a = b
            b = r
            r = a%b
        x = r
        y = ZZ(sqrt((M-x**2)//d))
        if(x**2 + d*y**2 == M):
            return(x,y)
        else:
            return(False)

def Corvo_imperiale(M):
    '''
    Extended cornacchia (signature is a joke)
    '''
    bad = Integer(9758835872005600424541432018720696484164386911251445952809873289484557900064682994157437552186699082331048989)
    if(gcd(bad,M)!=1):
        return(False)

    good_primes=[2,5,13,17,29, 37, 41, 53, 61, 73, 89, 97, 101, 109, 113, 137, 149, 157, 173, 181, 193, 197, 229, 233, 241, 257, 269, 277, 281, 293, 313, 317, 337, 349, 353, 373, 389, 397, 401, 409, 421, 433, 449, 457, 461]
    n = Integer(M)
    n_val = [0 for i in range (len(good_primes))]
    for k in range(len(good_primes)):
        p = Integer(good_primes[k])
        n_val[k] = n.valuation(p)
        n = n // (p**n_val[k])

    if(n %4 == 3):
        return(False)
    if not ((ZZ(n)).is_prime() or n == 1):
        return False
    x,y = Cornacchia_prime(n,1)
    c = ZZ[i](x+i*y)
    for k in range(len(good_primes)):
        if(n_val[k]!= 0):
            x_p,y_p = Cornacchia_prime(good_primes[k],1)
            c_p = ZZ[i](x_p+i*y_p)
            c = c*c_p**n_val[k]
    return(abs(c[0]),abs(c[1]))


### Ideal tools
def ideal_inverse(I):
    J = I.conjugate()
    J = J*(1/I.norm())
    return(J)

#to be tested
def ideal_isomorphism(I,J):
    i,j,k = I.quaternion_algebra().gens()
    trans = ideal_inverse(I)*J
    M = trans.basis_matrix()
    M = M.LLL()
    v = M[0]
    isom = v[0]+i*v[1]+j*v[2]+k*v[3]
    assert(I*isom == J)
    return(isom)


def ideal_element(x,I):
    xt = x.coefficient_tuple()
    v = vector([xt[0],xt[1],xt[2],xt[3]])
    M = I.basis_matrix()
    d = M.denominator()
    M = d*M
    M.change_ring(ZZ)
    p = I.quaternion_algebra().ramified_primes()[0]
    G = matrix.diagonal([1,1,p,p])
    L = IntegralLattice(G,M)
    v = d*v
    for w in v:
        if not w.is_integer():
            return(False)
    v.change_ring(ZZ)
    return(v in L)

def lattice_element(x,L):
    xt = x.coefficient_tuple()
    v = vector([xt[0],xt[1],xt[2],xt[3]])
    d = L.denominator()
    M = d*L
    M.change_ring(ZZ)
    G = matrix.diagonal([1,1,1,1])
    L = IntegralLattice(G,M)
    v = d*v
    for w in v:
        if not w.is_integer():
            return(False)
    v.change_ring(ZZ)
    return(v in L)

def ideal_add(I,M):
    i,j,k = I.quaternion_algebra().gens()
    M4 = I.basis_matrix()
    M4 = M4.transpose()
    M5 = M4.augment(M.transpose())
    d5 = M5.denominator()
    M5 = d5*M5
    M5 = M5.transpose()
    M5 = M5.change_ring(ZZ)
    M5 = M5.echelon_form()
    M4 = M5[:4]
    M4 = M4.change_ring(QQ)
    M4 = (1/d5)*M4
    B = [M4[x,0] + M4[x,1]*i+M4[x,2]*j+M4[x,3]*k for x in range(4)]
    J = I.quaternion_algebra().ideal(B)
    assert(ideal_element((M[0,0]+M[0,1]*i+M[0,2]*j+M[0,3]*k),J))
    assert all([ ideal_element(element,J) for element in I.basis()]), 'I not contained in the output'
    return(J)

#Here mostly for tests
#Argument n is the target integer, B the quaternion algebra
def representInteger(n, B, t=1000):
    i,j,k = B.gens()
    p = B.ramified_primes()[0]
    bi = floor((n//p).sqrt())
    if(bi<=0):
        return False
    for it in range(t):
        x = random.randrange(0,bi)
        bj = floor(((n-x*x)//p).sqrt())
        y = random.randrange(0,bj)
        r = Corvo_imperiale(n-p*(x*x+y*y))
        if(r!= False):
            res = r[0]+i*r[1]+j*x+k*y
            assert(res.reduced_norm()==n)
            if(gcd(ZZ(res.reduced_norm()),n*n) == n):
                return(res)
    return(False)

def intersection(M,L2):
    '''
        compute the intersection of two lattices of dimensions 4 and 2
    '''
    M4 = M
    M2 = L2
    d4 = M4.denominator()
    d2 = M2.denominator()
    M4 = d4*d2*M4
    M4.change_ring(ZZ)
    M2 = d4*d2*M2
    M2.change_ring(ZZ)
    M = [[0,0,0,0,0,0,0,0] for i in range(6)]
    M = matrix(ZZ,M)
    for i in range(4):
        for j in range(4):
            M[i,j] = M4[i,j]
            M[4+i%2,j] = M2[i%2,j]
            M[4+i%2,4+j] = M2[i%2,j]
    M = M.echelon_form()
    inter = M[4:6,4:8]
    G = matrix.identity(ZZ,4)
    assert(inter[0] in IntegralLattice(G,M2))
    assert(inter[0] in IntegralLattice(G,M4))
    assert(inter[1] in IntegralLattice(G,M2))
    assert(inter[1] in IntegralLattice(G,M4))
    inter.change_ring(QQ)
    return((1/(d2*d4))*inter)

#Intersection of 2 fractional ideals
def intersection_ideals(I1,I2):
    M1 = I1.basis_matrix().transpose().inverse()
    M2 = I2.basis_matrix().transpose().inverse()
    B = I1.quaternion_algebra()
    J1 = B.ideal(M1)
    J2 = B.ideal(M2)
    J = J1+J2
    M = J.basis_matrix().transpose().inverse()
    I = B.ideal(M)
    return(I)

# Connecting ideal from O1 to O2
def connectingIdeal(O1,O2):
    b1 = O1.basis()
    b2 = O2.basis()
    O1n2 = O1.free_module().intersection(O2.free_module())
    n12 = O1n2.index_in(O1.free_module())
    n = ZZ(n12)
    b = []
    for e1 in b1:
        for e2 in b2:
            b.append(n*e1*e2)
    I = O1.quaternion_algebra().ideal(b)
    return(I)

def idealPullback(Ic,I):
    assert(gcd(Ic.norm(),I.norm())==1)
    IcI = Ic*I
    Ip = Ic*I+Ic.left_order()*I.norm()
    return(Ip)

def idealPushforward(Ic,I):
    assert(gcd(Ic.norm(),I.norm())==1)
    IcI = I.quaternion_algebra().ideal(Ic.free_module().intersection(I.free_module()).basis())
    Ip = (1/Ic.norm())*Ic.conjugate()*IcI
    return(Ip)

#Argument I is the target ideal
def equivalentPrimeIdeal(I, bound=10,margin=False):
    p = I.quaternion_algebra().ramified_primes()[0]
    M = I.basis_matrix()
    M = M.LLL()
    B = I.basis()
    G = I.gram_matrix()
    U = G.LLL_gram().transpose()
    M = [sum(c * g for c, g in zip(row, B)) for row in U]
    n_I = I.norm()
    i,j,k = I.quaternion_algebra().gens()
    for x in range(bound):
        for a in range (0,x):
            for b in range (-x+a,x-a):
                for c in range (-x+a+b,x-a-b):
                    for d in range (-x+a+b+c,x-a-b-c):
                        v = vector([a,b,c,d])
                        el = sum( v[i]*M[i] for i in range(4))
                        n = el.reduced_norm()/n_I
                        assert(n.is_integer())
                        n = ZZ(n)
                        if(is_prime(n)):
                            #should assert that el is in I
                            assert(ideal_element(el,I))
                            if(margin!=False):
                                t = polygen(ZZ)
                                K = NumberField(t**2 + 1,'a')
                                R = K.ring_of_integers()
                                if(n*n<margin*p and (not (R(n).is_square()))):
                                    return(I*(el.conjugate()/n_I))
                            else:
                                return(I*(el.conjugate()/n_I))
    return(False)


#Argument I is the target ideal, Ntor the torsion level
#Not sure this is actually sufficient, since el/nI is not in Z+Ntor*O
def equivalentPrimeIdeal_constraint(I, Ntor, bound=10,margin=False):
    p = I.quaternion_algebra().ramified_primes()[0]
    B = I.quaternion_algebra()
    Otor = B.quaternion_order(B.ideal(list((I.left_order()*Ntor).basis())+[B(1)]).basis())
    Itor = B.ideal(Otor.free_module().intersection(I.free_module()).basis())
    M = Itor.basis_matrix()
    M = M.LLL()
    Ba = Itor.basis()
    G = Itor.gram_matrix()
    U = G.LLL_gram().transpose()
    M = [sum(c * g for c, g in zip(row, Ba)) for row in U]
    n_I = I.norm()
    assert(n_I==Itor.norm())
    for x in range(bound):
        for a in range (0,x):
            for b in range (-x+a,x-a):
                for c in range (-x+a+b,x-a-b):
                    for d in range (-x+a+b+c,x-a-b-c):
                        v = vector([a,b,c,d])
                        el = sum( v[i]*M[i] for i in range(4))
                        n = el.reduced_norm()/n_I
                        assert(n.is_integer())
                        n = ZZ(n)
                        if(is_prime(n)):
                            #should assert that el is in I
                            assert(ideal_element(el,I))
                            if(margin!=False):
                                t = polygen(ZZ)
                                K = NumberField(t**2 + 1,'a')
                                R = K.ring_of_integers()
                                if(n*n<margin*p*Ntor and (not (R(n).is_square()))):
                                    return(el,I*(el.conjugate()/n_I))
                            else:
                                return(el,I*(el.conjugate()/n_I))
    return(False)












#Tests
def test_corvo(n, expected):
    r = Corvo_imperiale(n)
    if(r != False):
        x,y = r
        return(x**2+y**2 == n)
    return(expected==False)

def test_cornacchia(n, d, expected):
    r = Cornacchia_prime(n,d)
    if(r != False):
        x,y = r
        return(x**2+d*y**2 == n)
    return(expected==False)

def test_cornacchia_sequence():
    r = True
    r = r and test_corvo(5*5*4*13*101, True)
    r = r and test_corvo(5*5*4*13*103, False)
    r = r and test_corvo(5*5*4*13*13*13*401*101*997, True)
    r = r and test_corvo(13*101*991, False)
    r = r and test_cornacchia(7,1, False)
    r = r and test_cornacchia(2,1, True)
    r = r and test_cornacchia(7,3, True)
    return(r)

def test_isomorphism(I,J):
    assert(I.left_order()==J.left_order())
    a = ideal_isomorphism(I,J)
    return(I*a == J)

def test_element(a,I,contained=True):
    return(contained==ideal_element(a,I))

def test_ideal_add(I,M):
    n = len(M)
    m = matrix(QQ,M)
    J = ideal_add(I,m)
    i,j,k = I.quaternion_algebra().gens()
    B = I.basis()
    res = True
    for e in B:
        res = res and ideal_element(e,J)
    for x in range (n):
        e = m[x,0]+m[x,1]*i+m[x,2]*j+m[x,3]*k
        res = res and ideal_element(e,J)
    return(res)

def test_add_sequence():
    p = 11
    B = QuaternionAlgebra(p)
    i,j,k = B.gens()
    O = B.quaternion_order([1,1*i,(j+i)/2,(1+k)/2])
    x = 1+2*i+k
    I = O*x+O*8
    J = I*(i+k)
    x2 = (Integer(1)/Integer(5))*(3+21*i+15*k)
    res = test_ideal_add(I,[x2.coefficient_tuple()])
    res = res and test_ideal_add(J,[x2.coefficient_tuple()])
    res = res and test_ideal_add(J,[x.coefficient_tuple()])
    return(res)

def test_tools_sequence():
    p = 11
    B = QuaternionAlgebra(p)
    i,j,k = B.gens()
    O = B.quaternion_order([1,1*i,(j+i)/2,(1+k)/2])
    x = 1+2*i+k
    I = O*x+O*8
    J = I*(i+k)
    assert(J.norm().is_integer())
    r = test_isomorphism(I,J)
    r = r and test_isomorphism(J,I)
    r = r and test_element(1+2*i+k,I,True)
    r = r and test_element((i-5*j+2*k),I,True)
    r = r and test_element((i+k)/3,I,False)
    return(r)

def test_representInteger(p,n):
    assert(p %4 == 3)
    B = QuaternionAlgebra(p)
    r = representInteger(n,B)
    if(r!=False):
        return((r.reduced_norm()==n) and ((r[0]+r[3]).is_integer()) and ((r[1]+r[2]).is_integer()))
    
def test_representInteger_sequence():
    r = True
    r = r and test_representInteger(11,25)
    r = r and test_representInteger(103,1000)
    r = r and test_representInteger(103,1000)
    r = r and test_representInteger(19,100)
    return(r)

#Currently self-testing through asserts
def test_intersection(I):
    M4 = I.basis_matrix()
    L2 = Matrix(QQ,[[0,0,1,0],[0,0,0,1]])
    I2 = intersection(M4,L2)
    return(True)

def test_intersection_sequence():
    try: 
        p = 11
        B = QuaternionAlgebra(p)
        i,j,k = B.gens()
        O = B.quaternion_order([1,1*i,(j+i)/2,(1+k)/2])
        x = 1+2*i+k
        I = O*x+O*8
        r = test_intersection(I)
        x = 1 + i + 3*j + 3*k
        g = 7 + 24*i
        I = O*x+O*8
        I = (1/g)*I
        r = r and test_intersection(I)
        x = representInteger(200,B)
        g = representInteger(773,B)
        I = O*x+O*25
        I = (1/g)*I
        r = r and test_intersection(I)
        I = O*x+O*8
        g=representInteger(625,B)
        I = (1/g)*I
        r = r and test_intersection(I)
        return(r)
    except Exception:
        return False
    

#Currently self-testing through asserts
def test_intersection_ideals(I1,I2):
    I = intersection_ideals(I1,I2)
    b = I.basis()
    res = True
    for e in b:
        res = res and ideal_element(e,I1)
        res = res and ideal_element(e,I2)
    return(res)

def test_intersection_ideals_sequence():
    p = 11
    B = QuaternionAlgebra(p)
    i,j,k = B.gens()
    O = B.quaternion_order([1,1*i,(j+i)/2,(1+k)/2])
    N1 = 8
    x1 = representInteger(N1*13,B)
    I1 = O*x1+O*N1
    N2 = 7
    x2 = representInteger(N2*13,B)
    I2 = O*x2+O*N2
    r = test_intersection_ideals(I1,I2)
    return(r)

def test_connectingIdeal(I):
    r = True
    J = connectingIdeal(I.left_order(),I.right_order())
    r = r and I.right_order()==J.right_order()
    r = r and I.left_order()==J.left_order()
    return(r)

def test_connectingIdeal_sequence():
    p = 11
    B = QuaternionAlgebra(p)
    i,j,k = B.gens()
    O = B.quaternion_order([1,1*i,(j+i)/2,(1+k)/2])
    x = 1+2*i+k
    I = O*x+O*8
    r = test_connectingIdeal(I)
    x = 4+3*i+k
    I = O*x+O*3
    r = r and test_connectingIdeal(I)
    p = 103
    B = QuaternionAlgebra(p)
    i,j,k = B.gens()
    O = B.quaternion_order([1,1*i,(j+i)/2,(1+k)/2])
    x = 1+2*i+j
    I = O*x+O*4
    r = r and test_connectingIdeal(I)
    return(r)

def test_idealPullPush(Ic,I):
    Ipf = I
    r = True
    if(I.left_order()==Ic.left_order()):
        Ipf = idealPushforward(Ic,I)
        r = r and Ipf.norm()==I.norm()
        r = r and Ipf.left_order()==Ic.right_order()
    else:
        assert(I.left_order()==Ic.right_order())
    Ipl = idealPullback(Ic,Ipf)
    r = r and Ipl.norm()==I.norm()
    r = r and Ipl.left_order()==Ic.left_order()
    Ipf2 = idealPushforward(Ic,I)
    r = r and Ipf2.norm()==I.norm()
    r = r and Ipf2.left_order()==Ic.right_order()
    r = r and Ipf2 == Ipf
    Ipl2 = idealPullback(Ic,Ipf)
    r = r and Ipl2.norm()==I.norm()
    r = r and Ipl2.left_order()==Ic.left_order()
    r = r and Ipl2 == Ipl
    return(r)

def test_idealPullPush_sequence():
    p = 11
    B = QuaternionAlgebra(p)
    i,j,k = B.gens()
    O = B.quaternion_order([1,1*i,(j+i)/2,(1+k)/2])
    x = 1+2*i+k
    I = O*x+O*8
    x2 = 4+3*i+k
    I2 = O*x+O*3
    r = test_idealPullPush(I,I2)
    r = r and test_idealPullPush(I2,I)
    p = 103
    B = QuaternionAlgebra(p)
    i,j,k = B.gens()
    O = B.quaternion_order([1,1*i,(j+i)/2,(1+k)/2])
    x = 2+2*i+2*j
    I = O*x+O*10
    I1 = I+O*2
    I2 = I1.conjugate()*I*(1/I1.norm())
    assert(I1*I2==I)
    r = r and test_idealPullPush(I1,I2)
    return(r)
    

def test_equivalentPrimeIdeal(I):
    r = True
    J = equivalentPrimeIdeal(I)
    r = r and J.norm().is_integer()
    r = r and ZZ(J.norm()).is_prime()
    r = r and (I.left_order() == J.left_order())
    return(r)


def test_equivalentPrimeIdeal_sequence():
    p = 11
    B = QuaternionAlgebra(p)
    i,j,k = B.gens()
    O = B.quaternion_order([1,1*i,(j+i)/2,(1+k)/2])
    x = 1+2*i+k
    I = O*x+O*8
    r = test_equivalentPrimeIdeal(I)
    x = 4+3*i+k
    I = O*x+O*3
    r = r and test_equivalentPrimeIdeal(I)
    p = 103
    B = QuaternionAlgebra(p)
    i,j,k = B.gens()
    O = B.quaternion_order([1,1*i,(j+i)/2,(1+k)/2])
    x = 1+2*i+j
    I = O*x+O*4
    r = r and test_equivalentPrimeIdeal(I)
    return(r)



def test_equivalentPrimeIdeal_constraint(I, Ntor):
    r = True
    B = I.quaternion_algebra()
    el,J = equivalentPrimeIdeal_constraint(I,Ntor)
    r = r and J.norm().is_integer()
    r = r and ZZ(J.norm()).is_prime()
    r = r and (I.left_order() == J.left_order())
    r = r and ideal_element(el,I)
    Otor = B.quaternion_order(B.ideal(list((I*Ntor).basis())+[B(1)]).basis())
    #r = r and ideal_element((1/I.norm())*el.conjugate(),Otor*1)
    r = r and ideal_element(el,Otor*1)
    return(r)


def test_equivalentPrimeIdeal_constraint_sequence():
    p = 11
    B = QuaternionAlgebra(p)
    i,j,k = B.gens()
    O = B.quaternion_order([1,1*i,(j+i)/2,(1+k)/2])
    x = 1+2*i+k
    I = O*x+O*8
    r = test_equivalentPrimeIdeal_constraint(I, 5)
    x = 4+3*i+k
    I = O*x+O*3
    r = r and test_equivalentPrimeIdeal_constraint(I, 5)
    p = 103
    B = QuaternionAlgebra(p)
    i,j,k = B.gens()
    O = B.quaternion_order([1,1*i,(j+i)/2,(1+k)/2])
    x = 1+2*i+j
    I = O*x+O*4
    r = r and test_equivalentPrimeIdeal_constraint(I,17)
    return(r)


if __name__=="__main__":
    res = True
    res = res and test_cornacchia_sequence()
    res = res and test_add_sequence()
    res = res and test_tools_sequence()
    res = res and test_representInteger_sequence()
    res = res and test_intersection_sequence()
    res = res and test_intersection_ideals_sequence()
    res = res and test_connectingIdeal_sequence()
    res = res and test_idealPullPush_sequence()
    res = res and test_equivalentPrimeIdeal_sequence()
    res = res and test_equivalentPrimeIdeal_constraint_sequence()
    print(res)

