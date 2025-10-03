from KLPT import *
from BorelKLPT import SigningBorelKLPT
from ScalarKLPT import SigningScalarKLPT
from setup import *


def KLPT_setup(all_prime=False):
    """
    Helper for tests
    """
    Ntor =3**(log(p,3).round()//2)# 3**(log(p,3).round()//4)*5**(log(p,5).round()//4) #product can be solved, but the mod sth are often badly behaved then
    Ntor_factor = list(Ntor.factor())
    if all_prime:
        Ntor = sqrt(p).round().next_prime()
        Ntor_factor=[(Ntor,1)]
    for _ in range(100):
        c = randint(p+2, p^2)
        γ = RepresentIntegerHeuristic(Ntor*c)
        if γ: break
    Itorp = O0*γ+O0*Ntor
    Mc=0
    M=0
    for _ in range(100):
        Mc = Ntor
        while not gcd(Ntor,Mc)==1:
            Mc = (Integer(randint(p+2,p^2)))
        if all_prime:
            Mc = (Integer(randint(2^(log(p,2)/2-2).round(),2^(log(p,2)/2).round()))).next_prime()
        c = 0
        while(c*Mc)<Mc*log(Mc).round():
            c = randint(p+2, p^2)
            c = c//gcd(Mc,c)
        γ = RepresentIntegerHeuristic(Mc*c)
        if γ:# and (Mc%4==1): 
            break
    Ic = O0*γ+O0*Mc
    Ic,_ = make_cyclic(Ic)
    Mc = ZZ(Ic.norm())
    for t in range(100):
        M = randint(p+2, p^2)
        if all_prime:
            M = Integer(randint(5, p.sqrt().round())).next_prime()
        c=0
        while(c*M)<M*log(M).round():
            c = randint(p+2, p^2)
            c = c//gcd(M,c)
        γ = RepresentIntegerHeuristic(M*c)
        if γ and (gcd(Mc,M)==1): 
            print(str(t)+" attempts for Ic")
            break
    Ip = O0*γ+O0*M
    I = pushforward_ideal(O0,Ic.right_order(),Ip,Ic)
    I,_ = make_cyclic(I)
    Itor = pushforward_ideal(O0,Ic.right_order(),Itorp,Ic)
    #assert(gcd(ZZ(Itor.norm()),ZZ(I.norm())==1))
    #assert(gcd(ZZ(Itor.norm()),ZZ(Ic.norm())==1))
    #assert(gcd(ZZ(Ic.norm()),ZZ(I.norm())==1))
    return(I,Ic,Itor, Ntor_factor)

def test_SigningKLPT(self):
    I,Ic,Itor = self.KLPT_setup()
    J = SigningKLPT(I,Ic)
    print(log(J.norm(),2),J.is_left_equivalent(I))

def test_SigningBorelKLPT():
    I,Ic,Itor, Ntor_factor = KLPT_setup(all_prime=False)
    print(I,Ic,Itor)
    print(I.norm(),Ic.norm(),Itor.norm(),Ntor_factor)
    J = SigningBorelKLPT(I,Ic,Itor,Ntor_factor)
    print(log(J.norm(),2))
    print("order test", J.left_order()==I.left_order())

def test_SigningScalarlKLPT():
    I,Ic,Itor, Ntor_factor = KLPT_setup(all_prime=False)
    Ntor = ZZ(Itor.norm())
    #Ntor = (Ntor//27).next_prime()
    #Ntor_factor = [(Ntor, 1)]
    print(I,Ic,Ntor)
    print(I.norm(),Ic.norm(),Ntor,Ntor_factor)
    J = SigningScalarKLPT(I,Ic,Ntor,Ntor_factor)
    print(log(J.norm(),2))
    print("order test", J.left_order()==I.left_order())


if __name__ == '__main__' and '__file__' in globals():
    test_SigningBorelKLPT()
    test_SigningScalarlKLPT()

