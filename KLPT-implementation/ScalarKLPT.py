from sage.all import Zmod, ZZ, inverse_mod, gcd, QQ, vector, CRT, factor, floor, log, Zmod, ceil


# Local imports
from ideals import (
    chi,
    is_cyclic,
    pullback_ideal,
    pushforward_ideal,
    equivalent_left_ideals,
    quaternion_basis_gcd,
    eichler_order_from_ideal,
    make_cyclic,
    quadratic_norm,
)
from setup import (
    O0,
    p,
    q,
    j,
    ω,
    l,
    eScalar,
    logp,
    loglogp
)
import random

from rerandomization import KLPT_input_rerandomization

from KLPT import (eichler_order_from_ideal)

from BorelKLPT import (prod, 
                       Borel_small_equivalent_ideal,
                       EquivalentPrimeIdealConstraint, 
                       EquivalentPrimeIdealfromLattice, 
                       full_strong_approximation_lattice_heuristic)


from KLPT import (eichler_order_from_ideal,
                  IdealModConstraint,
                  EichlerModConstraint)

#Given a maximal order O and a non-zero integer N, compute Z+NO as quaternion order
def compute_scalar_order(Ntor,maximal_order):
    M = (maximal_order*Ntor).basis_matrix()
    d = M.denominator()
    M5 = M.transpose().augment(vector([1,0,0,0])).transpose()
    Mi = M5*d
    Mi = Mi.change_ring(ZZ)
    Mi = Mi.echelon_form()[:4]
    Mq = Mi.change_ring(QQ)*(1/d)
    b = [Mq[x,0]+Mq[x,1]*ω + Mq[x,2]*j + Mq[x,3]*ω*j for x in range(4)]
    O = maximal_order.quaternion_algebra().quaternion_order(b)
    assert((maximal_order*Ntor).basis()[0] in O)
    assert((maximal_order*Ntor).basis()[1] in O)
    assert((maximal_order*Ntor).basis()[2] in O)
    assert((maximal_order*Ntor).basis()[3] in O)
    return(O)


#EquivalentPrimeIdeal in Z+NO, with norm square mod N if needed
def EquivalentPrimeIdeal_ZNO(I,N, square_constraint = True):
    # the max(N*N,p**(2/3)) is there for making tests with extremely small N run
    # normally it should equal N*N
    norm_bound = round((p**(6/7)))*I.norm()*max(N*N,round(p**(2/3)))
    norm_bound = norm_bound+10*(log(norm_bound,2).round())
    squareMod = False
    if(square_constraint):
        squareMod=N
    L,NL,δ = EquivalentPrimeIdealfromLattice(I,I.intersection(compute_scalar_order(N,I.left_order())),norm_bound, squareMod=squareMod)
    if(square_constraint):
        while(not Zmod(N)(NL).is_square()):
            L,NL,δ = EquivalentPrimeIdealfromLattice(I,I.intersection(compute_scalar_order(N,I.left_order())),norm_bound)
    return (L,ZZ(NL),δ)

#Pullback to O0 and equivalentPrimeIdeal in the order Z+NO0.
def derive_L_ScalarKLPT(I, Iτ, Nτ, Ntor,O0, O1, allow_nonsquare=True):
    for _ in range(20):

        # I' = [Iτ]^* I
        K_prime = pullback_ideal(O0, O1, I, Iτ)
        assert(K_prime.left_order()==O0)

        L,N,δ = EquivalentPrimeIdeal_ZNO(K_prime,Ntor,not allow_nonsquare)
        assert(L.left_order()==O0)

        # Bad delta, this will cause EichlerModConstraint to break
        if gcd(δ.reduced_norm(), Nτ) != 1 and gcd(δ.reduced_norm(), Ntor) != 1:
            print(f"DEBUG [SigningScalarKLPT]: Not sure why this is happening...")
            print(f"{factor(δ.reduced_norm()) = }")
            print(f"{factor(Nτ.reduced_norm()) = }")
            print(f"{gcd(δ.reduced_norm(), Nτ) = }")
            continue

        if L is not None:
            return L, N, δ

    # If we get here, something is likely broken
    raise ValueError(f"Never found an equivalent prime norm ideal")


#Compute L2 as 2^(e_Scalar-e1+extra)
#This gives the output size for StrongApproximation
def derive_L2_ScalarSigningKLPT(γ, L1, e1):
    """
    Given L1 = l^e1 and γ try and compute L2
    so that the output of SigningKLPT has norm
    exactly 2^e
    """
    g = quaternion_basis_gcd(γ, O0)
    extra = 2 * (floor(loglogp / 4) + ZZ(gcd(g, L1).valuation(l)))
    e2 = eScalar - e1 + extra
    return l**e2


def get_lambda(factors,mod,C,D, L2):
    K = Zmod(mod)
    try:
        λλ = K(L2) / K(p * quadratic_norm(C, D))
    except Exception as e:
        print(f"ERROR [get_lambda]: lambda failure e {e}")
        return None
    λs = []
    for x in factors:
        if not Zmod(x)(λλ).is_square():
            return None
        λs.append(int(Zmod(x)(λλ).sqrt()))
    λ = CRT(λs, factors)
    if λ == 0:
        return None
    return(λ)


def derive_C_and_D_SigningScalarKLPT(L, N, Iτ, Nτ, EichlerIτ, Ntor, Ntor_factor, Otor, γ, δ, L2, C1, D1):
    Cc, Dc =  EichlerModConstraint(Iτ,EichlerIτ,γ,δ)
    CI, DI =  IdealModConstraint(L,γ)
    D2 = random.randrange(0,Ntor)
    C2 = D2*C1*inverse_mod(D1,Ntor)
    C = CRT([CI,Cc,C2],[N,Nτ,Ntor])
    D = CRT([DI,Dc,D2],[N,Nτ,Ntor])
    μ_norm = p * quadratic_norm(C, D)

    KL = Zmod(N)
    Kτ = Zmod(Nτ)
    Ktor = Zmod(Ntor)

    try:
        square_mod_N = (KL(L2) / KL(μ_norm)).is_square()
        square_mod_Nτ = (Kτ(L2) / Kτ(μ_norm)).is_square()
        square_mod_Ntor = (Ktor(L2) / Ktor(μ_norm)).is_square()
    except:
        print("DEBUG [derive_C_and_D_SigningScalarKLPT]: is_square exception")
        return None, None
    # To compute the square root mod N*Nτ*Ntor
    # All of these must be true. Will happen
    # about 1/8 of the time.
    assert((γ*j*(C+D*ω)) in compute_scalar_order(Ntor,O0))
    assert((γ*j*(C+D*ω)) in L)
    assert((γ*j*(C+D*ω)*δ) in compute_scalar_order(Ntor,O0))
    assert((γ*j*(C+D*ω)*δ) in EichlerIτ)
    if square_mod_N and square_mod_Nτ and square_mod_Ntor:
        if  Ktor(D).is_unit() and Kτ(D).is_unit() and KL(D).is_unit():
            return C, D
    return None, None



#Assume Ntor (power of) prime for the moment, and N(I) must be a square mod Ntor
def SigningScalarKLPT_with_square_NI_mod_Ntor(I0, Iτ0, Ntor, Ntor_factor,allow_2_not_square_mod_Ntor=True):
    assert(Zmod(Ntor)(ZZ(I0.norm())).is_square())

    assert is_cyclic(I0), "I is not cyclic"
    assert is_cyclic(Iτ0), "Iτ is not cyclic"

    assert prod([a**b for (a,b) in Ntor_factor])==Ntor
    # Make I as small as possible
    Iτ1,alpha0 = Borel_small_equivalent_ideal(Iτ0)
    Iτ,_,alpha = EquivalentPrimeIdealConstraint(Iτ1,Ntor)
    ac0 = alpha0.conjugate()/(Iτ0.norm())
    ac1 = alpha.conjugate()/(Iτ1.norm())
    ac = ac0*ac1
    I1 = (1/ac)*I0*ac
    I, alphaI = I1,1+j-j
    # Prime norm ideal
    Nτ = ZZ(Iτ.norm())
    Nτtor = ZZ(Ntor)*ZZ((Iτ).norm())
    #I,_,alphaI = EquivalentPrimeIdeal_ZNO(I1,Ntor)
    I,_,alphaI = KLPT_input_rerandomization(I1,compute_scalar_order(Ntor,I1.left_order()),Nτ)
    aI = alphaI.conjugate()/(I1.norm())


    N=ZZ(I.norm())
    assert Nτ.is_prime()
   # assert Ntor.is_prime()
    assert gcd(N,Nτ)==1

    # Orders needed for pushback and pullforward
    O1 = I.left_order()
    assert Iτ.left_order() == O0
    assert Iτ.right_order() == O1

    # EichlerIτ = ℤ + Iτ = OL(I) ∩ OR(I)
    EichlerIτ = eichler_order_from_ideal(Iτ)
    EichlerIτtor = eichler_order_from_ideal(Iτ).intersection(compute_scalar_order(Ntor,O0))
    Otor = compute_scalar_order(Ntor,Iτ.right_order())
    #assert N.is_prime()
    assert Nτ.is_prime()
    L, N, δ = derive_L_ScalarKLPT(I, Iτ, Nτ, Ntor,O0, O1, False)
    assert L.left_order()==O0
    assert N.is_prime()

    print("INFO [SigningScalarKLPT]: Start loop")
    for _ in range(8000):
        # Compute the pullback K of I with left order O0, and
        # find an equivalent prime norm ideal L ~ K.
        # Only really random step in the algorithm, since RepresentInteer is deterministic
        # Equivalence after pullback to O0 is taken so that torsion influence is not impacted.
        # Not yet clear how exactly to set e1. For now just added a factor 3 and doubled the part of logp in the LearningToSQI formula
        #e1 = 4*(floor(2*logp - log(N, l) + 1.74 * loglogp)+4)
        C1,D1 = random.randrange(0,Ntor+1),1
        #this would be better, but for now strongApproximation does not handle D1 non coprime to Ntor
        #if(C1==Ntor):
        #    C1,D1 = 1,0
        
        e1 = floor(2*logp+2*loglogp+2*log(Ntor,l)-log(N,l))
        L1 = l**e1
        l1 = get_lambda([a**b for a,b in Ntor_factor],Ntor,C1,D1, N*L1)
        if l1 is None:
            print("INFO [SigningScalarKLPT]: No l1 found")
            continue
        assert((l1*l1*p*(C1*C1+D1*D1))%(Ntor)==L1*N%(Ntor))
        γ = full_strong_approximation_lattice_heuristic(Ntor,C1,D1,l1,N*L1, small_power_of_two=False)
        if γ is None:
            print("DEBUG [SigningScalarKLPT]: No γ found")
            continue
        
        
        # Compute the output size for StrongApproximation so that SigningKLPT output has size 2^e_Scalar
        # At least approximatel,as we do not enforce strict equality yet
        L2 = derive_L2_ScalarSigningKLPT(γ, L1, e1)
        # Find W,X such that L2 / norm(W + ωXNtor) is a square mod N*Nτ*p
        # And also  γ(W + ωXNtor)δ in EichlerIτ and Z+NtorO0, and γ(W + ωXNtor) in I
        # And X invertible mod Ntor, N, Nτ and p
        # Returns None,None if no such pair is found
        C,D = derive_C_and_D_SigningScalarKLPT(L, N, Iτ, Nτ, EichlerIτ, Ntor, Ntor_factor, Otor, γ, δ, L2, C1, D1)
        if C is None:
            print(f"DEBUG [SigningScalarKLPT]: Unable to compute a C,D given γ.")
            continue
        # Compute the square root λ of L2 / norm(Cj + jωD) mod N*Nτ*p
        λ = get_lambda([a**b for a,b in Ntor_factor]+[N,Nτ],N*Nτ*Ntor,C,D, L2)
        if( λ is None):
            print(f"DEBUG [SigningScalarKLPT]: No λ given C and D.")
            continue
        assert((λ*λ*p*(C*C+D*D))%(N*Ntor*Nτ)==L2%(N*Ntor*Nτ))

        # StrongApproxinmation. L2 is a power of 2, C,D and λ as computed above
        assert((γ*j*(C+D*ω)*δ) in EichlerIτ)
        assert((γ*λ*j*(C+D*ω)*δ) in EichlerIτ)
        μ = full_strong_approximation_lattice_heuristic(Ntor*N*Nτ,C,D,λ,L2)
        # No solution found, try another γ
        if not μ:
            print(f"DEBUG [SigningScalarKLPT]: StrongApproximation failed!", N, Nτ, Ntor, p, L2)
            continue
        assert((γ*μ *δ) in EichlerIτ)
        # Strong approximation norm check
        assert μ.reduced_norm() == L2
        print(f"INFO [SigningScalarKLPT]: Found a solution to the StrongApproximation!")
        # Now construct the equivalent ideal J
        β = γ * μ
        assert β in compute_scalar_order(Ntor,O0)
        assert β in L
        J_prime = chi(β, L)
        # J = [Iτ]_* J'
        J = pushforward_ideal(O0, O1, J_prime, Iτ)
        # Check things are actually equivalent
        assert equivalent_left_ideals(I, J)
        assert equivalent_left_ideals(I0, ac*J*(1/ac))
        # Check torsion is right
        assert β*δ in compute_scalar_order(Ntor,O0)
        assert β*δ in EichlerIτ
        assert β*δ in EichlerIτtor
        assert ac*β*δ*(1/ac) in compute_scalar_order(Ntor,I0.left_order())
        # Make sure the output is cyclic
        J, _ = make_cyclic(J)
        #Now the norm might not be exactly 2^e but a few bits less (or more). We don't care for now
        print(f"INFO [SigningScalarKLPT]: Preliminary solution of norm {factor(J.norm()) = }")
        print(f"INFO [SigningScalarKLPT]: Returning solution of norm {factor(J.norm()) = }")
        return ac*J*(1/ac)

# Assume Ntor (power of) prime for the moment, and requires 2 a non-square mod Ntor
def SigningScalarKLPT(I0, Iτ0, Ntor, Ntor_factor,allow_2_not_square_mod_Ntor=True):
    #Find gamma invertible mod Ntor
    # find n invertible mod Ntor, get endo gamma of norm n, find n
    assert(gcd(ZZ(Ntor),ZZ(I0.norm()))==1)
    J0 = None
    O = I0.left_order()
    if(not Zmod(Ntor)(ZZ(I0.norm())).is_square()):
        assert(not Zmod(ZZ(Ntor))(2).is_square())
        gamma = I0.left_order().random_element()
        Igamma = O0*1
        while(Igamma.norm()!=2):
            while(ZZ(gamma.reduced_norm()%2!=0)):
                gamma = O.random_element()
            Igamma = gamma*O+O*2
        assert Igamma.norm()==2
        #assert I0.left_order()==I.left_order()
        Imul = Igamma*I0
        assert Zmod(ZZ(Ntor))(ZZ(Imul.norm())).is_square()
        Ic = Iτ0*Igamma.conjugate()
        Ic,_ = make_cyclic(Ic)
        assert Ic.right_order()==Igamma.left_order()
        J = SigningScalarKLPT_with_square_NI_mod_Ntor(Imul,Ic,Ntor,Ntor_factor)
        assert(Igamma.conjugate().right_order()==J.left_order())
        J0 = Igamma.conjugate()*J
        assert(J0.norm()==J.norm()*2)
        assert(J0.left_order()==I0.left_order())
        J0, _ = make_cyclic(J0)
        print(f"INFO [SigningScalarKLPT]: Success for non-square mod Ntor")
    else :
        J0 = SigningScalarKLPT_with_square_NI_mod_Ntor(I0,Iτ0,Ntor,Ntor_factor)
        print(f"INFO [SigningScalarKLPT]: Success for square mod Ntor")
    return(J0)

