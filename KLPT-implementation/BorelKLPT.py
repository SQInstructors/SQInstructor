# Python imports
from random import choice

# Sage imports
from sage.all import (
    gcd,
    ZZ,
    Zmod,
    log,
    ceil,
    floor,
    sqrt,
    flatten,
    factor,
    Matrix,
    vector,
    randint,
    inverse_mod,
    prod,
    choice,
    CRT,
    Integer,
)

# Local imports
from ideals import (
    chi,
    is_cyclic,
    reduced_basis,
    pullback_ideal,
    pushforward_ideal,
    equivalent_left_ideals,
    quaternion_basis_gcd,
    eichler_order_from_ideal,
    make_cyclic,
    quadratic_norm,
)
from utilities import Cornacchia
from lattices import generate_small_norm_quat, generate_close_vectors
from setup import (
    O0,
    p,
    q,
    j,
    ω,
    l,
    eBorel,
    logp,
    loglogp,
    prime_norm_heuristic,
    represent_heuristic,
)

from KLPT import (eichler_order_from_ideal,
                  EquivalentPrimeIdealHeuristic,
                  IdealModConstraint,
                  EichlerModConstraint,
                  strong_approximation_construct_lattice)

def prod(liste):
    n = 1
    for elem in liste:
        n = n*elem
    return(n)

def derive_L2_BorelSigningKLPT(γ, L1, e1):
    """
    Given L1 = l^e1 and γ try and compute L2
    so that the output of SigningKLPT has norm
    exactly 2^e
    """
    g = quaternion_basis_gcd(γ, O0)
    extra = 2 * (floor(loglogp / 4) + ZZ(gcd(g, L1).valuation(l)))
    e2 = eBorel - e1 + extra
    return l**e2

def Borel_small_equivalent_ideal(I, reduced_basis_elements=None):
    """
    Computes the Minkowski reduced basis of the
    ideal and returns the smallest equivalent
    ideal J = Iα ~ I.
    """
    nI = I.norm()

    if not reduced_basis_elements:
        reduced_basis_elements = reduced_basis(I)

    b0 = reduced_basis_elements[0]
    if b0.reduced_norm() == nI**2:
        return I,b0/b0
    return I * I.right_order().left_ideal([b0.conjugate() / nI]),b0



def prime_norm_algebra_element_constraint(
    nI,
    Ibasis,
    Ntor,
    coeff_bound,
    search_bound,
    previous=set(),
    allowed_factors=None,
    random_elements=False,
):
    """
    Find an element α ∈ B_{0, ∞} with small,
    prime scaled norm.

    Optional: `allowed_factors` allows the norm to
               be composite, where it is expected that
               the result is a large prime multiplied by
               small factors dividing allowed_factors.
    """

    if random_elements:
        small_elements = generate_small_norm_quat_random(
            Ibasis, coeff_bound, search_bound*2
        )
    else:
        max_norm_bound = prime_norm_heuristic
        small_elements = generate_small_norm_quat(
            Ibasis, max_norm_bound, count=search_bound*2
        )

    for α in small_elements:
        α_norm = ZZ(α.reduced_norm()) // nI

        # Even norms can be rejected early
        # as allowed_factors is either odd
        # or None.
        if α_norm % 2 == 0:
            continue

        # We can allow α to have composite norm
        # if small factors are within the KLPT
        # target norm T.
        α_norm_reduced = α_norm
        if allowed_factors:
            g = gcd(α_norm, allowed_factors)
            if g != 1:
                α_norm_reduced //= g

        # If we've failed with this norm before
        # continue
        if α_norm_reduced in previous:
            continue

        # Check if the element has prime norm
        if α_norm_reduced.is_pseudoprime():
            # Check if the prime small enough
            if α_norm_reduced < 2*prime_norm_heuristic:
                if(Zmod(α_norm_reduced)(p).is_square()==Zmod(Ntor)(p).is_square()):
                    return α, α_norm
    return None, None


def EquivalentPrimeIdealConstraint(
    I, Ntor, previous=set(), allowed_factors=None, random_elements=False
):
    """
    Given an ideal I with norm nI, attempts
    to find an equivalent ideal J with prime norm.

    If unsuccessful, returns None
    """
    # TODO: what's a good initial small search?
    coeff_bound = max((floor(logp / 10)), 7)

    # TODO: what's a good search bound?
    search_bound = max(coeff_bound**4, 4096)

    # Norm of Ideal
    nI = ZZ(I.norm())

    # Compute the Minkowski reduced basis
    Ibasis = reduced_basis(I)

    # Find an element with small prime norm
    α, N = prime_norm_algebra_element_constraint(
        nI,
        Ibasis,
        Ntor,
        coeff_bound,
        search_bound,
        previous=previous,
        allowed_factors=allowed_factors,
        random_elements=random_elements,
    )

    if α is None:
        print(f"DEBUG [EquivalentPrimeIdealHeuristic] No equivalent prime found")
        return None, None, None

    assert ZZ(α.reduced_norm()) // nI == N
    assert α in I

    # Compute the ideal given α
    J = chi(α, I)

    return J, N, α


def EquivalentPrimeIdealfromLattice(I, lattice,norm_bound=False, squareMod=False):
    """
    Given an ideal I with norm nI, attempts
    to find an equivalent ideal J with prime norm.
    Equivalence is taken in lattice (intersection of I with some order)

    If unsuccessful, returns None
    """

    # Norm of Ideal
    nI = ZZ(I.norm())

    # Compute the Minkowski reduced basis
    latbasis = reduced_basis(lattice)
    if(norm_bound==False):
        norm_bound = 2**(logp+floor(logp*2/5))
        norm_bound = norm_bound*nI
    
    abound = floor(((norm_bound/latbasis[3].reduced_norm())/4).sqrt())
    assert(abound>0)
    assert(abound>log(latbasis[3].reduced_norm(),2))

    a,b,c,d = randint(0,abound),randint(0,abound),randint(0,abound),randint(0,abound)
    α = a*latbasis[0]+b*latbasis[1]+c*latbasis[2]+d*latbasis[3]
    counter = min(2*abound, 10000)
    while( (counter > 1) and (α.reduced_norm()>norm_bound) or (not ZZ(α.reduced_norm()/I.norm()).is_prime()) or (not((not squareMod) or (Zmod(squareMod)(ZZ(α.reduced_norm()/I.norm())).is_square())))):
        a,b,c,d = randint(0,abound),randint(0,abound),randint(0,abound),randint(0,abound)
        α = a*latbasis[0]+b*latbasis[1]+c*latbasis[2]+d*latbasis[3]

    if α is None:
        print(f"DEBUG [EquivalentPrimeIdealfromLattice] No equivalent prime found")
        return None, None, None

    N = ZZ(α.reduced_norm() // nI)

    if (squareMod):
        assert Zmod(squareMod)(N).is_square()
    assert ZZ(α.reduced_norm()) // nI == N
    assert ZZ(α.reduced_norm()) // nI == N
    assert α in I

    # Compute the ideal given α
    J = chi(α, I)
    assert J.norm()==N
    assert ZZ(N).is_prime()

    return J, N, α


def BorelDerive_L(I, Iτ, Nτ,Itor, Ntor,O0, O1):
    """
    Given an ideal I with left order O1 and an ideal
    Iτ with left order O0 of norm Nτ computes an ideal
    equivalent to the pullback of I under Iτ with prime
    norm.

    Input: I with left order O1
           Iτ with left order O0 and norm Nτ

    Output L ~ [Iτ]^* I with prime norm
           N = n(L)
           δ such that L = χ(K', δ)
    """
    for _ in range(20):

        # I' = [Iτ]^* I
        K_prime = pullback_ideal(O0, O1, I, Iτ)
        assert K_prime.left_order()==O0

        Ksmall, δsmall = Borel_small_equivalent_ideal(K_prime)
        assert(Ksmall==chi(δsmall,K_prime))
        L, N, δ = EquivalentPrimeIdealHeuristic(Ksmall)
        assert(L==chi(δ,Ksmall))
        δ = δ*(1/δsmall.conjugate())*K_prime.norm()
        assert(chi(δ,K_prime)==L)
        assert L.left_order()==O0

        # Bad delta, this will cause EichlerModConstraint to break
        if gcd(δ.reduced_norm(), Nτ) != 1 and gcd(δ.reduced_norm(), Ntor) != 1:
            print(f"DEBUG [SigningKLPT]: Not sure why this is happening...")
            print(f"{factor(δ.reduced_norm()) = }")
            print(f"{factor(Nτ.reduced_norm()) = }")
            print(f"{gcd(δ.reduced_norm(), Nτ) = }")
            continue

        if L is not None:
            return L, N, δ

    # If we get here, something is likely broken
    raise ValueError(f"Never found an equivalent prime norm ideal")




def FullRepresentIntegerHeuristic(M, parity=False):
    """
    Algorithm 1 (Page 8)

    Given an integer M, with M > p, attempts to
    find a random element γ with norm M.

    If no element is found after `bound` tries,
    returns none
    """

    def RepresentInteger(M, z, t, parity=False):
        M_prime = M - p * quadratic_norm(z, t)
        two_squares = Cornacchia(M_prime, -ZZ(ω**2))
        if two_squares:
            x, y = two_squares
            if parity and (x + t) % 2 == 0 and (y + z) % 2 == 0:
                return None
            return x + ω * y + j * (z + ω * t)
        # No solution for the given M
        return None

    if M <= p:
        raise ValueError(f"Can only represent integers M > p.")
    m = max(floor(sqrt(M / (p * (1 + q)))), 5)

    M = 4*M
    for _ in range(m**2):
        z = randint(-m, m)
        t = randint(-m, m)
        γ = RepresentInteger(M, z, t, parity=parity)
        if(γ is not None):
            a,b,c,d = γ.coefficient_tuple()
            if not ((ZZ(a+d)%2==0)and(ZZ(b+c)%2==0)):
                γ =None
            else:
                γ = γ/2

        if γ is not None:
            # Found a valid solution, return
            assert γ.reduced_norm() == M/4, "The norm is incorrect"
            assert γ in O0, "The element is not contained in O0"
            return γ

    # No solution found, return None
    print(f"DEBUG [RepresentIntegerHeuristic]: No solution found")
    return None



def derive_C_and_D_SigningBorelKLPT(L, N, Iτ, Nτ, EichlerIτ, Itor, Ntor, Ntor_factor, EichlerItor, γ, δ, L2):
    """
    Solves IdealModConstraint and EichlerModConstraint
    for a given γ and returns when we find an element
    μ such that L2 / Nrd(μ) is a square mod N*Nτ*Ntor

    Input: Ideals L, Iτ, Itor of norm N, Nτ, Ntor prime
           EichlerIτ the Eichler order of Iτ
           EichlerItor the Eichler order of Itor
           Iτ and L left O0-ideals
           Iτ*Itor well-defined
           γ, δ, elements of B
           L2, a divisor of T

    Output C, D such that L2 / p*(C^2 + D^2) is a square mod N*Nτ*Ntor
    and Cj+Dij in L and in EichlerIτ intersection EichlerItor
    """

    Etors = []
    Iconn = Itor.left_order()*1
    for a,b in Ntor_factor:
        shortened = Itor+(a**b)*Itor.left_order()
        ps = pushforward_ideal(Iconn.left_order(), Iconn.right_order(),shortened,Iconn)
        Etors.append(eichler_order_from_ideal(ps))
        Iconn = Iconn*ps


    C0, D0 = IdealModConstraint(L, γ)
    C1, D1 = EichlerModConstraint(Iτ, EichlerIτ, γ, δ)
    CDs = [EichlerModConstraint(Iτ, Et, γ, δ) for Et in Etors]
    Cs, Ds, Ns = [a for a,_ in CDs], [b for _,b in CDs], [a**b for a,b in Ntor_factor]

    # Compute CRT
    C = CRT([C0, C1] + Cs, [N, Nτ] + Ns)
    D = CRT([D0, D1] + Ds, [N, Nτ] + Ns)

    # We need to take a sqrt of this
    μ_norm = p * quadratic_norm(C, D)

    KL = Zmod(N)
    Kτ = Zmod(Nτ)
    Ktor = Zmod(Ntor)

    try:
        square_mod_N = (KL(L2) / KL(μ_norm)).is_square()
        square_mod_Nτ = (Kτ(L2) / Kτ(μ_norm)).is_square()
        square_mod_Ntor = (Ktor(L2) / Ktor(μ_norm)).is_square()
    except:
        print("DEBUG [derive_C_and_D_SigningBorelKLPT]: is_square exception")
        return None, None
    if square_mod_N and square_mod_Nτ and square_mod_Ntor:
        if  Ktor(D).is_unit() and Kτ(D).is_unit() and KL(D).is_unit():
            return C, D
    return None, None




def full_strong_approximation_lattice_heuristic(N, C, D, λ, L, small_power_of_two=False):
    """
    Constructs a lattice basis and then looks for
    close vectors to the target.

    Allows for optimising output from pN^4 to pN^3,
    which helps keep the norm small and hence the
    degree of the isogenies small
    """
    L = 4*L
    λ = 2*λ
    # We really only expect this for the case when N = 2^k
    swap = False
    if D == 0 or gcd(D, N) != 1:
        C, D = D, C
        swap = True

    # Construct the lattice
    lattice_basis, target, zp0, tp0 = strong_approximation_construct_lattice(
        N, C, D, λ, L, small_power_of_two=small_power_of_two
    )

    # Generate vectors close to target
    close_vectors = generate_close_vectors(lattice_basis, -target, p, L)

    xp, yp = None, None
    for close_v in close_vectors:
        zp, tp = close_v
        assert zp % N == 0, "Can't divide zp by N"
        assert tp % N == 0, "Can't divide tp by N"

        zp = ZZ(zp / N) + zp0
        tp = ZZ(tp / N) + tp0
        M = L - p * quadratic_norm(λ * C + zp * N, λ * D + tp * N)
        M, check = M.quo_rem(N**2)
        assert check == 0, "Cant divide by N^2"

        if M < 0:
            continue

        # Try and find a solution to
        # M = x^2 + y^2
        two_squares = Cornacchia(ZZ(M), -ZZ(ω**2))
        if two_squares:
            xp, yp = two_squares
            fz,ft = (λ * C + N * zp) ,(λ * D + N * tp)
            if(swap):
                fz,ft = ft, fz
            if((ZZ(ft+xp)%2==0)and(ZZ(fz+yp)%2==0)):
                break

    if xp is None:
        # Never found vector which had a valid solution
        return None

    # Use solution to construct element μ
    # μ = λ*j*(C + D*ω) + N*(xp + ω*yp + j*(zp + ω*tp))

    # If we swapped earlier, swap again!
    if swap:
        C, D = D, C

    μ = N * xp + N * yp * ω + (λ * C + N * zp) * j + (λ * D + N * tp) * j * ω
    assert ZZ(μ.reduced_norm())==L
    μ = μ/2

    # Check that Nrd(μ) == L
    # and that μ is in O0
    assert ZZ(μ.reduced_norm())%4==0
    assert μ.reduced_norm() == L//4
    assert μ in O0
    return μ


#assumes Ic and It of prime norm and multiplyable
def SigningBorelKLPT(I0, Iτ0, Itor0, Ntor_factor):
    assert is_cyclic(I0), "I is not cyclic"
    assert is_cyclic(Iτ0), "Iτ is not cyclic"
    assert is_cyclic(Itor0), "Itor is not cyclic"


    Ntor=ZZ(Itor0.norm())
    assert(prod([a**b for (a,b) in Ntor_factor])==Ntor)
    #Adapt inputs
    Iτ1,alpha0 = Borel_small_equivalent_ideal(Iτ0)
    Iτ,_,alpha = EquivalentPrimeIdealConstraint(Iτ1,Ntor)
    ac0 = alpha0.conjugate()/(Iτ0.norm())
    ac1 = alpha.conjugate()/(Iτ1.norm())
    ac = ac0*ac1
    I1 =  (1/ac)*I0*ac
    Itor = (1/ac)*Itor0*ac
    I,_,alphaI = EquivalentPrimeIdealfromLattice(I1,I1.intersection(eichler_order_from_ideal(Itor*1)))
    aI = alphaI.conjugate()/(I1.norm())

    assert Iτ.left_order() == O0

    # Prime norm ideals
    Ntor=ZZ(Itor.norm())
    Nτ = ZZ(Iτ.norm())
    Nτtor=ZZ((Iτ*Itor).norm())
    N=ZZ(I.norm())
    assert(Nτ.is_prime())
   # assert(Ntor.is_prime())
    assert(N.is_prime())

    # Orders needed for pushback and pullforward
    O1 = I.left_order()
    assert Iτ.right_order() == O1

    # Prime norm ideal
    Ntor  = ZZ(Itor.norm())
    Nτ    = ZZ(Iτ.norm())
    Nτtor = ZZ((Iτ*Itor).norm())
    N     = ZZ(I.norm())
    assert Nτ.is_prime()
    #assert Ntor.is_prime()
    assert N.is_prime()


    # Compute the pullback K of I with left order O0, and
    # find an equivalent prime norm ideal L ~ I (no torsion constraint).
    L, N, δ = BorelDerive_L(I, Iτ, Nτ, Itor, Ntor,O0, O1) 
    assert L.left_order()==O0
    assert N.is_prime()

    # Set e1
    e1 = floor(logp - log(N, l) + 1.74 * loglogp)+4 
    L1 = l**e1

    # Compute all Eichler orders
    EichlerIτ = eichler_order_from_ideal(Iτ)
    EichlerIτtor = eichler_order_from_ideal(Iτ*Itor)
    EichlerItor = eichler_order_from_ideal(Itor)
    assert(N.is_prime())
    assert(Nτ.is_prime())
    #assert(Ntor.is_prime())

    seen_γ = set()

    #main Loop adapted from LearningToSQI
    for _ in range(8000):
        γ = FullRepresentIntegerHeuristic(N * L1)
        if γ is None:
            print(f"DEBUG [SigningBorelKLPT]: Unable to compute a γ, trying again.")
            continue

        if γ in seen_γ:
            print(f"DEBUG [SigningBorelKLPT]: Already tried γ, trying again.")
            continue
        seen_γ.add(γ)

        # If this GCD is non-trivial, EichlerModConstraint will break
        if gcd(γ.reduced_norm(), Nτtor) != 1:
            continue

        # Set e2
        L2 = derive_L2_BorelSigningKLPT(γ, L1, e1)
        C, D = derive_C_and_D_SigningBorelKLPT(L, N, Iτ, Nτ, EichlerIτ, Itor, Ntor, Ntor_factor, EichlerItor, γ, δ, L2)
        if C is None:
            print(f"DEBUG [SigningBorelKLPT]: Unable to compute a C,D, given γ.")
            continue
        # Get λ = sqrt(p*(C**2+D**2)*L2^(-1)) mod Ntor*Nτ*N
        K = Zmod(N*Nτtor)
        try:
            λλ = K(L2) / K(p * quadratic_norm(C, D))
        except Exception as e:
            # p*(C^2 + D^2) is not invertible mod N
            print(f"ERROR [SigningBorelKLPT]: lambda failure e {e}")
            return None
        λs = []
        for x in [N,Nτ]+[a**b for a,b in Ntor_factor]:
            λs.append(int(Zmod(x)(λλ).sqrt()))
        λ = CRT(λs, [N,Nτ]+[a**b for a,b in Ntor_factor])
        if λ == 0:
            return None
        
        #StrongApproximaztion using full variant
        μ = full_strong_approximation_lattice_heuristic(N * Nτtor, C, D,λ, L2,False)

        # No solution found, try another γ
        if not μ:
            print(f"DEBUG [SigningBorelKLPT]: StrongApproximation failed!", N, Nτ, Ntor, p, L2)
            continue

        # Strong approximation norm check
        assert μ.reduced_norm() == L2
        print(f"INFO [SigningBorelKLPT]: Found a solution to the StrongApproximation!")

        # Now construct the equivalent ideal J
        β = γ * μ
        J_prime = chi(β, L)

        # J = [Iτ]_* J'
        J = pushforward_ideal(O0, O1, J_prime, Iτ)


        # Check things are actually equivalent
        assert equivalent_left_ideals(I, J)
        assert equivalent_left_ideals(I0, ac*J*(1/ac))
        # Check torsion is right
        assert β*(δ) in EichlerIτtor
        assert ac*β*(δ)*(1/ac) in eichler_order_from_ideal(Itor0)

        # Make sure the output is cyclic
        J, _ = make_cyclic(J)

        # S far we do not enforce to hit a specific power of  bu are happy with any.
        print(f"INFO [SigningBorelKLPT]: Preliminary solution of norm {factor(J.norm()) = }")
        #if J.norm() != l**e:
        #    print(f"DEBUG [SigningBorelKLPT]: J has the wrong norm, trying again")
        #    print(f"DEBUG [SigningBorelKLPT]: {factor(J.norm()) = }")
        #    continue

        print(f"INFO [SigningBorelKLPT]: Returning solution of norm {factor(J.norm()) = }")
        return ac*J*(1/ac)
