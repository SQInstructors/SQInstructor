
from KLPT import RepresentIntegerHeuristic
from deuring import *
from setup import *
from isogenies import EllipticCurveIsogenyFactored, generate_random_point,generate_point_order_D
from BorelKLPT import SigningBorelKLPT
from ScalarKLPT import SigningScalarKLPT
from utilities import inert_prime
from ideals import pushforward_ideal

def gen_isogeny_power_l():
    """
    Procedure generating an isogeny power of l, as of now it uses the efficient 
    keygen as described in Appendix D of the SQISign paper

    Input: None
    Output: 
        - EA: the codomain of the isogeny τ_prime
        - τ_prime: the secret isogeny from E0 → EA
        - Jτ: ideal with smooth norm, equivalent to Iτ and to
                τ_prime under the Deuring correspondence
        - Iτ: ideal with prime norm equivalent to Jτ
    """
    # Compute a random prime ≤ Bτ which is inert
    # in R[ω].
    # Note: this is the same as picking p ≡ 3 mod 4
    Nl = l**eτ
    for _ in range(1000):
        Nτ = inert_prime(ceil(sqrt(p)), -ZZ(ω**2))
        # We need the product to be large enough for
        # RepresentIntegerHeuristic.
        if Nτ * Nl > 2 * p:
            break
    assert Nτ * Nl > 2 * p
    # Compute an endomorphism γ of norm Nτ l^eτ
    # Nτ < Bτ
    γ = None

    # Stop infinite loops
    for _ in range(1000):
        γ = RepresentIntegerHeuristic(Nτ * Nl, parity=True)
        if γ is not None:
            break

    if γ is None:
        exit("Never found an alg element with norm (Nτ * Nl), Exiting...")

    Iτ = O0 * γ + O0 * Nτ
    Jτ = O0 * γ.conjugate() + O0 * Nl

    # Iτ has prime norm
    assert Iτ.norm() == Nτ, f"{Iτ.norm() = }, {Nτ = }, {Nl = }"
    # Jτ has smooth norm l^e
    assert Jτ.norm() == Nl, f"{Jτ.norm() = }, {Nτ = }, {Nl = }"

    # Iτ is an integral ideal: Iτ ⊂ O0
    # assert is_integral(Iτ), "Iτ is not integral"

    # Jτ is an integral ideal: Iτ ⊂ O0
    # assert is_integral(Jτ), "Jτ is not integral"

    # Jτ is a cyclic isogeny
    # assert is_cyclic(Jτ), "Jτ is not cyclic"

    # Compute the secret isogeny τ
    I_trivial = O0.unit_ideal()
    ϕ_trivial = E0.isogeny(E0(0))
    τ_prime = IdealToIsogenyFromKLPT(
        Jτ, I_trivial, ϕ_trivial, end_close_to_E0=True
    )
    EA = τ_prime.codomain()

    # The isogeny τ_prime should have degree = n(Jτ)
    assert (
        τ_prime.degree() == Jτ.norm()
    ), f"{factor(τ_prime.degree()) = } {factor(Jτ.norm()) = }"

    EA.set_order((p**2 - 1) ** 2)
    return EA, τ_prime, Jτ, Iτ

def setup():
    """
    Helper for tests
    """
    Tf = list(T.factor())

    prod1, prod2 = 1,1
    prod1 = 2**((E0.order()).valuation(2)-5)
    for m,_ in Tf[:len(Tf)//4]:
        prod2 = prod2*m
    print("prods done")
    # Pc = generate_point_order_D(E0,prod1)
    # print("Pc done")
    # Cphi = EllipticCurveIsogenyFactored(E0,Pc)
    E1, Cphi, Ic, Iτ = gen_isogeny_power_l()
    # this assertion fails, should we care?
    # assert Iτ.right_order() == Ic.right_order()
    print("Cphi done")
    # print(Cphi)
    # print(Cphi.degree().factor())
    Cphid = Cphi.dual()
    print("Dual done")
    # Ic = kernel_to_ideal(Pc,prod1)
    PI_0 = generate_point_order_D(E0,prod2)
    # PI = generate_point_order_D(E1,prod2)
    PI = Cphi(PI_0)
    Iphi = EllipticCurveIsogenyFactored(E1,PI)
    # print(Iphi)
    print("Iphi done")
    I_0 = kernel_to_ideal(PI_0, prod2)
    I = pushforward_ideal(O0, Ic.right_order(), I_0, Ic)
    assert I.left_order() == Ic.right_order()
    print("I done")
    print('We still have to test that I is the ideal corresponding to Iphi')

    #Ntor = sqrt(p).round().next_prime()
    #PI = generate_point_order_D(E1,Ntor)
    Itor = 0
    Ntor = list(E0.order().factor())[-2][0]
    Ntor_factor = [(Ntor,1)]
    print(Ntor)
    Ptor0 = generate_point_order_D(E0,Ntor)
    Ptor = Cphi(Ptor0)
    PtorIm = Iphi(Ptor)

    Itor0 = kernel_to_ideal(Ptor0, Ntor)
    Itor = pushforward_ideal(O0, Ic.right_order(), Itor0, Ic)


    assert Itor.left_order() == Ic.right_order()
    assert Itor.left_order() == I.left_order()

    #assert(gcd(ZZ(Itor.norm()),ZZ(I.norm())==1))
    #assert(gcd(ZZ(Itor.norm()),ZZ(Ic.norm())==1))
    #assert(gcd(ZZ(Ic.norm()),ZZ(I.norm())==1))
    return I,Ic,Itor,Cphi,Iphi,Ptor, PtorIm, Ntor_factor

def id2iso_borel_test():
    I,Ic,Itor,Cphi,Iphi,Ptor, PtorIm, Ntor_factor = setup()
    print("setup done")

    J = SigningBorelKLPT(I,Ic,Itor, Ntor_factor)
    print("KLPT done")
    print(log(J.norm(),2))
    print(J.norm()%2**log(J.norm(),2))
    phiJ = IdealToIsogenyFromKLPT(J, Ic, Cphi)
    print("id2iso done")
    PtorImJ = phiJ(Ptor)
    print("eval done")
    EI = Iphi.codomain()
    EJ = phiJ.codomain()
    assert EI.is_isomorphic(EJ)
    iso = EI.isomorphism_to(EJ)
    PtorImI = iso(PtorIm)
    return PtorImI.weil_pairing(PtorImJ,ZZ(Itor.norm())).is_one()


def id2iso_scalar_test():
    I,Ic,Itor,Cphi,Iphi,Ptor, PtorIm, Ntor_factor = setup()
    print("setup done")
    Ntor = Ntor_factor[0][0]
    Qtor0 = generate_point_order_D(E0,Ntor)
    Qtor = Cphi(Qtor0)
    while (Ptor.weil_pairing(Qtor,ZZ(Itor.norm())).is_one()):
        Qtor0 = generate_point_order_D(E0,Ntor)
        Qtor = Cphi(Qtor0)
    print("points ok\n")

    QtorIm = Iphi(Qtor)

    J = SigningScalarKLPT(I,Ic,Ntor, Ntor_factor,True)
    print("KLPT done")
    print(log(J.norm(),2))
    print(J.norm()%2**log(J.norm(),2))
    phiJ = IdealToIsogenyFromKLPT(J, Ic, Cphi)
    print("id2iso done")
    PtorImJ = phiJ(Ptor)
    QtorImJ = phiJ(Qtor)
    print("eval done")
    EI = Iphi.codomain()
    EJ = phiJ.codomain()
    assert EI.is_isomorphic(EJ)
    iso = EI.isomorphism_to(EJ)
    PtorImI = iso(PtorIm)
    QtorImI = iso(QtorIm)
    #return PtorImI.weil_pairing(PtorImJ,ZZ(Itor.norm())).is_one()
    try:
        a = PtorImI.log(PtorImJ)
        b = QtorImI.log(QtorImJ)
        print("diag")
        return(a%Ntor==b%Ntor)
    except ValueError:
        print("wrong")
        return(False)

if __name__=="__main__":
    print(id2iso_scalar_test())
    print(id2iso_borel_test())
