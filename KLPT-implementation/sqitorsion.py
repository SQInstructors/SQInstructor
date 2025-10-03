from sage.all import *
from KLPT import RepresentIntegerHeuristic
from deuring import *
from setup import *
from isogenies import EllipticCurveIsogenyFactored, generate_random_point,generate_point_order_D
from sage.schemes.elliptic_curves.weierstrass_morphism import WeierstrassIsomorphism
from BorelKLPT import SigningBorelKLPT
from ScalarKLPT import SigningScalarKLPT
from utilities import inert_prime
from ideals import pushforward_ideal
from torsion import equivalentTorsionPrimeIdeal, torsion_matrices_basis

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

def gen_isogeny_coprime(T_prime = T_prime):
    """
    Compute an isogeny starting from E0 and corresponding ideal
    of degree / norm T_prime

    Input: None
    Output: 
        - E1: the codomain of the commitment isogeny
        - ψ_ker: the kernel of ψ
        - ψ: the secret commitment isogeny ψ : E0 → E1
        - Iψ: the ideal equivalent to ψ.
    """
    # Generate a random kernel
    # of order T_prime
    P, Q = torsion_basis(E0, T_prime)
    x = randint(1, T_prime)
    ψ_ker = P + x * Q

    # Generate the ideal Iψ from ψ_ker
    Iψ = kernel_to_ideal(ψ_ker, T_prime)
    assert Iψ.norm() == T_prime, "Iψ has the wrong norm"

    # Generate the ideal ψ and its codomain
    ψ = EllipticCurveIsogenyFactored(E0, ψ_ker, order=T_prime)
    E1 = ψ.codomain()

    E1.set_order((p**2 - 1) ** 2)
    return E1, ψ_ker, ψ, Iψ


class SQItorsion():

    def __init__(self, Scalar=False):
        """
        We only use the object to store intermediate values,
        we could also init this with parameters and stop using
        globals.
        """
        # Key pair
        # pk = EA
        # sk = (τ_prime, Iτ, Jτ)
        self.pk = None
        self.sk = None

        # Secret commitment values
        # commitment_secrets = (ψ_ker, ψ, Iψ)
        self.commitment_secrets = None

        self.Ntor = list(E0.order().factor())[-1][0]
        self.Ntor_factor = [(self.Ntor,1)]
        # Take into account current limitation of ScalarKLPT
        if(Scalar):
            cntr = 1
            while(Zmod(self.Ntor)(2).is_square()):
                cntr = cntr+1
                self.Ntor = list(E0.order().factor())[-cntr][0]
                self.Ntor_factor = [(self.Ntor,1)]
        self.E0basis = torsion_basis(E0, self.Ntor, canonical=True)
        self.action_matrices = torsion_matrices_basis(self.E0basis[0],
                                                      self.E0basis[1],
                                                      self.Ntor)

        prod2 = 1
        Tf = list(T.factor())
        for m,_ in Tf[:len(Tf)//4]:
            prod2 = prod2*m
        self.T_commit = prod2

    def keygen(self, debug = False):
        """
        Compute the public and secret keys.

        Stores to self:
            - self.pk = EA: the public key
    
            - self.sk = Iτ, τ, τ_ker
                Iτ: the ideal equivalent to τ of norm Tprime.
                τ: the secret commitment isogeny τ : E0 → E1
                τ_ker: the kernel of τ
        """
        if debug:
            print("Generating trivial keypair")
            EA = E0
            τ_prime = E0.isogeny(E0(0))
            Iτ = O0.unit_ideal()
            Jτ = O0.unit_ideal()
            self.pk = EA
            self.sk = (τ_prime, Iτ, Jτ)
            self.pkbasis = (self.E0basis[0], self.E0basis[1])

            N = self.Ntor
            Ptor0 = self.E0basis[0]
            self.pkpoint = Ptor0

            Itor0 = kernel_to_ideal(Ptor0, N)
            self.skpoint_ideal = Itor0
            return None
        EA, τ_prime, Jτ, Iτ = gen_isogeny_power_l()
        assert Jτ.norm().valuation(2) >= 10, f"Jτ norm is not divisible enough by 2 {Jτ.norm()}"
        self.pk = EA
        self.sk = (τ_prime, Jτ, Iτ)
        # self.pkbasis = torsion_basis(self.pk, self.Ntor, canonical=True)
        self.pkbasis = (τ_prime(self.E0basis[0]), τ_prime(self.E0basis[1]))
        N = self.Ntor
        Ptor0 = self.E0basis[0]
        Ptor = τ_prime(Ptor0)
        assert Ptor == self.pkbasis[0]
        self.pkpoint = Ptor

        Itor0 = kernel_to_ideal(Ptor0, N)
        Itor = pushforward_ideal(O0, Jτ.right_order(), Itor0, Jτ)
        self.skpoint_ideal = Itor

        return None

    def export_public_key(self):
        """
        Helper function to return the public key
        """
        if self.pk is None:
            raise ValueError(f"Must first generate a keypair with `self.keygen()`")
        return self.pk

    def commitment(self):
        """
        Compute the commitment isogeny and curve

        Input: None
        Output: E1: the codomain of the commitment isogeny

        Stores to self:

        self.commitment_secrets = (ψ_ker, ψ, Iψ))
            ψ_ker: the kernel of ψ
            ψ: the secret commitment isogeny ψ : E0 → E1
            Iψ: the ideal equivalent to ψ.


        self.cmtbasis = ( P1, Q1 )
            the image of the basis of E0[N] under ψ

        self.commit_curve = E1
            the commitment curve E1

        """
        E1, ψ_ker, ψ, Iψ = gen_isogeny_coprime(T_prime = self.T_commit)
        self.commitment_secrets = ψ_ker, ψ, Iψ
        self.cmtbasis = ( ψ(self.E0basis[0]), ψ(self.E0basis[1]) )

        self.commit_curve = E1
        return E1

    def response(self, challenge):
        """
        Given a challenge, compute the response isogeny

        Input:
            - challenge: the challenge, either a matrix in GL(2,Z/NZ)
                or an integer in Z/NZ depending on the subclass

        Output:
            - phiJ: the isogeny ɸ:E_pk -> E_com
            of smooth norm acting as M on the N-torsion.
        """
        raise NotImplementedError("This is implemented in subclasses")

    def verify(self, challenge, phiJ):
        """
        Verify the response isogeny

        Input:
            - challenge: the challenge, either a matrix in GL(2,Z/NZ)
                or an integer in Z/NZ depending on the subclass
            - phiJ: the isogeny ɸ:E_pk -> E_com of smooth norm
                acting as M on the N-torsion.
        Output:
            - True if the isogeny acts as M on the N-torsion, False otherwise
        """
        raise NotImplementedError("This is implemented in subclasses")

class SQIscalar(SQItorsion):
    """
    Wrapper class to use the Scalar KLPT.
    """
    def __init__(self):
        super().__init__(Scalar=True)
        

    def challenge(self, debug = True):
        """
        Generate a random challenge matrix M in GL(2,Z/NZ)
        Input: None
        Output: M: a random matrix in GL(2,Z/NZ)
        """

        N = self.Ntor
        ZN = Zmod(N)
        if debug:
            print(f"Using trivial challenge matrix in GL(2,Z/{N}Z)")
            return matrix(ZN, [[1,0],[0,1]])
        invertible = False
        while not invertible:
            M = matrix(ZN,
                       [
                           [ZN.random_element(),ZN.random_element()],
                           [ZN.random_element(),ZN.random_element()]])
            invertible = M.is_invertible()
        return M

    def _equivalent_correct_ideal(self, M, debug = True):
        """
        Given an integer m ideal I
        equivalent to Iψ * Iτ and such that the action of
        I on E0[N] is given by .

        Input:
            - M: a matrix in GL(2,Z/NZ)

        Output:
            - I: an ideal equivalent to Iψ * Iτ and such that
                the action of I on E0[N] is given by M.
        """
        if self.sk is None or self.commitment_secrets is None:
            raise ValueError(f"Must first generate a keypair with `self.keygen()` and a commitment with `self.commitment()`")
        τ_prime, Iτ, Jτ = self.sk
        ψ_ker, ψ, Iψ = self.commitment_secrets
        assert Iτ.left_order() == Iψ.left_order(), "Ideals are not compatible"
        if debug:
            print("DEBUG: skipping equivalent ideal computation")
            assert M.is_one()
            if Iτ.norm().is_one():
                print("DEBUG: Iτ is trivial")
                return Iψ
            return Iτ.conjugate()*Iψ
        # Compute the ideal I equivalent to Iψ * Iτ, only implemented for Scalar
        I = equivalentTorsionPrimeIdeal(Iτ.conjugate()*Iψ, M, self.Ntor,
                                        base = self.pkbasis,
                                        A = self.action_matrices)
        assert I.norm() in ZZ
        # check that the ideal starts at the pk
        assert I.left_order() == Iτ.right_order()
        # check that the ideal ends at the commitment
        # assert I.right_order() == Iψ.right_order()
        return I

    def _equivalent_smooth_isogeny_scalar(self, I):
        """
        Given an input ideal I find an equivalent ideal J of smooth norm and 
        find the associated isogeny
        """
        τ_prime, Iτ, Jτ = self.sk
        # breakpoint()
        print(f"Computing scalar KLPT with connecting ideal of norm {Iτ.norm().factor()}")
        assert I.left_order() == Iτ.right_order(), "Ideals are not compatible"
        J = SigningScalarKLPT(I, Iτ, self.Ntor, list(self.Ntor_factor))
        assert J.norm().valuation(2) >= 100
        print(f"Found ideal J of norm {J.norm().factor()}")
        phiJ = IdealToIsogenyFromKLPT(J, Iτ, τ_prime)
        return phiJ, J

    def response(self, challenge):
        """
        Given a challenge matrix M in GL(2,Z/NZ), compute an ideal I
        equivalent to Iψ * Iτ and such that the action of
        I on E0[N] is given by M.

        Input:
            - M: a matrix in GL(2,Z/NZ)

        Output:
            - phiJ: the isogeny ɸ:E_pk -> E_com
            of smooth norm acting as M on the N-torsion.
        """
        if self.sk is None or self.commitment_secrets is None:
            raise ValueError(f"Must first generate a keypair with `self.keygen()` and a commitment with `self.commitment()`")
        I = self._equivalent_correct_ideal(challenge)
        phiJ, _ = self._equivalent_smooth_isogeny_scalar(I)
        return phiJ

    def verify(self, M, phiJ):
        """
        Verify that the isogeny phiJ acts as M on the N-torsion

        Input:
            - M: a matrix in GL(2,Z/NZ)
            - phiJ: the isogeny ɸ:E_pk -> E_com
            of smooth norm acting as M on the N-torsion.

        Output:
            - True if the isogeny acts as M on the N-torsion, False otherwise
        """
        if self.pk is None or self.commit_curve is None or self.cmtbasis is None:
            raise ValueError(f"Must first generate a keypair with `self.keygen()` and a commitment with `self.commitment()`")
        E1 = self.commit_curve
        Epk = self.pk
        if not phiJ.domain().is_isomorphic(Epk):
            raise ValueError(f"The domain of the response isogeny must be the public key")
        if not phiJ.codomain() == E1:
            if not phiJ.codomain().is_isomorphic(E1):
                raise ValueError(f"The codomain of the response isogeny must be the commitment")
            else:
                iso = phiJ.codomain().isomorphism_to(E1)
                old_phiJ = phiJ
                # phiJ = iso * phiJ
                phiJ = WeierstrassIsomorphism(phiJ.codomain(), None, E1) * phiJ

        Pk1, Pk2 = self.pkbasis
        Cmt1, Cmt2 = self.cmtbasis
        N = self.Ntor
        ZN = Zmod(N)
        # Compute the image of the basis of E0[N] under phiJ
        Im1 = phiJ(Pk1)
        Im2 = phiJ(Pk2)
        from torsion import TorsionMatrix
        MM = TorsionMatrix(Im1, Im2, Cmt1, Cmt2, N)
        MM = MM.inverse()
        print(MM)
        Ctm1M = ZN(M[0,0])*Cmt1 + ZN(M[0,1])*Cmt2
        Ctm2M = ZN(M[1,0])*Cmt1 + ZN(M[1,1])*Cmt2
        out = True
        if not Im1.weil_pairing(Ctm1M, N, algorithm = 'pari').is_one():
            out = False
        if not Im2.weil_pairing(Ctm2M, N, algorithm = 'pari').is_one():
            out = False
        # breakpoint()
        return out


class SQIborel(SQItorsion):
    """
    Wrapper class to use the Borel KLPT instead of the Scalar KLPT
    """

    def __init__(self):
        super().__init__()

    def challenge(self, debug = False):
        """
        Generate a random challenge of the form P + m Q

        Input: None
        Output: m: a random element in Z/NZ
        """
        N = self.Ntor
        ZN = Zmod(N)
        if debug:
            return ZN(0)
        return ZN.random_element()

    def _equivalent_correct_ideal(self, challenge, debug = False):
        """
        Given an integer challenge compute an ideal I
        equivalent to Iψ * Iτ and such that the action of
        I on τ(P0 + m Q0) is Ψ(P0) and on τ(Q0) is Ψ(Q0).

        Input:
            - challenge: an integer in Z/NZ

        Output:
            - I: an ideal equivalent to Iψ * Iτ and such that
                the action of I on E0[N] is given by M.
        """
        if self.sk is None or self.commitment_secrets is None:
            raise ValueError((
                f"Must first generate a keypair with `self.keygen()`"
                "and a commitment with `self.commitment()`"))
        τ_prime, Iτ, Jτ = self.sk
        ψ_ker, ψ, Iψ = self.commitment_secrets
        assert Iτ.left_order() == Iψ.left_order(), "Ideals are not compatible"

        # finding ideal that maps P + m Q to P and Q to Q
        M = matrix(Zmod(self.Ntor), [[1,challenge],[0,1]])
        if debug:
            print("DEBUG: skipping equivalent ideal computation")
            assert M.is_one()
            if Iτ.norm().is_one():
                print("DEBUG: Iτ is trivial")
                return Iψ
            return Iτ.conjugate()*Iψ
        # Compute the ideal I equivalent to Iψ * Iτ, only implemented for Scalar
        I = equivalentTorsionPrimeIdeal(Iτ.conjugate()*Iψ, M, self.Ntor,
                                        base = self.pkbasis,
                                        A = self.action_matrices)
        assert I.norm() in ZZ
        # check that the ideal starts at the pk
        assert I.left_order() == Iτ.right_order()
        # check that the ideal ends at the commitment
        # assert I.right_order() == Iψ.right_order()
        return I

    def _equivalent_smooth_isogeny_borel(self, I, challenge):
        """
        Given an input ideal I find an equivalent ideal J of smooth norm and 
        find the associated isogeny
        """
        τ_prime, Iτ, _ = self.sk
        # Itor = self.skpoint_ideal
        Ppk, Qpk = self.pkbasis

        # Ptor is the point to be preserved
        Ptor = Ppk + challenge*Qpk

        Ptor0 = τ_prime.dual()(Ptor)
        Itor0 = kernel_to_ideal(Ptor0, self.Ntor)
        Itor  = pushforward_ideal(O0, Iτ.right_order(), Itor0, Iτ)
        assert Itor.left_order() == Iτ.right_order()

        print(f"Computing borel KLPT with connecting ideal of norm {Iτ.norm().factor()}")
        assert I.left_order() == Iτ.right_order(), "Ideals are not compatible"
        J    = SigningBorelKLPT(I, Iτ, Itor, self.Ntor_factor)
        assert J.norm().valuation(2) >= 100
        print(f"Found ideal J of norm {J.norm().factor()}")
        phiJ = IdealToIsogenyFromKLPT(J, Iτ, τ_prime)
        return phiJ, J

    def response(self, challenge):
        """
        Given a challenge integer m \in Z/NZ, compute an ideal I
        equivalent to Iψ * Iτ and such that phi applied on
        P0 + m Q0 give the committed point, for now Ψ(P0).

        Input:
            - m: an integer in Z/NZ

        Output:
            - phiJ: the isogeny ɸ:E_pk -> E_com
            of smooth norm acting as M on the N-torsion.
        """
        if self.sk is None or self.commitment_secrets is None:
            raise ValueError(
                    (f"Must first generate a keypair with"
                     "`self.keygen()` and a commitment with "
                     "`self.commitment()`"))
        I = self._equivalent_correct_ideal(challenge)
        phiJ, _ = self._equivalent_smooth_isogeny_borel(I, challenge)
        return phiJ

    def verify(self, m, phiJ, debug = False):
        """
        Verify that the isogeny phiJ acts as M on the N-torsion

        Input:
            - M: a matrix in GL(2,Z/NZ)
            - phiJ: the isogeny ɸ:E_pk -> E_com
            of smooth norm acting as M on the N-torsion.

        Output:
            - True if the isogeny acts as M on the N-torsion, False otherwise
        """
        if self.pk is None or self.commit_curve is None or self.cmtbasis is None:
            raise ValueError(f"Must first generate a keypair with `self.keygen()` and a commitment with `self.commitment()`")
        if m.is_zero():
            print("WARNING: m is zero, most likely this is a debug run")
        E1 = self.commit_curve
        Epk = self.pk
        if not phiJ.domain().is_isomorphic(Epk):
            raise ValueError(f"The domain of the response isogeny must be the public key")
        if not phiJ.codomain() == E1:
            if not phiJ.codomain().is_isomorphic(E1):
                raise ValueError(f"The codomain of the response isogeny must be the commitment")
            else:
                print("[WARNING] Isomorphism needed to link the codomains")
                phiJ = WeierstrassIsomorphism(phiJ.codomain(), None, E1) * phiJ

        Pk1, Pk2 = self.pkbasis
        Cmt1, Cmt2 = self.cmtbasis
        N = self.Ntor
        ZN = Zmod(N)
        # Compute the image of the basis of E0[N] under phiJ
        Im1 = phiJ(Pk1 + m*Pk2)
        out = True
        if not Im1.weil_pairing(Cmt1, N, algorithm = 'pari').is_one():
            out = False
        if debug:
            if not out:
                print("DEBUG: verification failed, entering debugger")
            print(phiJ(Pk1).weil_pairing(Cmt1 + m*Cmt2, N, algorithm = 'pari'))
            # breakpoint()
        return out

