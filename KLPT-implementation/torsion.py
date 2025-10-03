from sage.all import vector,  matrix, ZZ, identity_matrix, Zmod, is_prime, QQ, log
from isogenies import torsion_basis
from helpers import ideal_add
from pari_interface import discrete_log_pari
from setup import *
from deuring import E01, E0ι, E0π, E0ιπ


def torsion_matrices_basis(P0, Q0, N):
    '''
    Return the torsion matrix associated to the endomorphisms
    1, ɩ, π and ɩπ, definined on the quaternion algebra B(-1,-p)
    and acting on points in E0 (j = 1728)
    '''
    if P0.curve().j_invariant() == 1728:
        print('Attention points defined on a different curve')
    assert P0.curve().base_field().characteristic() % 4 == 3, 'characteristic is not 3 mod 4'
    assert P0.curve().cardinality() % (N**2) == 0, 'torsion points not avaible' 
    A1  = identity_matrix(Zmod(N),2)
    Aπ  = TorsionMatrix(E0π(P0), E0π(Q0),P0,Q0,N)
    Aɩ  = TorsionMatrix(E0ι(P0), E0ι(Q0),P0,Q0,N)
    Aɩπ = TorsionMatrix(E0ιπ(P0), E0ιπ(Q0),P0,Q0,N)
    assert Aɩπ == Aπ*Aɩ, 'composition of matrices not behaving as expected'
    return [A1,Aɩ,Aπ,Aɩπ]


def endomorphism_to_torsion(x, N, A = None, base = None, debug = True):
    '''
    Return the corresponding linear map associated to 
    the endomorphism on the N torsion

    - A are the matrix associated to 1, ɩ, π and ɩπ
    - base is the torsion basis used.

    '''
    assert x in B, 'x is not a quaternion'
    if base:
        P0, Q0 = base
    else:
        P0, Q0 = torsion_basis(E0, N, canonical=True)
    if not A:
        A = torsion_matrices_basis(P0, Q0, N)
    return sum([Zmod(N)(xcoo) * A[idx] for (idx,xcoo) in enumerate(x.coefficient_tuple())])


def endomorphism_with_torsion(I, M, A = None, base = None, debug = True):
    '''
    Return an element y in I such that  the corresponding linear map associated to 
    the endomorphism on the N torsion is M, where N is the cardinality of 
    the base ring of M
    '''
    R = M.base_ring()
    N = R.cardinality()
    # constrain on avaible torsion, in future we may change it to avaible_torsion % N == 0
    assert (p**2 - 1 ) % N == 0, 'N torsion not avaible'
    if base:
        P0, Q0 = base
    else:
        P0, Q0 = torsion_basis(E0, N, canonical=True)
    if A:
        assert all([R == AA.base_ring() for AA in A]), 'input matrix with wrong base ring'
    else:
        A = torsion_matrices_basis(P0, Q0, N)

    if debug:
        print(f'Using torsion basis:\nP0 = {P0}\nQ0 = {Q0}')
        print(f'Using torsion matrices:\n{A[0]}\n{A[1]}\n{A[2]}\n{A[3]}')
        print(f'Searching for endomorphism in ideal with basis:\n{I.basis()}')
        print(f'and torsion matrix:\n{M}')
    base = matrix(R,[(endomorphism_to_torsion(x,N, A = A, base = (P0,Q0)).list()) for x in I.basis()])
    sol = vector(M.list()) / base
    if debug:
        print(f'basis of N-torsion used:\n{base}')
        print(f'solution vector: {sol}')
    # ZZ coercion necessary since quaternions are over QQ
    return sum( [ ZZ(s)*b for (s,b) in zip(sol,I.basis())])


def ideal_to_torsion(I, N, base):
    '''
        return the N-torsion of I in O_0
    '''
    pass


def equivalentTorsionPrimeIdeal(I0, M, N_tor, bound = 100, A = None, base = None, debug = True):
    """
        Let ɸ0 the isogeny associated to I0, if P0,Q0 = base then
        we return an ideal I of prime norm such that the assicated 
        isogeny ɸ satisfies [ɸ(P0),ɸ(Q0)] = M * [ɸ0(P0), ɸ0(Q0)]
    """
    i0 = endomorphism_with_torsion(I0, M, A = A, base = base, debug = debug)
    S0 = ideal_add(I0*N_tor,matrix(QQ,[i0.coefficient_tuple()]))
    assert i0 in S0
    B = S0.basis()
    G = S0.gram_matrix()
    U = G.LLL_gram().transpose()
    M = [sum(c * g for c, g in zip(row, B)) for row in U]
    def size(x):
        return log(x,2).numerical_approx()
    # print([size(m.reduced_norm()) for m in M])
    n_I = I0.norm()
    for x in range(bound):
        for a in range (0,x):
            for b in range (-x+a,x-a):
                for c in range (-x+a+b,x-a-b):
                    for d in range (-x+a+b+c,x-a-b-c):
                        v = vector([d,c,b,a])
                        el = sum( v[i]*M[i] for i in range(4))
                        n = el.reduced_norm()/n_I
                        assert(n.is_integer())
                        n = ZZ(n)
                        if(is_prime(n)):
                            #should assert that el is in I
                            assert el in S0
                            assert el in I0
                            # assert all([(el[i] - i0[i] )% N_tor == 0 for i in [1,2,3]])
                            return I0*(el.conjugate()/n_I)

def TorsionMatrix(R1, R2, P, Q, D):
    '''
    Given a basis P,Q of E[D] finds
    a 2 by 2 matrix A such that 
    (R1, R2)^t = A * (P,Q)^t.
    '''
    A = [ list(BiDLP(R, P, Q, D)) for R in [R1,R2]]
    return matrix(Zmod(D), A)


def BiDLP(R, P, Q, D):
    """
    Given a basis P,Q of E[D] finds
    a,b such that R = [a]P + [b]Q.

    Uses the fact that
        e([a]P + [b]Q, [c]P + [d]Q) = e(P,Q)^(ad-bc)
    """
    # e(P,Q)
    pair_PQ = P.weil_pairing(Q, D, algorithm="pari")

    # Write R = aP + bQ for unknown a,b
    # e(R, Q) = e(P, Q)^a
    pair_a = R.weil_pairing(Q, D, algorithm="pari")

    # e(R,-P) = e(P, Q)^b
    pair_b = R.weil_pairing(-P, D, algorithm="pari")

    # Now solve the dlog on Fq
    a = discrete_log_pari(pair_a, pair_PQ, D)
    b = discrete_log_pari(pair_b, pair_PQ, D)

    return a, b

