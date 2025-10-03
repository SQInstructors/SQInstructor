#!/usr/bin/env python3
from sage.all import *
proof.all(False)
set_random_seed()

def serialize(O):
    s = repr([O.discriminant(), list(O.basis_matrix())])
    assert parse(s) == O
    return s

def parse(Ostring):
    p, basis = sage_eval(Ostring)
    if not is_prime(p):
        raise ValueError(f'{p:#x} is not prime')
    B = QuaternionAlgebra(p)
    O = B.maximal_order(order_basis=tuple(B(g) for g in basis))
    if p != O.discriminant():
        raise ValueError(f'{O} is not maximal')
    return O

def random_ideal(O0, deg=None):
    p = O0.discriminant()
    if deg is None:
        deg = p**2
    l = 257   # or something
    steps = ceil(log(deg, l))
    I = O0*1
    for t in reversed(range(steps)):
        O = I.right_order()
        while True:
            α = sum(randrange(l)*g for g in O.basis())
            if α.reduced_norm().valuation(l) == 1:
                break
        assert α in O
        J = O*l + O*α
        assert J.left_order() == O
        assert J.norm() == l
        I *= J
    #γ = I.minimal_element()
    #I *= γ.conjugate() / I.norm()
    return I

def random_ideals(O):
    print(f'{O = }', file=sys.stderr)

    while True:
        yield random_ideal(O)

def start_order(p):
    A = QuaternionAlgebra(p)
    O0 = A.maximal_order()
    return O0

def random_order(p):
    O0 = start_order(p)
    O = random_ideal(O0).right_order()
    return O

def minima(I, normalize=1):
    basis = I.reduced_basis()
    Obasis = (I.right_order() * 1).reduced_basis()
    minima = sorted([g.reduced_norm() / I.norm() / normalize**(1/2)
                     for g in basis])
    Ominima = sorted([g.reduced_norm() / normalize**(2/3)
                      for g in Obasis])[1:]
    return (minima, Ominima)

if __name__ == '__main__':
    # args are either a prime characteristic,
    # or the serialization of a maximal order
    try:
        p = int(sys.argv[1])
        O = start_order(p)
    except ValueError:
        O = parse(sys.argv[1])
        p = O.discriminant()

    # Print forever the minima of random left O-ideas classes and
    # their right orders, scaled down by factors of p^½ and p^⅔
    # respectively
    for n,I in enumerate(random_ideals(O)):
        m, Om = minima(I, p)
        s = ','.join(f'{float(m):17.15f}' for m in m + Om)
        print(s, flush=True)
