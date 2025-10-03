#!/usr/bin/env python3
import sample

def random_scalar_ideal(O0, N):
    assert N % 257 != 0
    p = O0.discriminant()
    B = O0.quaternion_algebra()
    O = sample.random_ideal(O0).right_order()
    I = sample.random_ideal(O, p**2 * N**3)
    d = I.norm()
    IN = B.ideal( (N*I).basis() + (d,) )
    assert IN.norm() == d
    assert not IN.left_order().is_maximal()
    assert not IN.right_order().is_maximal()
    return IN

def minima(I, N, normalize=1):
    basis = I.reduced_basis()
    normalize = I.norm() * normalize**(1/2)
    minima = sorted([g.reduced_norm() for g in basis])
    parity = [n % N for n in minima]
    normal = [n / normalize for n in minima]
    return normal, parity

if __name__ == '__main__':
    import sys
    # args are either a prime characteristic and a level
    p = int(sys.argv[1])
    O0 = sample.start_order(p)
    N = sample.ZZ(sys.argv[2])

    # Print forever the minima of random left O-ideas classes and
    # their right orders, scaled down by factors of p^½ and p^⅔
    # respectively
    while True:
        I = random_scalar_ideal(O0, N)
        nor, par = minima(I, N, p * N**3)
        print(','.join(f'{float(m):17.15f}' for m in nor)
              + ','
              + ','.join(str(p) for p in par),
              flush=True)
