#!/usr/bin/env sage
proof.all(False)

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-v', '--verbose', action='store_true')
parser.add_argument('-p', '--parameter', choices=['128', '256'], default='128')
parser.add_argument('-d', '--dlogs', action='store_true', help='Recompute discrete logarithms')
args = parser.parse_args()

if args.verbose:
    print(f'Generating constants for parameter {args.parameter}', file=sys.stderr)

################################################################

if args.parameter == '128':
    ALPHA1 = -0x52fd2eac925dea794e9b5def96e5004fad803ceac43e1f4565cc8f57e442a6cc6
    ALPHA2 =  0x707613852b1c015d265033276a78eef61d5efbe124c530324118e1aea7c5ee713

    NUMBER_OF_PRIMES = 65
    c = 335
    
    A = [0xb34c24addb87baf2dd78ef7454900860624f861940859753002e91b5e0b01ed218230b5467fb67bd24461d74898c052bfab47ca5d33d2b8f2aad21d25d65ac92bfc3, 0xb1a27ca357c0cd93119731939777be6140ca1eaf3c8421395c10dcdb291bdfc4f7da48e5662ad134b72e4792aecb1a97bfba63a1fba087306224ae16367c4fafe3dc]
    xP = [0xada260e4f64847c143cf51deab4dac2ff228f6a9d5ae72aa40280fe3849a01dec65c5a9c726040593d3b771254acf073ecbcabb29a975e4d0a95abf548ff5e2ac967, 0x37841734b835d3f5c7c1b4e54093b4681b5655daaef380df5f21f48912bae6abb49bc7e2a60563ea337c027bd63a4e705e11d819d9276ab22bf0c234c8fb058d1e20]
    xQ = [0xb1fe2a19950bfeee3535d7230640d3ba02cf81802b6aa0817e9103c0fc22236925085e4404be5e7fa178410892f315defa55d7f8728f6b4e3bde0c446633ceeeb0ad, 0x1afdd6d4212ee8b10dc0b4d1557f80015726d62f131d2d003ec7e8848353193da3883143d0c1e88481d3288875ea6348f67043ec5b9a6561c49f4e7df4d4efb7ad39]

elif args.parameter == '256':
    ALPHA1 = -0xcec79b435a1f29171b67f33f94738cdb99ff79ab94cc3a58cf15dfb3ed1529f0a8b2eb48bf9a00e2304a6cedaf56114babf9f699859928cfe3b2274dd0f1333f6
    ALPHA2 = 0x848a5a0bca8547efbf2557524bc9feccb98aca333ecd39619366888446bf8892dd6eb129114a0668457dbbc4aa53da5b1355699fe0d091eed7654ac90fd5a4f5f

    NUMBER_OF_PRIMES = 75
    c = 256
    
    A = [0x2bc99d9a4119dc11a821d4d84dd8b0207a808b287447db5ac00d406c94d06e9bbfb8d46b4ea56b0d1186241bba207987780d4646ab2ac47b0c59349d60298a0895e962c10c6d3b95ca4441d2cebc, 0x3628318aa326979f50f58ca44ed56e3f0e7304b9a05f0f74c2bb559bedf03f38c238ce825ec3d145de3b8d2bd474b6e5a29b0df60ac41b2235e12e69a69cbc120de0cbbff92c0e70ec094085fcfa]
    xP = [0xe3fa137669028467ad42e84baab75c6a55a83024208869949a52523b84be73037f0d85852252a823ac4986615398d15a237bcd7462f2bc8931bdc78dd557cad3fb7f2ef8afb359ff16cd91f575fa, 0x694fce5012d6aa5e245276c6b1bddceb018e2e88de6add76f5f8805cbd139b60c64eeed3f2b5fc64802ff92c83c9ec7e43fd81283e79fc10f4e1d6bc690a7d705ea1d449882329bd5331031c58d5]
    xQ = [0xa70699854803b067d8604b41a460d24c5333d6d9839b047c291f34ccef9faa58835b34cdac255e0c02f2d8c1ba94bb39305aaddfcbaa556745fb4f144efc4ad1052f2edba407ad6885f1985b83bc, 0x11e147b757a35f4c4ef26d2bbb875878a899e50c52210ae4b05b544966e0def2f54de2a7c4acb399c016ba2d67b642e315d38acc87a685aab6709b52030bc0848067877b3b27ade220c736ed7fb7d]

else:
    raise RuntimeError(f"Unknown parameter { args.parameter }")

################################################################


# computes the remaining primes such that L is the product of the n smallest split primes in frakO_0

def completeL(n, L):
    ls = L.prime_factors()
    q = ls[-1]
    while len(ls) < n:
        q = next_prime(q)
        if q % 4 == 1:     # in our case, q splits iff q = 1 mod 4
            ls.append(q)
    return prod(ls)


######## constructing the parameters from alpha

NORM_ALPHA = ALPHA1^2 + ALPHA2^2
TRACE_ALPHA = 2*ALPHA1
L2 = squarefree_part(NORM_ALPHA)
L1 = isqrt(NORM_ALPHA/L2)


# L = L1*L2*L3 is the product of the NUMBER_OF_PRIMES smallest primes that split in frakO_0.

L = completeL(NUMBER_OF_PRIMES, L1*L2)
L3 = L // (L1*L2)

p = 4*c*L - 1
Fp2.<i> = GF((p,2), modulus=[1,0,1])
assert i^2 == -1

A, xP, xQ = map(Fp2, (A, xP, xQ))

# starting curve and points

E0 = EllipticCurve([0,A,0,1,0])
E0.set_order((p+1)^2)
PE0 = E0.lift_x(xP)
QE0 = E0.lift_x(xQ)

iso = E0.isomorphism_to(E0.short_weierstrass_model())
E0 = iso.codomain()
PE0, QE0 = map(iso, (PE0, QE0))
E0.set_order(iso.domain().order())

assert PE0 in E0 and QE0 in E0
assert PE0.order() == L1*L2
assert QE0.order() == L1

_, _, _, a4, a6 = E0.a_invariants()
assert E0 == EllipticCurve([a4,a6])

xP = PE0.xy()[0]
xQ = QE0.xy()[0]
assert E0.lift_x(xP) in (+PE0, -PE0)
assert E0.lift_x(xQ) in (+QE0, -QE0)

if args.verbose:
    print(E0, file=sys.stderr)
    print(PE0, file=sys.stderr)
    print(QE0, file=sys.stderr)


######## constructing the relation lattice

f = ALPHA2      # conductor
hh = (f+1)//2   # class number

# characteristic polynomial of α
ZZx.<x> = ZZ[]
charpol = x^2 - TRACE_ALPHA*x + NORM_ALPHA

ells = sorted(L1.prime_factors() + L2.prime_factors() + L3.prime_factors())

# eigenvalues mod ℓ (arbitrarily sorted by magnitude)
if args.verbose:
    print("Computing eigenvalues", file=sys.stderr)

ev = [sorted(charpol.roots(GF(ell), multiplicities=False)) for ell in ells]

if args.verbose:
    for ell, lams in zip(ells, ev):
        print(ell, lams, file=sys.stderr)

# the generator of the "positive direction" ideal splitting each ℓ.
# "positive direction" defined by the arbitrary sorting above
if args.verbose:
    print("Computing ideal generators", file=sys.stderr)

O.<I> = GaussianIntegers()
gens = [O.ideal(λ.parent().order(), (ALPHA1 + I*ALPHA2) - λ.lift()).gens_reduced()[0]
        for λ, _ in ev]

# The class group, mapped to the subgroup of the algebraic torus
# T₂(GF(f))² ⊂ T₂(GF(f)) ⊂ GF(f²)^×
if args.verbose:
    print("Computing class group structure", file=sys.stderr)

F.<ibar> = GF((f, 2), modulus=[1,0,1])
clsgrp = [l.polynomial().change_ring(GF(f))(x=ibar)^(2-2*f)
          for l in gens]

assert all(g^hh == 1
           for g in clsgrp)

# Find a class generating clsgrp
facts = hh.factor()
geni, gen = next((i, g)
                 for i, g in enumerate(clsgrp)
                 if all(g^(hh/r[0]) != 1
                        for r in facts))

import os
dlogs_filename = f'dlogs--{args.parameter}.txt'
if args.dlogs or not os.path.exists(dlogs_filename):
    import time
    # Compute the dlogs of all classes to base gen
    dur = {
        '128': '2 minutes',
        '256': 'forever',
        }
    if args.verbose:
        print(f'''(Re)computing the discrete logarithms.
This is going to take about {dur[args.parameter]}... ''', end='', file=sys.stderr)

    t = -time.monotonic()
    ord = [hh, pari.factor(hh)]
    dlogs = [ZZ(pari.fflog(l, gen, ord))
             for l in clsgrp]
    t += time.monotonic()

    if args.verbose:
        print(f'done in {round(t)} seconds!', file=sys.stderr)

    with open(dlogs_filename, 'wb') as df:
        df.write(f'h: {hh :#x}\n'.encode())
        for ell, (e, _), d in zip(ells, ev, dlogs):
            df.write(f'({ell}, θ-{e}): {d:#x}\n'.encode())
else:
    with open(dlogs_filename, 'rb') as df:
        import re
        line = df.readline().decode()
        assert ZZ(re.match('h: (0x[0-9a-f]+)\n', line).group(1)) == hh
        dlogs = []
        for i, (line, ell, (e, _)) in enumerate(zip(df, ells, ev)):
            m = re.match('\(([0-9]+), θ-([0-9]+)\): (0x[0-9a-f]+)\n', line.decode())
            assert int(m.group(1)) == ell, f'{m.group(1)} ≠ {ell}'
            assert int(m.group(2)) == e,  f'{m.group(2)} ≠ {e}'
            if i == geni:
                assert ZZ(m.group(3)) == 1
            dlogs.append(ZZ(m.group(3)))

assert all(gen^dlog == g and 0 <= dlog < hh
           for dlog, g in zip(dlogs, clsgrp))

# Construct and BKZ-reduce the relation matrix
if args.verbose:
    print("Reducing relation matrix", file=sys.stderr)

import fpylll

block_size = { '128' : 30, '256' : 40 }[args.parameter]   # adjust this if needed

M = fpylll.IntegerMatrix(len(dlogs), len(dlogs))
for i in range(len(dlogs)):
    M[i,geni] = dlogs[i]
    M[i,i] -= 1
M[geni,geni] = hh
Mbkz = fpylll.BKZ.reduction(M, fpylll.BKZ.Param(block_size))

assert all(prod(g^e for g, e in zip(clsgrp, row)) == 1
           for row in Mbkz)

################################################################

sz = lambda t: t.bit_length()//64*64 + 64

print('\n#pragma once')
print()
print('#include "bigint.hpp"')
print()

print(f'constexpr int_t<{sz(ALPHA1):3}> ALPHA1 = {ALPHA1:+#x}_cppi;')
print(f'constexpr int_t<{sz(ALPHA2):3}> ALPHA2 = {ALPHA2:+#x}_cppi;')
print(f'constexpr int_t<{sz(NORM_ALPHA):3}> NORM_ALPHA = {NORM_ALPHA:+#x}_cppi;')       #XXX derived
print(f'constexpr int_t<{sz(TRACE_ALPHA):3}> TRACE_ALPHA = {TRACE_ALPHA:+#x}_cppi;')    #XXX derived
print()

print(f'constexpr uint_t<{sz(L1):3}> L1 = {L1:#x}_cppui;')                              #XXX derived
print(f'constexpr uint_t<{sz(L2):3}> L2 = {L2:#x}_cppui;')                              #XXX derived
print(f'constexpr uint_t<{sz(L3):3}> L3 = {L3:#x}_cppui;')                              #XXX derived
print(f'constexpr uint_t<{sz(L):3}> L = {L:#x}_cppui;')                                 #XXX derived
print()

print(f'constexpr uint_t<{sz(2*p):3}> p = {p:#x}_cppui;')                               #XXX derived
print()

print(f'#include "fp.hpp"')
print(f'#include "fq.hpp"')
print()

print(f'fq const a4 = fq::make({a4[0].lift():#x}_cppui, {a4[1].lift():#x}_cppui);')
print(f'fq const a6 = fq::make({a6[0].lift():#x}_cppui, {a6[1].lift():#x}_cppui);')
print(f'fq const xP = fq::make({xP[0].lift():#x}_cppui, {xP[1].lift():#x}_cppui);')
print(f'fq const xQ = fq::make({xQ[0].lift():#x}_cppui, {xQ[1].lift():#x}_cppui);')
print()

print(f'constexpr uint_t<{sz(hh):3}> hh = {hh:#x}_cppui;')
print(f'constexpr matrix<int16_t, {len(dlogs)}, {len(dlogs)}> relations = {{{{')
for i in range(len(dlogs)):
    print(f'    {{{{ { ", ".join(str(e) for e in Mbkz[i] )} }}}},')
print('}};')

print(f'constexpr int geni = { geni };')
