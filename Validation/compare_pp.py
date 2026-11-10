#!/usr/bin/env python3
"""Compares a band structure (.agr) or density of states (.gp) file
with a reference file, for the Validation/verify_pp script.

  python3 compare_pp.py new_file reference_file [nband]

nband (optional): compare only the lowest nband bands at each k-point
(the highest bands of an interpolation are the least accurate).

Band files (xmgrace .agr): the data sets are separated by blank lines.
The band sets are the sets with the most common number of points (the
2-point sets for the zero line and the vertical lines are ignored).
The energies are sorted at each k-point before comparing, so a different
band ordering (state matching across k, which depends on the LAPACK
library and on the number of threads) does not count as a difference.
A third column (unfolding weight or atomic orbital character) is compared
after pairing it with the sorted energies; degenerate states may give a
large difference there, so it is printed only as information.

Density of states files (gnuplot .gp): all columns after the energy are
compared point by point.

Written 10 November 2026. claude
"""

import sys


def read_sets(filename):
    """Returns the header lines and the data sets (lists of rows of floats)."""
    sets, current, header = [], [], []
    with open(filename) as f:
        for line in f:
            fields = line.split()
            row = None
            if fields and line.lstrip()[0] not in '@#&':
                try:
                    row = [float(x) for x in fields]
                except ValueError:
                    row = None
            if row is None:
                if current:
                    sets.append(current)
                    current = []
                if not sets:
                    header.append(line)
            else:
                current.append(row)
    if current:
        sets.append(current)
    return header, sets


def band_matrix(sets):
    """Returns x(k) and e(band, k) [and w(band, k) if present]."""
    lengths = [len(s) for s in sets if len(s) > 2]
    if not lengths:
        return None
    nk = max(set(lengths), key=lengths.count)
    bands = [s for s in sets if len(s) == nk]
    x = [row[0] for row in bands[0]]
    e = [[row[1] for row in b] for b in bands]
    w = None
    if all(len(row) > 2 for b in bands for row in b):
        w = [[row[2] for row in b] for b in bands]
    return x, e, w


def compare_bands(new, ref, nband=None):
    a, b = band_matrix(new), band_matrix(ref)
    if a is None or b is None:
        print('   no band data found')
        return 2
    xa, ea, wa = a
    xb, eb, wb = b
    print('   bands: {} sets of {} k-points (reference {} of {})'.format(
          len(ea), len(xa), len(eb), len(xb)))
    if len(ea) != len(eb) or len(xa) != len(xb):
        print('   DIFFERENT number of bands or k-points')
        return 2

    dx = max(abs(u - v) for u, v in zip(xa, xb))
    raw = max(abs(u - v) for ba, bb in zip(ea, eb) for u, v in zip(ba, bb))
    nb = len(ea)
    if nband:
        nb = min(nband, nb)
        print('   only the lowest {} bands are compared'.format(nb))
    srt, wdif, kmax, nmax = 0.0, 0.0, 0, 0
    for k in range(len(xa)):
        pa = sorted((ea[n][k], wa[n][k] if wa else 0.0) for n in range(len(ea)))[:nb]
        pb = sorted((eb[n][k], wb[n][k] if wb else 0.0) for n in range(len(eb)))[:nb]
        for n, (u, v) in enumerate(zip(pa, pb)):
            if abs(u[0] - v[0]) > srt:
                srt, kmax, nmax = abs(u[0] - v[0]), k, n
        if wa and wb:
            wdif = max(wdif, max(abs(u[1] - v[1]) for u, v in zip(pa, pb)))

    print('   maximum difference in data (energies sorted at each k-point) {:12.6f}'.format(srt))
    print('      band {} of {} (sorted) at k-point {} (x = {:.5f})'.format(nmax + 1, len(ea), kmax + 1, xa[kmax]))
    print('   maximum difference in the k-path coordinate                  {:12.6f}'.format(dx))
    print('   maximum difference with the bands as ordered in the files    {:12.6f}'.format(raw))
    if wa and wb:
        print('   maximum difference in the third column (paired by energy)    {:12.6f}'.format(wdif))
    return 0


def compare_dos(new, ref):
    """The energy range may differ by a grid point; rows are matched by energy."""
    a = {round(row[0], 5): row for s in new for row in s}
    b = {round(row[0], 5): row for s in ref for row in s}
    common = sorted(set(a) & set(b))
    print('   dos: {} energies (reference {}), {} in common'.format(len(a), len(b), len(common)))
    if not common:
        print('   no common energies')
        return 2
    for e in sorted(set(a) ^ set(b)):
        row = a.get(e, b.get(e))
        print('      only in {}: {}'.format('new' if e in a else 'reference',
              '  '.join('{:.6f}'.format(v) for v in row)))
    ncol = min(len(a[common[0]]), len(b[common[0]]))
    for j in range(1, ncol):
        d = max(abs(a[e][j] - b[e][j]) for e in common)
        big = max(abs(b[e][j]) for e in common)
        print('   maximum difference in column {}  {:12.6f}   (largest value {:10.4f})'.format(j + 1, d, big))
    return 0


def main():
    if len(sys.argv) not in (3, 4):
        print(__doc__)
        return 1
    nband = int(sys.argv[3]) if len(sys.argv) == 4 else None
    hnew, snew = read_sets(sys.argv[1])
    href, sref = read_sets(sys.argv[2])
    print('   comparing {} with {}'.format(sys.argv[1], sys.argv[2]))
    if any('plot' in line for line in href):
        return compare_dos(snew, sref)
    return compare_bands(snew, sref, nband)


if __name__ == '__main__':
    sys.exit(main())
