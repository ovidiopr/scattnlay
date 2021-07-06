#!/usr/bin/env python3
import mpmath as mp


# APPLIED OPTICS / Vol. 53, No. 31 / 1 November 2014, eq(13)
def LeRu_cutoff(z):
    x = mp.fabs(z)
    return int(x + 11 * x**(1/3) + 1)


# Wu, Wang, Radio Science, Volume 26, Number 6, Pages 1393-1401, November-December 1991,
# after eq 13f.
# Riccati-Bessel z*j_n(z)
def psi(n,z):
    return mp.sqrt( (mp.pi * z)/2 ) * mp.autoprec(mp.besselj)(n+1/2,z)
# Riccati-Bessel -z*y_n(z)
def xi(n,z):
    return -mp.sqrt( (mp.pi * z)/2 ) * mp.autoprec(mp.bessely)(n+1/2,z)
# Riccati-Bessel psi - i* xi        z*j_n(z) + i* z*j_n(z)
def ksi(n,z):
    return psi(n,z) - 1.j * xi(n,z)

def psi_div_ksi(n,z):
    return psi(n,z)/ksi(n,z)

def psi_mul_ksi(n,z):
    return psi(n,z)*ksi(n,z)

def psi_div_xi(n,z):
    return psi(n,z)/xi(n,z)

# to compare r(n,z) with Wolfram Alpha
# n=49, z=1.3-2.1i,  SphericalBesselJ[n-1,z]/SphericalBesselJ[n,z]
def r(n,z):
    if n > 0:
        return psi(n-1,z)/psi(n,z)
    return mp.cos(z)/mp.sin(z)


def D(n, z, f):
    return f(n-1,z)/f(n,z) - n/z

def D1(n,z):
    if n == 0: return mp.cos(z)/mp.sin(z)
    return D(n, z, psi)

# Wolfram Alpha example D2(10, 10-10j): SphericalBesselY[9, 10-10i]/SphericalBesselY[10,10-10i]-10/(10-10i)
def D2(n,z):
    if n == 0: return -mp.sin(z)/mp.cos(z)
    return D(n, z, xi)

def D3(n,z):
    if n == 0: return 1j
    return D(n, z, ksi)


# bulk sphere
# Ovidio
def an(n, x, m):
    # print(f'D1 = {D1(n, m*x)}\n\
    # psi_n = {psi(n,x)}\n\
    # psi_nm1 = {psi(n-1,x)}\n\
    # ksi_n = {ksi(n,x)}\n\
    # ksi_nm1 = {ksi(n-1,x)}\n\
    # ')
    return (
            ( ( D1(n, m*x)/m + n/x )*psi(n,x) - psi(n-1,x) ) /
            ( ( D1(n, m*x)/m + n/x )*ksi(n,x) - ksi(n-1,x) )
    )
def bn(n, x, m):
    # print(f'D1 = {D1(n, m*x)}\n\
    # psi_n = {psi(n,x)}\n\
    # psi_nm1 = {psi(n-1,x)}\n\
    # ksi_n = {ksi(n,x)}\n\
    # ksi_nm1 = {ksi(n-1,x)}\n\
    # ')
    return (
            ( ( D1(n, m*x)*m + n/x ) * psi(n,x) - psi(n-1,x) ) /
            ( ( D1(n, m*x)*m + n/x ) * ksi(n,x) - ksi(n-1,x) )
    )

# # Du
# def an(n, x, m):
#     mx = m*x
#     return (
#         ((r(n,mx)/m + n*(1-1/m**2)/x)*psi(n,x)-psi(n-1,x))/
#         ((r(n,mx)/m + n*(1-1/m**2)/x)*ksi(n,x)-ksi(n-1,x))
#     )
# def bn(n,x,m):
#     mx = m*x
#     return (
#         (r(n,mx)*m*psi(n,x)-psi(n-1,x))/
#         (r(n,mx)*m*ksi(n,x)-ksi(n-1,x))
#     )

def Qext_diff(n,x,m):
    return ((2*n +1)*2/x**2) * (an(n,x,m) + bn(n,x,m)).real

def Qsca_diff(n,x,m):
    return ((2*n +1)*2/x**2) * (abs(an(n,x,m))**2 + abs(bn(n,x,m))**2)

def Qany(x, m, nmax, output_dps, Qdiff):
    Qany = 0
    Qlist = []
    Qprev=""
    convergence_max_length = 35
    i = 0
    for n in range(1, nmax+1):
        diff = Qdiff(n,x,m)
        Qany += diff
        Qnext = mp.nstr(Qany, output_dps)
        Qlist.append(Qany)
        i += 1
        if Qprev !=Qnext: i = 0
        Qprev = Qnext
        print(n, ' ', end='', flush=True)
        if i >= convergence_max_length: break
    print('')
    Qlist.append(Qany)
    return Qlist

def Qsca(x, m, nmax, output_dps):
    return Qany(x, m, nmax, output_dps, Qsca_diff)

def Qext(x, m, nmax, output_dps):
    return Qany(x, m, nmax, output_dps, Qext_diff)
