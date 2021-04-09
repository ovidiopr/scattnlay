#!/usr/bin/env python3
import mpmath as mp
import numpy as np
Du_test = [
# // x, [Re(m), Im(m)], Qext, Qsca, test_name
[0.099, [0.75,0], 7.417859e-06, 7.417859e-06, 'a'],
[0.101, [0.75,0], 8.033538e-06, 8.033538e-06, 'b'],
[10,    [0.75,0],     2.232265, 2.232265, 'c'],
[1000,  [0.75,0],     1.997908, 1.997908, 'd'],
[100,   [1.33,-1e-5], 2.101321, 2.096594, 'e'],
[10000, [1.33,-1e-5], 2.004089, 1.723857, 'f'],
[0.055, [1.5, -1],    0.101491, 1.131687e-05, 'g'],
[0.056, [1.5, -1],   0.1033467, 1.216311e-05, 'h'],
[100,   [1.5, -1],    2.097502, 1.283697, 'i'],
[10000, [1.5, -1],    2.004368, 1.236574, 'j'],
[1,     [10,  -10],   2.532993, 2.049405, 'k'],
[100,   [10,  -10,],  2.071124, 1.836785, 'l'],
[10000, [10,  -10],   2.005914, 1.795393, 'm'],
[80, [1.05,  1],   0, 0, 'Yang'],
[1, [mp.pi,  1],   0, 0, 'pi'],
[1, [mp.pi,  -1],   0, 0, 'pi'],
[1, [mp.pi,  mp.pi],   0, 0, 'pi'],
[1, [2*mp.pi,  -1],   0, 0, 'pi'],
[1, [2*mp.pi,  mp.pi],   0, 0, 'pi'],
[1, [2*mp.pi,  1],   0, 0, 'pi'],
]
# // Dtest refractive index is m={1.05,1}, the size parameter is x = 80
n_list = [0,1,30,50,60,70,75,80,85,90,99,116,130];

def get_z_values(du_list):
    zlist = []
    for record in du_list:
        x = mp.mpf(str(record[0]))
        m = mp.mpc(str(record[1][0]), str(record[1][1]))
        z = x*m
        zlist.append(z)
    return zlist


# APPLIED OPTICS / Vol. 53, No. 31 / 1 November 2014, eq(13)
def LeRu_cutoff(z):
    x = mp.fabs(z)
    return int(x + 11 * x**(1/3) + 1)


def get_n_list(z, max_number_of_elements = 10):
    nmax = LeRu_cutoff(z)
    factor = nmax**(1/(max_number_of_elements-2))
    n_list = [int(factor**i) for i in range(max_number_of_elements-1) ]
    n_list.append(0)
    n_set = set(n_list)
    return sorted(n_set)

# Riccati-Bessel z*j_n(z)
def psi(n,z):
    return mp.sqrt( (mp.pi * z)/2 ) * mp.autoprec(mp.besselj)(n+1/2,z)


# to compare r(n,z) with Wolfram Alpha
# n=49, z=1.3-2.1i,  SphericalBesselJ[n-1,z]/SphericalBesselJ[n,z]
def r(n,z):
    if n > 0:
        return psi(n-1,z)/psi(n,z)
    return mp.cos(z)/mp.sin(z)


def D1(n,z):
    return r(n,z) - n/z


def get_test_data_nlist(z_record, output_dps, n):
    x = str(z_record[0])
    mr = str(z_record[1][0])
    mi = str(z_record[1][1])
    z_str = ''
    try:
        z = mp.mpf(x)*mp.mpc(mr, mi)
        D1nz = D1(n,z)
        z_str=('{{'+
                 mp.nstr(z.real, output_dps*2)+','+
                 mp.nstr(z.imag, output_dps*2)+'},'+
                 str(n)+',{'+
                 mp.nstr(D1nz.real, output_dps)+','+
                 mp.nstr(D1nz.imag, output_dps)+'},'+
                 mp.nstr(mp.fabs(D1nz.real* 10**-output_dps),2)+',' +
                 mp.nstr(mp.fabs(D1nz.imag* 10**-output_dps),2)+
                 '},\n')
    except:
        pass
    return z_str

def get_test_data(Du_test, output_dps, max_num_elements_of_n_list):
    output_str = ('// complex(z), n, complex(D1_n(z)), abs_err_real, abs_err_imag\n'+
    'std::vector< std::tuple< std::complex<double>, int, std::complex<double>, double, double > >'+
    'D1_test_'+str(output_dps)+'digits = {\n')
    for z_record in Du_test:
        x = str(z_record[0])
        mr = str(z_record[1][0])
        mi = str(z_record[1][1])
        mp.mp.dps = 20
        z = mp.mpf(x)*mp.mpc(mr, mi)
        n_list = get_n_list(z, max_num_elements_of_n_list)
        print(z, n_list)
        for n in n_list:
            mp.mp.dps = 20
            old_z_string = get_test_data_nlist(z_record, output_dps, n)
            mp.mp.dps = 37
            new_z_string = get_test_data_nlist(z_record, output_dps, n)
            while old_z_string != new_z_string:
                new_dps = int(mp.mp.dps * 1.41)
                if new_dps > 300: break
                mp.mp.dps = new_dps
                print("New dps = ", mp.mp.dps)
                old_z_string = new_z_string
                new_z_string = get_test_data_nlist(z_record, output_dps, n)

            output_str += new_z_string
    output_str += '};\n'
    return output_str


def main ():
    output_dps = 10
    max_num_elements_of_nlist = 15

    out_filename = 'test_spec_functions_data.h'
    output_str = get_test_data(Du_test, output_dps, max_num_elements_of_nlist)
    with open(out_filename, 'w') as out_file:
        out_file.write(output_str)
main()
