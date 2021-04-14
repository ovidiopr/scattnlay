#!/usr/bin/env python3
import mpmath as mp
import mpmath_riccati_bessel as mrb
import mpmath_input_arguments as mia

class test_data:
    def __init__(self, list_to_parse):
        self.evaluated_data = []
        self.comment = list_to_parse[0]
        if self.comment[:2] != '//': raise ValueError('Not a comment')
        self.ending = list_to_parse[-1]
        if self.ending != '};': raise ValueError('For C++ we expect closing };')
        self.list_to_parse = list_to_parse


class update_special_functions_evaluations:
    def __init__(self, filename='default_out.h', complex_arguments = []):
        self.evaluated_data = []
        self.test_setup = []
        self.filename = filename
        self.complex_arguments = complex_arguments

    def read_evaluated_data(self):
        with open(self.filename, 'r') as in_file:
            content = in_file.readlines()
        content = [x.strip() for x in content]
        while '' in content:
            record_end_index = content.index('')
            new_record = content[:record_end_index]
            self.evaluated_data.append(test_data(new_record))
            content = content[record_end_index+1:]
        self.evaluated_data.append(test_data(new_record))


def get_n_list(z, max_number_of_elements = 10):
    nmax = mrb.LeRu_cutoff(z)
    factor = nmax**(1/(max_number_of_elements-2))
    n_list = [int(factor**i) for i in range(max_number_of_elements-1) ]
    n_list.append(0)
    n_set = set(n_list)
    return sorted(n_set)


def get_test_data_nlist(z_record, output_dps, n):
    x = str(z_record[0])
    mr = str(z_record[1][0])
    mi = str(z_record[1][1])
    z_str = ''
    try:
        z = mp.mpf(x)*mp.mpc(mr, mi)
        D1nz = mrb.D1(n,z)
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
    output_str = ('// complex(z), n, complex(D1(n,z)), abs_err_real, abs_err_imag\n'+
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
    sf_evals = update_special_functions_evaluations(filename='test_spec_functions_data.hpp',
                                                    complex_arguments = mia.complex_arguments)
    sf_evals.read_evaluated_data()
    print(sf_evals.evaluated_data[1].list_to_parse)
    output_dps = 16
    max_num_elements_of_nlist = 51

    # out_filename = 'test_spec_functions_data.hpp'
    # output_str = get_test_data(mia.complex_arguments, output_dps, max_num_elements_of_nlist)
    # with open(out_filename, 'w') as out_file:
    #     out_file.write(output_str)

main()
