#!/usr/bin/env python3
import mpmath as mp
import mpmath_riccati_bessel as mrb
import mpmath_input_arguments as mia
import os.path


class TestData:
    def __init__(self, list_to_parse, filetype):
        self.filetype = filetype
        if self.filetype == 'c++':
            self.cpp_parse(list_to_parse)
        else:
            raise NotImplementedError("Only C++ files *.hpp parsing was implemented")

    def cpp_parse(self, list_to_parse):
        self.comment = list_to_parse[0]
        if self.comment[:2] != '//': raise ValueError('Not a comment')
        self.typeline = list_to_parse[1]
        if 'std::vector' not in self.typeline: raise ValueError('Unexpected C++ container')
        self.testname = list_to_parse[2]
        self.opening = list_to_parse[3]
        if self.opening != '= {': raise ValueError('For C++ we expect opeing with = {');
        self.ending = list_to_parse[-1]
        if self.ending != '};': raise ValueError('For C++ we expect closing };')
        self.evaluated_data = list_to_parse[4:-1]

    def get_string(self):
        out_sting = self.comment + '\n' + self.typeline + '\n' + self.testname + '\n' + self.opening + '\n'
        for result in self.evaluated_data:
            out_sting += result + '\n'
        out_sting += self.ending + '\n'
        return out_sting


class UpdateSpecialFunctionsEvaluations:
    def __init__(self, filename='default_out.hpp', complex_arguments=[],
                 output_dps=16, max_num_elements_of_nlist=51):
        self.evaluated_data = []
        self.test_setup = []
        self.filename = filename
        self.read_evaluated_data()
        self.complex_arguments = complex_arguments
        self.output_dps = output_dps
        self.max_num_elements_of_nlist = max_num_elements_of_nlist

    def read_evaluated_data(self):
        self.filetype = 'undefined'
        if self.filename.endswith('.hpp'):
            self.filetype = 'c++'
        if self.filename.endswith('.f90'):
            self.filetype = 'fortran'
        if not os.path.exists(self.filename):
            print("WARNING! Found no data file:", self.filename)
            return
        with open(self.filename, 'r') as in_file:
            content = in_file.readlines()
        content = [x.strip() for x in content]
        while '' in content:
            record_end_index = content.index('')
            new_record = content[:record_end_index]
            content = content[record_end_index + 1:]
            self.add_record(new_record)
        self.add_record(content)

    def add_record(self, new_record):
        if len(new_record) == 0: return
        if len(new_record) < 6: raise ValueError('Not enough lines in record:', new_record)
        self.evaluated_data.append(TestData(new_record, self.filetype))


    def get_file_content(self):
        self.evaluated_data.sort(key=lambda x: x.testname)#, reverse=True)
        out_string = ''
        for record in self.evaluated_data:
            out_string += record.get_string() + '\n'
        return out_string[:-1]

    def remove(self, testname):
        for i, result in enumerate(self.evaluated_data):
            if result.testname == testname:
                del self.evaluated_data[i]

    def get_n_list(self, z, max_number_of_elements=10):
        nmax = mrb.LeRu_cutoff(z)
        factor = nmax ** (1 / (max_number_of_elements - 2))
        n_list = [int(factor ** i) for i in range(max_number_of_elements - 1)]
        n_list.append(0)
        n_set = set(n_list)
        return sorted(n_set)

    def get_test_data_nlist(self, z_record, output_dps, n, func):
        x = str(z_record[0])
        mr = str(z_record[1][0])
        mi = str(z_record[1][1])
        z_str = ''
        try:
            z = mp.mpf(x) * mp.mpc(mr, mi)
            D1nz = func(n, z)
            z_str = ('{{' +
                     mp.nstr(z.real, output_dps * 2) + ',' +
                     mp.nstr(z.imag, output_dps * 2) + '},' +
                     str(n) + ',{' +
                     mp.nstr(D1nz.real, output_dps) + ',' +
                     mp.nstr(D1nz.imag, output_dps) + '},' +
                     mp.nstr(mp.fabs(D1nz.real * 10 ** -output_dps), 2) + ',' +
                     mp.nstr(mp.fabs(D1nz.imag * 10 ** -output_dps), 2) +
                     '},')
        except:
            pass
        return z_str

    def get_test_data(self, Du_test, output_dps, max_num_elements_of_n_list, func, funcname):
        output_list = ['// complex(z), n, complex(f(n,z)), abs_err_real, abs_err_imag',
        'std::vector< std::tuple< std::complex<double>, int, std::complex<double>, double, double > >',
        str(funcname)+'_test_' + str(output_dps) + 'digits','= {']
        for z_record in Du_test:
            x = str(z_record[0])
            mr = str(z_record[1][0])
            mi = str(z_record[1][1])
            mp.mp.dps = 20
            z = mp.mpf(x) * mp.mpc(mr, mi)
            n_list = self.get_n_list(z, max_num_elements_of_n_list)
            if z_record[4] == 'Yang': n_list = [0,1,30,50,60,70,75,80,85,90,99,116,130]
            print(z, n_list)
            failed_evaluations = 0
            for n in n_list:
                mp.mp.dps = 20
                old_z_string = self.get_test_data_nlist(z_record, output_dps, n, func)
                mp.mp.dps = 37
                new_z_string = self.get_test_data_nlist(z_record, output_dps, n, func)
                while old_z_string != new_z_string:
                    new_dps = int(mp.mp.dps * 1.41)
                    if new_dps > 300: break
                    mp.mp.dps = new_dps
                    print("New dps = ", mp.mp.dps, 'n =', n, ' (max ',n_list[-1],') for z =', z, '     ', end='')
                    old_z_string = new_z_string
                    new_z_string = self.get_test_data_nlist(z_record, output_dps, n, func)

                if new_z_string != '':
                    output_list.append(new_z_string)
                else:
                    failed_evaluations += 1
                #     break
            print("\nFailed evaluations ", failed_evaluations, ' of ', len(n_list))
        output_list.append('};')
        return output_list

    def run_test(self, func, funcname):
        out_list_result = self.get_test_data(mia.complex_arguments, self.output_dps,
                                             self.max_num_elements_of_nlist,
                                             func, funcname)
        testname = str(funcname)+'_test_' + str(self.output_dps) + 'digits'
        self.remove(testname)
        self.add_record(out_list_result)


def main():
    sf_evals = UpdateSpecialFunctionsEvaluations(filename='test_spec_functions_data.hpp',
                                                 complex_arguments=mia.complex_arguments,
                                                 output_dps=16, max_num_elements_of_nlist=51)
                                                 # output_dps=5, max_num_elements_of_nlist=3)
    # sf_evals.run_test(mrb.D1, 'D1')
    sf_evals.run_test(mrb.D2, 'D2')
    # sf_evals.run_test(mrb.D3, 'D3')
    # sf_evals.run_test(mrb.psi_div_ksi, 'psi_div_ksi')
    # sf_evals.run_test(mrb.psi_div_xi, 'psi_div_xi')
    with open(sf_evals.filename, 'w') as out_file:
        out_file.write(sf_evals.get_file_content())


main()
