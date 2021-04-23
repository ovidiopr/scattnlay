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
        self.evaluated_data.sort(key=lambda x: x.testname)  # , reverse=True)
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

    def compose_result_string(self, mpf_x, mpf_m, n, mpf_value, output_dps):
        return ('{'+
                 mp.nstr(mpf_x, output_dps * 2) + ',{' +
                 mp.nstr(mpf_m.real, output_dps * 2) + ',' +
                 mp.nstr(mpf_m.imag, output_dps * 2) + '},' +
                 str(n) + ',{' +
                 mp.nstr(mpf_value.real, output_dps) + ',' +
                 mp.nstr(mpf_value.imag, output_dps) + '},' +
                 mp.nstr(mp.fabs(mpf_value.real * 10 ** -output_dps), 2) + ',' +
                 mp.nstr(mp.fabs(mpf_value.imag * 10 ** -output_dps), 2) +
                 '},')

    def get_test_data_nlist(self, z_record, output_dps, n, func):
        isNeedMoreDPS = False
        x = str(z_record[0])
        mr = str(z_record[1][0])
        mi = str(z_record[1][1])
        z_str = ''
        try:
            mpf_x = mp.mpf(x)
            mpf_m = mp.mpc(mr, mi)
            z = mpf_x*mpf_m
            if self.is_only_x: z = mp.mpf(x)
            if self.is_xm:
                mpf_value = func(n, mpf_x, mpf_m)
            else:
                mpf_value = func(n, z)
            z_str = self.compose_result_string(mpf_x, mpf_m, n, mpf_value, output_dps)
            if mp.nstr(mpf_value.real, output_dps) == '0.0' \
                    or mp.nstr(mpf_value.imag, output_dps) == '0.0':
                isNeedMoreDPS = True
        except:
            isNeedMoreDPS = True
        return z_str, isNeedMoreDPS

    def get_test_data(self, Du_test, output_dps, max_num_elements_of_n_list, func, funcname):
        output_list = ['// x, complex(m), n, complex(f(n,z)), abs_err_real, abs_err_imag',
                       'std::vector< std::tuple< nmie::FloatType, std::complex<nmie::FloatType>, int, std::complex<nmie::FloatType>, nmie::FloatType, nmie::FloatType > >',
                       str(funcname) + '_test_' + str(output_dps) + 'digits', '= {']
        for z_record in Du_test:
            x = str(z_record[0])
            mr = str(z_record[1][0])
            mi = str(z_record[1][1])
            mp.mp.dps = 20
            z = mp.mpf(x) * mp.mpc(mr, mi)
            n_list = self.get_n_list(z, max_num_elements_of_n_list)
            if z_record[4] == 'Yang': n_list = [0, 1, 30, 50, 60, 70, 75, 80, 85, 90, 99, 116, 130]
            print(z, n_list)
            failed_evaluations = 0
            for n in n_list:
                mp.mp.dps = output_dps
                old_z_string, isNeedMoreDPS = self.get_test_data_nlist(z_record, output_dps, n, func, )
                mp.mp.dps = int(output_dps*1.41)
                new_z_string, isNeedMoreDPS = self.get_test_data_nlist(z_record, output_dps, n, func)
                while old_z_string != new_z_string \
                        or isNeedMoreDPS:
                    new_dps = int(mp.mp.dps * 1.41)
                    if new_dps > 300: break
                    mp.mp.dps = new_dps
                    print("New dps = ", mp.mp.dps, 'n =', n, ' (max ', n_list[-1], ') for z =', z, '     ', end='')
                    old_z_string = new_z_string
                    new_z_string, isNeedMoreDPS = self.get_test_data_nlist(z_record, output_dps, n, func)

                if new_z_string != '':
                    output_list.append(new_z_string)
                else:
                    failed_evaluations += 1
                #     break
            result_str = "All done!"
            if failed_evaluations > 0: result_str = " FAILED!"
            print("\n", result_str, "Failed evaluations ", failed_evaluations, ' of ', len(n_list))
        output_list.append('};')
        return output_list

    def run_test(self, func, funcname, is_only_x=False, is_xm=False):
        self.is_only_x = is_only_x
        self.is_xm = is_xm
        self.remove_argument_duplicates()
        out_list_result = self.get_test_data(self.complex_arguments, self.output_dps,
                                             self.max_num_elements_of_nlist,
                                             func, funcname)
        testname = str(funcname) + '_test_' + str(self.output_dps) + 'digits'
        self.remove(testname)
        self.add_record(out_list_result)

    def remove_argument_duplicates(self):
        print("Arguments in input: ", len(self.complex_arguments))
        mp.mp.dps = 20
        self.complex_arguments.sort()
        filtered_list = []
        filtered_list.append(self.complex_arguments[0])
        for i in range(1, len(self.complex_arguments)):
            # if x and m are the same: continue
            if (filtered_list[-1][0] == self.complex_arguments[i][0] and
                    filtered_list[-1][1] == self.complex_arguments[i][1]):
                continue
            # argument list is sorted, so when only x is needed
            # keep the record with the largest m
            if (self.is_only_x
                    and filtered_list[-1][0] == self.complex_arguments[i][0]):
                # continue
                del filtered_list[-1]
            filtered_list.append(self.complex_arguments[i])
        self.complex_arguments = filtered_list
        # print(self.complex_arguments)
        print("Arguments after filtering: ", len(self.complex_arguments))
        # exit(0)


def main():
    sf_evals = UpdateSpecialFunctionsEvaluations(filename='test_spec_functions_data.hpp',
                                                 complex_arguments=mia.complex_arguments,
                                                 output_dps=30, max_num_elements_of_nlist=51)
                                                 # output_dps=7, max_num_elements_of_nlist=51)
                                                 # output_dps=5, max_num_elements_of_nlist=3)
    # sf_evals.run_test(mrb.D1, 'D1')
    # sf_evals.run_test(mrb.D2, 'D2')
    # sf_evals.run_test(mrb.D3, 'D3')
    # sf_evals.run_test(mrb.psi, 'psi', is_only_x=True)
    # sf_evals.run_test(mrb.xi, 'xi', is_only_x=True)
    # # In literature Zeta or Ksi denote the Riccati-Bessel function of third kind.
    # sf_evals.run_test(mrb.ksi, 'zeta', is_only_x=True)

    # sf_evals.run_test(mrb.an, 'an', is_xm=True)
    # sf_evals.run_test(mrb.bn, 'bn', is_xm=True)

    # sf_evals.run_test(mrb.psi, 'psi')
    # sf_evals.run_test(mrb.psi_div_ksi, 'psi_div_ksi')
    # sf_evals.run_test(mrb.psi_mul_ksi, 'psi_mul_zeta', is_only_x=True)
    # sf_evals.run_test(mrb.psi_div_xi, 'psi_div_xi')
    with open(sf_evals.filename, 'w') as out_file:
        out_file.write(sf_evals.get_file_content())

    for record in mia.complex_arguments:
        mp.mp.dps = 20
        output_dps = 16
        x = mp.mpf(str(record[0]))
        mr = str(record[1][0])
        mi = str(record[1][1])
        m = mp.mpc(mr, mi)
        Qext_ref = record[2]
        Qsca_ref = record[3]
        test_case = record[4]
        nmax = int(x + 4.05*x**(1./3.) + 2)+2+28
        print(f"\n ===== test case: {test_case} =====", flush=True)
        print(f"x={x}, m={m}, N={nmax} \nQsca_ref = {Qsca_ref}    \tQext_ref = {Qext_ref}", flush=True)
        Qext_mp = mrb.Qext(x,m,nmax, output_dps)
        Qsca_mp = mrb.Qsca(x,m,nmax, output_dps)
        print(f"Qsca_mp  = {mp.nstr(Qsca_mp[-1],output_dps)}    \tQext_mp  = {mp.nstr(Qext_mp[-1],output_dps)}", flush=True)
        print(mp.nstr(Qsca_mp,output_dps))
        print(mp.nstr(Qext_mp,output_dps))

        # n=1
        # print(f'n={n}, x={x}, m={m}\nbn[{n}]={mp.nstr(mrb.bn(n,x,m), output_dps)}')

main()
