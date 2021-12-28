#ifndef SCATTNLAY_TESTS_TEST_CASES_HPP_
#define SCATTNLAY_TESTS_TEST_CASES_HPP_
//
//
#include <vector>
#include <tuple>
#include <string>
#include <complex>


std::vector< std::tuple< double, std::complex<double>, std::string > >
parameters_bulk_sphere
{
// x, {Re(m), Im(m)}, test_name
{0.099, {0.75,0}, "Hong Du testcase:a"},
{0.101, {0.75,0}, "Hong Du testcase:b"},
{10,    {0.75,0}, "Hong Du testcase:c"},
{100,   {1.33,1e-5}, "Hong Du testcase:e"},
{0.055, {1.5, 1},    "Hong Du testcase:g"},
{0.056, {1.5, 1}, "Hong Du testcase:h"},
{100,   {1.5, 1}, "Hong Du testcase:i"},
{1,     {10,  10},   "Hong Du testcase:k"},
{1000,  {0.75,0}, "Hong Du testcase:d"},
//{10000, {1.33,1e-5}, "Hong Du testcase:f"}, // passes but takes too long
//{100,   {10,  10,},  "Hong Du testcase:l"}, // fails in any precision, TODO fixme
//{10000, {1.5, 1},    "Hong Du testcase:j"},
//{10000, {10,  10},   "Hong Du testcase:m"},
#ifdef MULTI_PRECISION
//{10000, {1.33,1e-5}, "Hong Du testcase:f"}, // passes but takes too long
//{100,   {10,  10,},  "Hong Du testcase:l"}, // fails in any precision
//{10000, {1.5, 1},    "Hong Du testcase:j"},
//{10000, {10,  10},   "Hong Du testcase:m"},
#endif
};

#endif //SCATTNLAY_TESTS_TEST_CASES_HPP_
