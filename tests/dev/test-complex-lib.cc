#include <complex>
#include <string>
#include <stdio.h>
#include "ucomplex.h"
void output(std::string operation, complex c1, std::complex<double> c2) {
    double diff_r = c1.r - c2.real();
    double diff_i = c1.i - c2.imag();
    printf("%s: \tdiff=(%g,%g), \tc1 = (%g,%g),  \tc2 = (%g,%g)\n",
	   operation.c_str(),
	   diff_r, diff_i,
	   c1.r, c1.i,
	   c2.real(), c2.imag());
    // Epsilon for double is 1e-16, check it 1.0 + 1e-16 = 1.0
    if (std::abs(diff_r) > 1e-16||std::abs(diff_i) > 1e-16) 
	//if (diff_r > 1e-15||diff_i > 1e-15) 
	printf("\n********\n\t\tWARNING!! Non-zero diff!!!\n********\n");
}
/********************************************************************/
/********************************************************************/
/********************************************************************/
void output(std::string operation, complex c1, complex c2) {
    double diff_r = c1.r - c2.r;
    double diff_i = c1.i - c2.i;
    printf("%s: \tdiff=(%g,%g), \tc1 = (%g,%g),  \tc2 = (%g,%g)\n",
	   operation.c_str(),
	   diff_r, diff_i,
	   c1.r, c1.i,
	   c2.r, c2.i);
    // Epsilon for double is 1e-16, check it 1.0 + 1e-16 = 1.0
    if (std::abs(diff_r) > 1e-16||std::abs(diff_i) > 1e-16) 
	//if (diff_r > 1e-15||diff_i > 1e-15) 
	printf("\n********\n\t\tWARNING!! Non-zero diff!!!\n********\n");
}
/********************************************************************/
/********************************************************************/
/********************************************************************/
void output(std::string operation, std::complex<double> c1, std::complex<double> c2) {
    double diff_r = c1.real() - c2.real();
    double diff_i = c1.imag() - c2.imag();
    printf("%s: \tdiff=(%g,%g), \tc1 = (%g,%g),  \tc2 = (%g,%g)\n",
	   operation.c_str(),
	   diff_r, diff_i,
	   c1.real(), c1.imag(),
	   c2.real(), c2.imag());
    // Epsilon for double is 1e-16, check it 1.0 + 1e-16 = 1.0
    if (std::abs(diff_r) > 1e-16||std::abs(diff_i) > 1e-16) 
	//if (diff_r > 1e-15||diff_i > 1e-15) 
	printf("\n********\n\t\tWARNING!! Non-zero diff!!!\n********\n");
}
/********************************************************************/
/********************************************************************/
/********************************************************************/
void output_double(std::string operation, double c1, double c2) {
    double diff = c1 - c2;
    printf("%s: \tdiff=(%g), \tc1 = (%g),  \tc2 = (%g)\n",
	   operation.c_str(),
	   diff, c1, c2);
    // Epsilon for double is 1e-16, check it 1.0 + 1e-16 = 1.0
    if (std::abs(diff) > 1e-16)
	printf("\n********\n\t\tWARNING!! Non-zero diff!!!\n********\n");
}
/********************************************************************/
/********************************************************************/
/********************************************************************/
int main() {
    double r = 1.412039487560983471903;
    complex a1 = {1.0012001301318988478, -1.0123131231241412987237},
	b1 = {0.42312341412412398237, 1.32312314124141232987987}, c1;
    std::complex<double> a2(a1.r, a1.i), b2(b1.r, b1.i), c2;    
    printf("Re(a1) = %g, Im(a1) = %g\n", a1.r, a1.i);    
    printf("Re(a2) = %g, Im(a2) = %g\n", a2.real(), a2.imag());    

    printf("Re(b1) = %g, Im(b1) = %g\n", b1.r, b1.i);    
    printf("Re(b2) = %g, Im(b2) = %g\n", b2.real(), b2.imag());
    printf("r = %g\n", r);

    c1 = Cadd(a1, b1);
    c2 = a2 + b2;
    output("Add", c1, c2);
    c1 = Cadd(b1, a1);
    c2 = b2 + a2;
    output("sAdd", c1, c2);

    c1 = Csub(a1, b1);
    c2 = a2 - b2;
    output("Sub", c1, c2);
    c1 = Csub(b1, a1);
    c2 = b2 - a2;
    output("sSub", c1, c2);

    c1 = Cmul(a1, b1);
    c2 = a2 * b2;
    output("Mul", c1, c2);
    c1 = Cmul(b1, a1);
    c2 = b2 * a2;
    output("sMul", c1, c2);

    c1 = RCmul(r, b1);
    c2 = r * b2;
    output("RCMul", c1, c2);

    c1 = Cdiv(a1, b1);
    c2 = a2 / b2;
    output("Div", c1, c2);
    c1 = Cdiv(b1, a1);
    c2 = b2 / a2;
    output("sDiv", c1, c2);

    c1 = Conjg(a1);
    c2 = std::conj(a2);
    output("Conj", c1, c2);
    c1 = Conjg(b1);
    c2 = std::conj(b2);
    output("sConj", c1, c2);

    c1 = Cinv(a1);
    c2 = 1.0/a2;
    output("Inv", c1, c2);
    c1 = Cinv(b1);
    c2 = 1.0/b2;
    output("sInv", c1, c2);

    output_double("Abs", Cabs(a1), std::abs(a2));
    output_double("sAbs", Cabs(b1), std::abs(b2));
    output_double("Arg", Carc(a1), std::arg(a2));
    output_double("sArg", Carc(b1), std::arg(b2));

    c1 = Cexp(a1);
    c2 = std::exp(a2);
    output("Exp", c1, c2);
    c1 = Cexp(b1);
    c2 = std::exp(b2);
    output("sExp", c1, c2);

    c1 = Clog(a1);
    c2 = std::log(a2);
    output("Log", c1, c2);
    c1 = Clog(b1);
    c2 = std::log(b2);
    output("sLog", c1, c2);

    c1 = Csqrt(a1);
    c2 = std::sqrt(a2);
    output("Sqrt", c1, c2);
    c1 = Csqrt(b1);
    c2 = std::sqrt(b2);
    output("sSqrt", c1, c2);

    c1 = Ccos(a1);
    c2 = std::cos(a2);
    output("Cos", c1, c2);
    c1 = Ccos(b1);
    c2 = std::cos(b2);
    output("sCos", c1, c2);

    c1 = Csin(a1);
    c2 = std::sin(a2);
    output("Sin", c1, c2);
    c1 = Csin(b1);
    c2 = std::sin(b2);
    output("sSin", c1, c2);

    c1 = Ctan(a1);
    c2 = std::tan(a2);
    output("Tan", c1, c2);
    c1 = Ctan(b1);
    c2 = std::tan(b2);
    output("sTan", c1, c2);

    c1 = Carc_cos(a1);
    c2 = std::acos(a2);
    output("aCos", c1, c2);
    c1 = Carc_cos(b1);
    c2 = std::acos(b2);
    output("saCos", c1, c2);

    c1 = Ccos(Carc_cos(b1));
    output("saCos_cos_b1", c1, b1);
    c2 = std::cos(std::acos(b2));
    output("saCos_cos_b2", c2, b2);


    c1 = Carc_sin(a1);
    c2 = std::asin(a2);
    output("aSin", c1, c2);
    c1 = Carc_sin(b1);
    c2 = std::asin(b2);
    output("saSin", c1, c2);

    c1 = Carc_tan(a1);
    c2 = std::atan(a2);
    output("aTan", c1, c2);
    c1 = Carc_tan(b1);
    c2 = std::atan(b2);
    output("saTan", c1, c2);

    // Hyperbolic
 
    c1 = Cch(a1);
    c2 = std::cosh(a2);
    output("Cosh", c1, c2);
    c1 = Cch(b1);
    c2 = std::cosh(b2);
    output("sCosh", c1, c2);

    c1 = Csh(a1);
    c2 = std::sinh(a2);
    output("Sinh", c1, c2);
    c1 = Csh(b1);
    c2 = std::sinh(b2);
    output("sSinh", c1, c2);

    c1 = Cth(a1);
    c2 = std::tanh(a2);
    output("Tanh", c1, c2);
    c1 = Cth(b1);
    c2 = std::tanh(b2);
    output("sTanh", c1, c2);

    c1 = Carc_ch(a1);
    c2 = std::acosh(a2);
    output("aCosh", c1, c2);
    c1 = Carc_ch(b1);
    c2 = std::acosh(b2);
    output("saCosh", c1, c2);

    c1 = Cch(Carc_ch(b1));
    output("saCosh_cosh_b1", c1, b1);
    c2 = std::cosh(std::acosh(b2));
    output("saCosh_cosh_b2", c2, b2);


    c1 = Carc_sh(a1);
    c2 = std::asinh(a2);
    output("aSinh", c1, c2);
    c1 = Carc_sh(b1);
    c2 = std::asinh(b2);
    output("saSinh", c1, c2);

    c1 = Carc_th(a1);
    c2 = std::atanh(a2);
    output("aTanh", c1, c2);
    c1 = Carc_th(b1);
    c2 = std::atanh(b2);
    output("saTanh", c1, c2);

    
    return 0;
}
