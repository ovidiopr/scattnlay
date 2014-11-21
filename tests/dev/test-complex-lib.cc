#include <complex>
#include <string>
#include <stdio.h>
#include "../../ucomplex.h"
void output(std::string operation, complex c1, std::complex<double> c2) {
    double diff_r = c1.r - c2.real();
    double diff_i = c1.i - c2.imag();
    printf("%s: \tdiff=(%g,%g), \tc1 = (%g,%g),  \tc2 = (%g,%g)\n",
	   operation.c_str(),
	   diff_r, diff_i,
	   c1.r, c1.i,
	   c2.real(), c2.imag());
    // Epsilon for double is 1e-16, check it 1.0 + 1e-16 = 1.0
    if (diff_r > 1e-16||diff_i > 1e-16) printf("WARNING!! Non-zero diff!!!\n");
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
    if (diff > 1e-16) printf("WARNING!! Non-zero diff!!!\n");
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
    c2 = conj(a2);
    output("Conj", c1, c2);
    c1 = Conjg(b1);
    c2 = conj(b2);
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

    c1 = Cinv(a1);
    c2 = 1.0/a2;
    output("Inv", c1, c2);

    
    return 0;
}
