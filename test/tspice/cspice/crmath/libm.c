/* The transcendental functions of the C library that SPICELIB and the test families
   call, from CORE-MATH (commit b464d7d0280f4c5091b8d2fff80b6cb03ab1d37f, MIT license,
   see LICENSE). A correctly rounded result is the same on every platform; those of the
   C library differ between its versions in the last bit. math.h is not included: these
   definitions replace its functions at the link. */
double cr_acos(double), cr_asin(double), cr_atan(double), cr_atan2(double, double),
    cr_cos(double), cr_cosh(double), cr_exp(double), cr_log(double),
    cr_pow(double, double), cr_sin(double), cr_sinh(double), cr_tan(double),
    cr_tanh(double);
void cr_sincos(double, double *, double *);

double acos(double x) { return cr_acos(x); }
double asin(double x) { return cr_asin(x); }
double atan(double x) { return cr_atan(x); }
double atan2(double y, double x) { return cr_atan2(y, x); }
double cos(double x) { return cr_cos(x); }
double cosh(double x) { return cr_cosh(x); }
double exp(double x) { return cr_exp(x); }
double log(double x) { return cr_log(x); }
double pow(double x, double y) { return cr_pow(x, y); }
double sin(double x) { return cr_sin(x); }
double sinh(double x) { return cr_sinh(x); }
double tan(double x) { return cr_tan(x); }
double tanh(double x) { return cr_tanh(x); }
void sincos(double x, double *s, double *c) { cr_sincos(x, s, c); }
