/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2025 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: MIT                                              */
/*----------------------------------------------------------------------------*/

/*
 *  libgrpp - a library for the evaluation of integrals over
 *            generalized relativistic pseudopotentials.
 *
 *  Copyright (C) 2021-2023 Alexander Oleynichenko
 */

#include <float.h> // required for LDBL_EPSILON
#include <math.h>  // required for fabsl()
#ifndef M_PI
#define M_PI 3.1415926535897932384626433
#endif

#include "grpp_specfunc.h"

/*
 * This code is taken from the website:
 * http://www.mymathlib.com/functions/dawsons_integral.html
 */

////////////////////////////////////////////////////////////////////////////////
// File: specfunc_dawson.c                                                    //
// Routine(s):                                                                //
//    Dawsons_Integral                                                        //
//    xDawsons_Integral                                                       //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  Description:                                                              //
//     Dawson's integral, Daw(x), is the integral from 0 to x of the          //
//     integrand:                                                             //
//                          exp(-x^2) exp(t^2) dt.                            //
//     I.e.                                                                   //
//                Daw(x) = exp(-x^2) * I[0,x] (exp(t^2)) dt,                  //
//     where I[0,x] indicates the integral from 0 to x.                       //
////////////////////////////////////////////////////////////////////////////////
#define Asymptotic_Expansion_Cutoff 50

//                         Externally Defined Routines                        //
static long double xChebyshev_Tn_Series(long double x, long double a[],
                                        int degree);

//                         Internally Defined Routines                        //
double libgrpp_Dawsons_Integral(double x);
long double xDawsons_Integral(long double x);

static long double Dawson_Power_Series(long double x);
static long double Dawson_Chebyshev_Expansion_1_175(long double x);
static long double Dawson_Chebyshev_Expansion_175_250(long double x);
static long double Dawson_Chebyshev_Expansion_250_325(long double x);
static long double Dawson_Chebyshev_Expansion_325_425(long double x);
static long double Dawson_Chebyshev_Expansion_425_550(long double x);
static long double Dawson_Chebyshev_Expansion_550_725(long double x);
static long double Dawson_Asymptotic_Expansion(long double x);

////////////////////////////////////////////////////////////////////////////////
// double libgrpp_Dawsons_Integral( double x ) //
//                                                                            //
//  Description:                                                              //
//     Dawson's integral, Daw(x), is the integral with integrand              //
//                          exp(-x^2) exp(t^2) dt                             //
//     where the integral extends from 0 to x.                                //
//                                                                            //
//  Arguments:                                                                //
//     double  x  The argument of Dawson's integral Daw().                    //
//                                                                            //
//  Return Value:                                                             //
//     The value of Dawson's integral, Daw(), evaluated at x.                 //
//                                                                            //
//  Example:                                                                  //
//     double y, x;                                                           //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Dawsons_Integral( x );                                             //
////////////////////////////////////////////////////////////////////////////////
double libgrpp_Dawsons_Integral(double x) {
  return (double)xDawsons_Integral((long double)x);
}

////////////////////////////////////////////////////////////////////////////////
// long double xDawsons_Integral( long double x )                             //
//                                                                            //
//  Description:                                                              //
//     Dawson's integral, Daw(x), is the integral with integrand              //
//                          exp(-x^2) exp(t^2) dt                             //
//     where the integral extends from 0 to x.                                //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of Dawson's integral Daw().               //
//                                                                            //
//  Return Value:                                                             //
//     The value of Dawson's integral, Daw(), evaluated at x.                 //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = xDawsons_Integral( x );                                            //
////////////////////////////////////////////////////////////////////////////////

long double xDawsons_Integral(long double x) {
  long double abs_x = fabsl(x);

  if (abs_x <= 1.0L)
    return Dawson_Power_Series(x);
  if (abs_x <= 1.75L)
    return Dawson_Chebyshev_Expansion_1_175(x);
  if (abs_x <= 2.50L)
    return Dawson_Chebyshev_Expansion_175_250(x);
  if (abs_x <= 3.25L)
    return Dawson_Chebyshev_Expansion_250_325(x);
  if (abs_x <= 4.25L)
    return Dawson_Chebyshev_Expansion_325_425(x);
  if (abs_x <= 5.50L)
    return Dawson_Chebyshev_Expansion_425_550(x);
  if (abs_x <= 7.25L)
    return Dawson_Chebyshev_Expansion_550_725(x);
  return Dawson_Asymptotic_Expansion(x);
}

////////////////////////////////////////////////////////////////////////////////
// static long double Dawson_Power_Series( long double x )                    //
//                                                                            //
//  Description:                                                              //
//     Evaluate Dawsons integral for -1 <= x <= 1 using the power series      //
//     expansion:                                                             //
//                  Daw(x) = x Sum [(-2x^2)^j / (2j + 1)!!,                   //
//     where the sum extends over j = 0, ... .                                //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of Dawson's integral Daw(), where         //
//                     where |x| <= 1.                                        //
//                                                                            //
//  Return Value:                                                             //
//     The value of Dawson's integral, Daw(), evaluated at x where |x| <= 1.  //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Dawson_Power_Series(x);                                            //
////////////////////////////////////////////////////////////////////////////////

static long double Dawson_Power_Series(long double x) {
  long double two_x2 = -2.0L * x * x;
  long double sum = 0.0;
  long double term = 1.0L;
  long double factorial = 1.0L;
  long double xn = 1.0L;
  const long double epsilon = LDBL_EPSILON / 2.0L;
  int y = 0;

  if (x == 0.0L)
    return 0.0L;
  do {
    sum += term;
    y += 1;
    factorial *= (long double)(y + y + 1);
    xn *= two_x2;
    term = xn / factorial;
  } while (fabsl(term) > epsilon * fabsl(sum));
  return x * sum;
}

////////////////////////////////////////////////////////////////////////////////
// static long double Dawson_Chebyshev_Expansion_1_175( long double x )       //
//                                                                            //
//  Description:                                                              //
//     Evaluate Dawsons integral for 1 <= |x| <= 1.75 using the Chebyshev     //
//     expansion of Daw(x) for x in the interval [1,1.75].                    //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of Dawson's integral Daw(), where         //
//                     where 1 <= |x| <= 1.75.                                //
//                                                                            //
//  Return Value:                                                             //
//     The value of Dawson's integral, Daw(), evaluated at x where            //
//     1 <= x <= 1.75.                                                        //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Dawson_Chebyshev_Expansion_1_175(x);                               //
////////////////////////////////////////////////////////////////////////////////

static long double Dawson_Chebyshev_Expansion_1_175(long double x) {
  static long double c[] = {
      +4.563960711239483142081e-1L,  -9.268566100670767619861e-2L,
      -7.334392170021420220239e-3L,  +3.379523740404396755124e-3L,
      -3.085898448678595090813e-4L,  -1.519846724619319512311e-5L,
      +4.903955822454009397182e-6L,  -2.106910538629224721838e-7L,
      -2.930676220603996192089e-8L,  +3.326790071774057337673e-9L,
      +3.335593457695769191326e-11L, -2.279104036721012221982e-11L,
      +7.877561633156348806091e-13L, +9.173158167107974472228e-14L,
      -7.341175636102869400671e-15L, -1.763370444125849029511e-16L,
      +3.792946298506435014290e-17L, -4.251969162435936250171e-19L,
      -1.358295820818448686821e-19L, +5.268740962820224108235e-21L,
      +3.414939674304748094484e-22L};

  static const int degree = sizeof(c) / sizeof(long double) - 1;
  static const long double midpoint = (1.75L + 1.0L) / 2.0L;
  static const long double half_length = (1.75L - 1.0L) / 2.0L;
  long double daw =
      xChebyshev_Tn_Series((fabsl(x) - midpoint) / half_length, c, degree);
  return (x > 0.0L) ? daw : -daw;
}

////////////////////////////////////////////////////////////////////////////////
// static long double Dawson_Chebyshev_Expansion_175_250( long double x )     //
//                                                                            //
//  Description:                                                              //
//     Evaluate Dawsons integral for 1.75 <= |x| <= 2.50 using the Chebyshev  //
//     expansion of Daw(x) for x in the interval [1.75,2.50].                 //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of Dawson's integral Daw(), where         //
//                     where 1.75 <= |x| <= 2.50.                             //
//                                                                            //
//  Return Value:                                                             //
//     The value of Dawson's integral, Daw(), evaluated at x where            //
//     1.75 <= x <= 2.50.                                                     //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Dawson_Chebyshev_Expansion_175_250(x);                             //
////////////////////////////////////////////////////////////////////////////////

static long double Dawson_Chebyshev_Expansion_175_250(long double x) {
  static long double c[] = {
      +2.843711194548592808550e-1L,  -6.791774139166808940530e-2L,
      +6.955211059379384327814e-3L,  -2.726366582146839486784e-4L,
      -6.516682485087925163874e-5L,  +1.404387911504935155228e-5L,
      -1.103288540946056915318e-6L,  -1.422154597293404846081e-8L,
      +1.102714664312839585330e-8L,  -8.659211557383544255053e-10L,
      -8.048589443963965285748e-12L, +6.092061709996351761426e-12L,
      -3.580977611213519234324e-13L, -1.085173558590137965737e-14L,
      +2.411707924175380740802e-15L, -7.760751294610276598631e-17L,
      -6.701490147030045891595e-18L, +6.350145841254563572100e-19L,
      -2.034625734538917052251e-21L, -2.260543651146274653910e-21L,
      +9.782419961387425633151e-23L};

  static const int degree = sizeof(c) / sizeof(long double) - 1;
  static const long double midpoint = (2.50L + 1.75L) / 2.0L;
  static const long double half_length = (2.50L - 1.75L) / 2.0;
  long double daw =
      xChebyshev_Tn_Series((fabsl(x) - midpoint) / half_length, c, degree);
  return (x > 0.0L) ? daw : -daw;
}

////////////////////////////////////////////////////////////////////////////////
// static long double Dawson_Chebyshev_Expansion_250_325( long double x )     //
//                                                                            //
//  Description:                                                              //
//     Evaluate Dawsons integral for 2.5 <= |x| <= 3.25 using the Chebyshev   //
//     expansion of Daw(x) for x in the interval [2.50,3.25].                 //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of Dawson's integral Daw(), where         //
//                     where 2.50 <= |x| <= 3.25.                             //
//                                                                            //
//  Return Value:                                                             //
//     The value of Dawson's integral, Daw(), evaluated at x where            //
//     2.50 <= x <= 3.25.                                                     //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Dawson_Chebyshev_Expansion_250_325(x);                             //
////////////////////////////////////////////////////////////////////////////////

static long double Dawson_Chebyshev_Expansion_250_325(long double x) {
  static long double c[] = {
      +1.901351274204578126827e-1L,  -3.000575522193632460118e-2L,
      +2.672138524890489432579e-3L,  -2.498237548675235150519e-4L,
      +2.013483163459701593271e-5L,  -8.454663603108548182962e-7L,
      -8.036589636334016432368e-8L,  +2.055498509671357933537e-8L,
      -2.052151324060186596995e-9L,  +8.584315967075483822464e-11L,
      +5.062689357469596748991e-12L, -1.038671167196342609090e-12L,
      +6.367962851860231236238e-14L, +3.084688422647419767229e-16L,
      -3.417946142546575188490e-16L, +2.311567730100119302160e-17L,
      -6.170132546983726244716e-20L, -9.133176920944950460847e-20L,
      +5.712092431423316128728e-21L, +1.269641078369737220790e-23L,
      -2.072659711527711312699e-23L};

  static const int degree = sizeof(c) / sizeof(long double) - 1;
  static const long double midpoint = (3.25L + 2.50L) / 2.0L;
  static const long double half_length = (3.25L - 2.50L) / 2.0L;
  long double daw =
      xChebyshev_Tn_Series((fabsl(x) - midpoint) / half_length, c, degree);
  return (x > 0.0L) ? daw : -daw;
}

////////////////////////////////////////////////////////////////////////////////
// static long double Dawson_Chebyshev_Expansion_325_425( long double x )     //
//                                                                            //
//  Description:                                                              //
//     Evaluate Dawsons integral for 3.25 <= |x| <= 4.25 using the Chebyshev  //
//     expansion of Daw(x) for x in the interval [3.25,4.75].                 //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of Dawson's integral Daw(), where         //
//                     where 3.25 <= |x| <= 4.25.                             //
//                                                                            //
//  Return Value:                                                             //
//     The value of Dawson's integral, Daw(), evaluated at x where            //
//     3.25 <= x <= 4.25.                                                     //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Dawson_Chebyshev_Expansion_325_425(x);                             //
////////////////////////////////////////////////////////////////////////////////

static long double Dawson_Chebyshev_Expansion_325_425(long double x) {
  static long double c[] = {
      +1.402884974484995678749e-1L,  -2.053975371995937033959e-2L,
      +1.595388628922920119352e-3L,  -1.336894584910985998203e-4L,
      +1.224903774178156286300e-5L,  -1.206856028658387948773e-6L,
      +1.187997233269528945503e-7L,  -1.012936061496824448259e-8L,
      +5.244408240062370605664e-10L, +2.901444759022254846562e-11L,
      -1.168987502493903926906e-11L, +1.640096995420504465839e-12L,
      -1.339190668554209618318e-13L, +3.643815972666851044790e-15L,
      +6.922486581126169160232e-16L, -1.158761251467106749752e-16L,
      +8.164320395639210093180e-18L, -5.397918405779863087588e-20L,
      -5.052069908100339242896e-20L, +5.322512674746973445361e-21L,
      -1.869294542789169825747e-22L};

  static const int degree = sizeof(c) / sizeof(long double) - 1;
  // static const long double lower_bound = 3.25L;
  // static const long double upper_bound = 4.25L;
  static const long double midpoint = (4.25L + 3.25L) / 2.0L;
  static const long double half_length = (4.25L - 3.25L) / 2.0L;
  long double daw =
      xChebyshev_Tn_Series((fabsl(x) - midpoint) / half_length, c, degree);
  return (x > 0.0L) ? daw : -daw;
}

////////////////////////////////////////////////////////////////////////////////
// static long double Dawson_Chebyshev_Expansion_425_550( long double x )     //
//                                                                            //
//  Description:                                                              //
//     Evaluate Dawsons integral for 4.25 <= |x| <= 5.50 using the Chebyshev  //
//     expansion of Daw(x) for x in the interval [4.25,5.50].                 //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of Dawson's integral Daw(), where         //
//                     where 4.25 <= |x| <= 5.50.                             //
//                                                                            //
//  Return Value:                                                             //
//     The value of Dawson's integral, Daw(), evaluated at x where            //
//     4.25 <= x <= 5.50.                                                     //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Dawson_Chebyshev_Expansion_425_550(x);                             //
////////////////////////////////////////////////////////////////////////////////

static long double Dawson_Chebyshev_Expansion_425_550(long double x) {
  static long double c[] = {
      +1.058610209741581514157e-1L,  -1.429297757627935191694e-2L,
      +9.911301703835545472874e-4L,  -7.079903107876049846509e-5L,
      +5.229587914675267516134e-6L,  -4.016071345964089296212e-7L,
      +3.231734714422926453741e-8L,  -2.752870944370338482109e-9L,
      +2.503059741885009530630e-10L, -2.418699000594890423278e-11L,
      +2.410158905786160001792e-12L, -2.327254341132174000949e-13L,
      +1.958284411563056492727e-14L, -1.099893145048991004460e-15L,
      -2.959085292526991317697e-17L, +1.966366179276295203082e-17L,
      -3.314408783993662492621e-18L, +3.635520318133814622089e-19L,
      -2.550826919215104648800e-20L, +3.830090587178262542288e-22L,
      +1.836693763159216122739e-22L};

  static const int degree = sizeof(c) / sizeof(long double) - 1;
  static const long double midpoint = (5.50L + 4.25L) / 2.0L;
  static const long double half_length = (5.50L - 4.25L) / 2.0L;
  long double daw =
      xChebyshev_Tn_Series((fabsl(x) - midpoint) / half_length, c, degree);
  return (x > 0.0L) ? daw : -daw;
}

////////////////////////////////////////////////////////////////////////////////
// static long double Dawson_Chebyshev_Expansion_550_725( long double x )     //
//                                                                            //
//  Description:                                                              //
//     Evaluate Dawsons integral for 5.50 <= |x| <= 7.25 using the Chebyshev  //
//     expansion of Daw(x) for x in the interval [5.50,7.25].                 //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of Dawson's integral Daw(), where         //
//                     where 5.50 <= |x| <= 7.25.                             //
//                                                                            //
//  Return Value:                                                             //
//     The value of Dawson's integral, Daw(), evaluated at x where            //
//     5.50 <= x <= 7.25.                                                     //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Dawson_Chebyshev_Expansion_550_725(x);                             //
////////////////////////////////////////////////////////////////////////////////

static long double Dawson_Chebyshev_Expansion_550_725(long double x) {
  static long double c[] = {+8.024637207807814739314e-2L,
                            -1.136614891549306029413e-2L,
                            +8.164249750628661856014e-4L,
                            -5.951964778701328943018e-5L,
                            +4.407349502747483429390e-6L,
                            -3.317746826184531133862e-7L,
                            +2.541483569880571680365e-8L,
                            -1.983391157250772649001e-9L,
                            +1.579050614491277335581e-10L,
                            -1.284592098551537518322e-11L,
                            +1.070070857004674207604e-12L,
                            -9.151832297362522251950e-14L,
                            +8.065447314948125338081e-15L,
                            -7.360105847607056315915e-16L,
                            +6.995966000187407197283e-17L,
                            -6.964349343411584120055e-18L,
                            +7.268789359189778223225e-19L,
                            -7.885125241947769024019e-20L,
                            +8.689022564130615225208e-21L,
                            -9.353211304381231554634e-22L +
                                9.218280404899298404756e-23L};

  static const int degree = sizeof(c) / sizeof(long double) - 1;
  // static const long double lower_bound = 5.50L;
  // static const long double upper_bound = 7.25L;
  static const long double midpoint = (7.25L + 5.50L) / 2.0L;
  static const long double half_length = (7.25L - 5.50L) / 2.0L;
  long double daw =
      xChebyshev_Tn_Series((fabsl(x) - midpoint) / half_length, c, degree);
  return (x > 0.0L) ? daw : -daw;
}

////////////////////////////////////////////////////////////////////////////////
// static long double Dawson_Asymptotic_Expansion( long double x )            //
//                                                                            //
//  Description:                                                              //
//     For a large magnitude of the argument x, Dawson's integral can be      //
//     expressed as the asymptotic series                                     //
//     Daw(x) ~ (1/2x) [ 1 + 1 / (2x^2) + ... + (2j - 1)!! / (2x^2)^j + ... ] //
//                                                                            //
//  Arguments:                                                                //
//     long double  x  The argument of Dawson's integral Daw(), where         //
//                     |x| > 7.                                               //
//                                                                            //
//  Return Value:                                                             //
//     The value of Dawson's integral, Daw(), evaluated at x where |x| > 7.   //
//                                                                            //
//  Example:                                                                  //
//     long double y, x;                                                      //
//                                                                            //
//     ( code to initialize x )                                               //
//                                                                            //
//     y = Dawson_Asymptotic_Expansion( x );                                  //
////////////////////////////////////////////////////////////////////////////////
static long double Dawson_Asymptotic_Expansion(long double x) {
  long double term[Asymptotic_Expansion_Cutoff + 1];
  long double x2 = x * x;
  long double two_x = x + x;
  long double two_x2 = x2 + x2;
  long double xn = two_x2;
  long double Sn = 0.0L;
  long double factorial = 1.0L;
  int n;

  term[0] = 1.0L;
  term[1] = 1.0L / xn;
  for (n = 2; n <= Asymptotic_Expansion_Cutoff; n++) {
    xn *= two_x2;
    factorial *= (long double)(n + n - 1);
    term[n] = factorial / xn;
    if (term[n] < LDBL_EPSILON / 2.0L)
      break;
  }

  if (n > Asymptotic_Expansion_Cutoff)
    n = Asymptotic_Expansion_Cutoff;
  for (; n >= 0; n--)
    Sn += term[n];
  return Sn / two_x;
}

////////////////////////////////////////////////////////////////////////////////
// File: xchebyshev_Tn_series.c                                               //
// Routine(s):                                                                //
//    xChebyshev_Tn_Series                                                    //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// long double xChebyshev_Tn_Series(long double x, long double a[],int degree)//
//                                                                            //
//  Description:                                                              //
//     This routine uses Clenshaw's recursion algorithm to evaluate a given   //
//     polynomial p(x) expressed as a linear combination of Chebyshev         //
//     polynomials of the first kind, Tn, at a point x,                       //
//       p(x) = a[0] + a[1]*T[1](x) + a[2]*T[2](x) + ... + a[deg]*T[deg](x).  //
//                                                                            //
//     Clenshaw's recursion formula applied to Chebyshev polynomials of the   //
//     first kind is:                                                         //
//     Set y[degree + 2] = 0, y[degree + 1] = 0, then for k = degree, ..., 1  //
//     set y[k] = 2 * x * y[k+1] - y[k+2] + a[k].  Finally                    //
//     set y[0] = x * y[1] - y[2] + a[0].  Then p(x) = y[0].                  //
//                                                                            //
//  Arguments:                                                                //
//     long double x                                                          //
//        The point at which to evaluate the polynomial.                      //
//     long double a[]                                                        //
//        The coefficients of the expansion in terms of Chebyshev polynomials,//
//        i.e. a[k] is the coefficient of T[k](x).  Note that in the calling  //
//        routine a must be defined double a[N] where N >= degree + 1.        //
//     int    degree                                                          //
//        The degree of the polynomial p(x).                                  //
//                                                                            //
//  Return Value:                                                             //
//     The value of the polynomial at x.                                      //
//     If degree is negative, then 0.0 is returned.                           //
//                                                                            //
//  Example:                                                                  //
//     long double x, a[N], p;                                                //
//     int    deg = N - 1;                                                    //
//                                                                            //
//     ( code to initialize x, and a[i] i = 0, ... , a[deg] )                 //
//                                                                            //
//     p = xChebyshev_Tn_Series(x, a, deg);                                   //
////////////////////////////////////////////////////////////////////////////////

static long double xChebyshev_Tn_Series(long double x, long double a[],
                                        int degree) {
  long double yp2 = 0.0L;
  long double yp1 = 0.0L;
  long double y = 0.0L;
  long double two_x = x + x;
  int k;

  // Check that degree >= 0.  If not, then return 0. //

  if (degree < 0)
    return 0.0L;

  // Apply Clenshaw's recursion save the last iteration. //

  for (k = degree; k >= 1; k--, yp2 = yp1, yp1 = y)
    y = two_x * yp1 - yp2 + a[k];

  // Now apply the last iteration and return the result. //

  return x * yp1 - yp2 + a[0];
}
