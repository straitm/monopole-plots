void make_sens_slowonly()
{
  const double plottop = 1.5e-11;

  const int n = 24;

  double logbeta[n] = {
    -3.7,
    -3.6,
    -3.5,
    -3.4,
    -3.3,
    -3.2,
    -3.1,
    -3.0,
    -2.9,
    -2.8,
    -2.7,
    -2.6,
    -2.5,
    -2.4,

    -2.3,
    -2.2218,
    -2.2,
    -2.1549,
    -2.1,
    -2.0457,
    -2.0,

    -1.99,
    -1.98,
    -1.97
  };

  double limit4pi[n] = {
    1.0277e-13,
    1.5143e-15,
    4.6744e-16,
    2.4261e-16,
    1.823e-16,
    1.6907e-16,
    1.6372e-16,
    1.6049e-16,
    1.6033e-16,
    1.5863e-16,
    1.5942e-16,
    1.5707e-16,
    1.5912e-16,
    1.6035e-16,
    1.8081e-16,
    3.9765e-16,
    4.7614e-16,
    6.9744e-16,
    1.1907e-15,
    2.2337e-15,
    6.3414e-15,

    6.3414e-14,
    6.3414e-13,
    6.3414e-12,
  };

  double limit1pi[n] = {
    1.2872e-12,
    6.4971e-15,
    2.034e-15,
    1.0103e-15,
    7.7187e-16,
    6.7797e-16,
    6.6578e-16,
    6.3281e-16,
    6.6146e-16,
    6.258e-16,
    6.597e-16,
    6.48e-16,
    6.4164e-16,
    6.4401e-16,
    7.0432e-16,
    8.3876e-16,
    9.889e-16,
    1.3999e-15,
    2.3832e-15,
    4.5455e-15,
    7.7471e-15,

    7.7471e-14,
    7.7471e-13,
    7.7471e-12,
  };

  double loglimit4pi[n] = { 0 }, loglimit1pi[n] = { 0 };

  for(int i = 0; i < n; i++){
    loglimit4pi[i] = log10(limit4pi[i]);
    loglimit1pi[i] = log10(limit1pi[i]);
  }

  TSpline3 spline1pi("whatever", logbeta, loglimit1pi, n);
  TGraph graph1pi(n, logbeta, loglimit1pi);

  TSpline3 spline4pi("whatever", logbeta, loglimit4pi, n);
  TGraph graph4pi(n, logbeta, loglimit4pi);

  const double minlogbeta = logbeta[0];
  const double maxlogbeta = logbeta[n-1];

  const int npoints = 200;

  printf("%f %g %g\n", logbeta[0], plottop, plottop);

  double final_analysis_mult = 0;
  {
    const double livetime_2025_analysis = 236.977204e6;

    const double additional_livetime =
      // From 2024-10-13 to 2028-12-31, optmistically assuming we get another
      // year of data because of LBNF delays and/or DARPA funding.
      // Assume 90% uptime for good data.
      + (1.5/12. + 4.0) * 365 * 86400 * 0.9;

     final_analysis_mult =
       (additional_livetime + livetime_2025_analysis)/livetime_2025_analysis;
  }

  for(int i = 0; i <= npoints; i++){
    const double lbeta = minlogbeta + (maxlogbeta - minlogbeta)*double(i)/npoints;

    // Hack to smooth out the curve
    const double the_limit1pi = pow(10, spline1pi.Eval(lbeta));
    const double the_limitg1pi = pow(10, graph1pi.Eval(lbeta));

    const double the_limit4pi = pow(10, spline4pi.Eval(lbeta));
    const double the_limitg4pi = pow(10, graph4pi.Eval(lbeta));

    const double limit1pi = lbeta < -2.5? the_limit1pi: the_limitg1pi;

    const double limit4pi = lbeta < -2.5? the_limit4pi: the_limitg4pi;

    //#define FINAL
    #ifdef FINAL
      printf("%f %g %g\n", lbeta,
             limit1pi/final_analysis_mult, limit4pi/final_analysis_mult);
    #else
      printf("%f %g %g\n", lbeta,
             limit1pi, limit4pi);
    #endif
  }

  printf("%f %g %g\n", logbeta[n-1], plottop, plottop);
}
