void make_sens_slowonly()
{

  const double x = 910.55, y = 910.55, z = 233.17;

  const double Aproj = (x + y + z)/2 * 100 * 100;

  // Limit for 100% efficiency for heavy monopoles
  const double livetime = 257.870154e6
  #if 0
    * 13/8. // Thirteen years instead of eight
  #endif
  ;
  const double baselimit = log(10)/ ( Aproj * livetime * 4*M_PI);

  const int n = 21;

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
    -2.2,
    -2.1,
    -2.0,

    -1.99,
    -1.98,
    -1.97
  };

  double eff[n] = {
    0.0023,
    0.1549,
    0.3408,
    0.4372,
    0.4811,
    0.4888,
    0.5032,
    0.5051,
    0.5081,
    0.5104,
    0.5149,
    0.5145,
    0.5078,
    0.5029,
    0.4524,
    0.2964,
    0.1754,
    0.0517,

    0.005,
    0.0005,
    0.00005,
  };

  double logeff[n] = { 0 };

  for(int i = 0; i < n; i++) logeff[i] = log10(eff[i]);

  TSpline3 spline("whatever", logbeta, eff, n);
  TGraph graph(n, logbeta, eff);

  const double minlogbeta = logbeta[0];
  const double maxlogbeta = logbeta[n-1];

  const int npoints = 200;

  printf("%f 1.125e-11 1.125e-11\n", logbeta[0]);

  for(int i = 0; i <= npoints; i++){
    const double lbeta = minlogbeta + (maxlogbeta - minlogbeta)*double(i)/npoints;

    // Hack to smooth out the curve
    const double the_eff = spline.Eval(lbeta);
    const double the_effg = graph.Eval(lbeta);

    const double limit = baselimit / (lbeta < -2.5? the_eff: the_effg);

    printf("%f %g %g\n", lbeta, limit*2, limit);
  }

  printf("%f 1.125e-11 1.125e-11\n", logbeta[n-1]);
}
