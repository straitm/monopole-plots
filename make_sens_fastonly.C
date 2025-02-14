{
  printf("-3.75 1.5e-11 1.5e-11\n");

  const double x = 910.55, y = 910.55, z = 233.17;

  const double Aproj = (x + y + z)/2 * 100 * 100;

  // 63% efficiency, 8.0 years live-time
  const double baselimit = log(10)/ ( Aproj * 8*3.15e7 * 2*M_PI * 0.63 );

  const int nbase = 42;

  double basebeta[nbase] = {
    0.000105,
    0.000125,
    0.000155,
    0.000195,
    0.000255,
    0.000315,
    0.000395,
    0.000505,
    0.000635,
    0.000795,
    0.001005,
    0.001255,
    0.001585,
    0.001995,
    0.002515,
    0.003165,
    0.003985,
    0.005015,
    0.006305,
    0.007945,
    0.010005,
    0.012584,
    0.015844,
    0.019954,
    0.025124,
    0.031623,
    0.039813,
    0.050122,
    0.063092,
    0.079431,
    0.1,
    0.125889,
    0.158487,
    0.199525,
    0.251192,
    0.316229,
    0.398105,
    0.50119,
    0.630953,
    0.794325,
    1
  };

  double baseeff[nbase] = {
    0.0000,
    0.0000,
    0.0000,
    0.0000,
    0.0000,
    0.0000,
    0.0000003,
    0.000003,
    0.00003,
    0.0003,
    0.0054,
    0.0243,
    0.0552,
    0.1045,
    0.1656,
    0.2412,
    0.3177,
    0.3886,
    0.4566,
    0.5135,
    0.5498,
    0.5709,
    0.5817,
    0.5900,
    0.5965,
    0.5978,
    0.5999,
    0.6068,
    0.6108,
    0.6143,
    0.6143,
    0.6153,
    0.6187,
    0.6175,
    0.6179,
    0.6166,
    0.6148,
    0.6132,
    0.6155,
    0.6288,
    0.6288
  };

  TGraph base(nbase, basebeta, baseeff);

  const double crosstalkeff = (4519+4487+4596.)/3./6702;

  const int ndelta = 7;
  double deltabeta[ndelta] = {
    0,
    0.75,
    0.8,
    0.85,
    0.9,
    0.99,
    1
  };

  double deltaeff[ndelta] = {
    1,
    1,
    0.57/0.63,
    0.53/0.63,
    0.45/0.63,
    0.03/0.63,
    0
  };

  TGraph deltarayeff(ndelta, deltabeta, deltaeff);

  const double minlogbeta = -3.75;
  const double maxlogbeta = -0.022;

  const int N = 30;
  for(int i = 0; i <= N; i++){
    const double logbeta = minlogbeta + i*(maxlogbeta-minlogbeta)/N;
    const double beta = pow(10., logbeta);

    if(base.Eval(beta) <= 0) continue;

    const double limit = baselimit / base.Eval(beta)
      / (beta > 0.68? crosstalkeff: 1) / deltarayeff.Eval(beta);

    printf("%f %g %g\n", logbeta, limit, limit/2);
  }

  printf("-0.022 1.125e-11 1.125e-11\n");
}
