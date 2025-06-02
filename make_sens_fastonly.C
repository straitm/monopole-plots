void make_sens_fastonly()
{
  printf("-3.75 1.5e-11 1.5e-11\n");

  const double x = 910.55, y = 910.55, z = 233.17;

  const double Aproj = (x + y + z)/2 * 100 * 100;

  // Limit for 100% efficiency for heavy monopoles
  const double livetime = 234498017.706;

  const double baselimit = log(10)/ ( Aproj * livetime * 4*M_PI);

  const int nbase = 43;

  double basebeta[nbase] = {
    1e-5,
    0.000105,
    0.000125,
    0.000155,
    0.000195,
    0.000255,
    0.000315,
    0.000395,
    0.000505,
    0.000635,
    0.000700,
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
    1e-10,
    1e-10,
    1e-10,
    1e-10,
    1e-10,
    1e-10,
    1e-10,
    1e-10,
    1e-10,
    0.0000003,
    0.00003,
    0.0003,
    0.0053,
    0.0236,
    0.0530,
    0.1007,
    0.1588,
    0.2308,
    0.3034,
    0.3708,
    0.4362,
    0.4898,
    0.5244,
    0.5445,
    0.5536,
    0.5607,
    0.5678,
    0.5679,
    0.5703,
    0.5763,
    0.5797,
    0.5841,
    0.5844,
    0.5845,
    0.5875,
    0.5870,
    0.5882,
    0.5867,
    0.5857,
    0.5831,
    0.5855,
    0.5941,
    0.5941,
  };

  for(int i = 0; i < nbase; i++){
    basebeta[i] = log10(basebeta[i]);
    if(baseeff[i] > 0) baseeff[i] = log10(baseeff[i]);
  }

  TGraph base(nbase, basebeta, baseeff);

  const double crosstalkeff = (4519+4487+4596.)/3./6702.;

  const int ndelta = 7;
  double deltabeta[ndelta] = {
    1e-5,
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
    0.57/pow(10, baseeff[nbase-1]),
    0.53/pow(10, baseeff[nbase-1]),
    0.45/pow(10, baseeff[nbase-1]),
    0.03/pow(10, baseeff[nbase-1]),
    1e-8
  };

  for(int i = 0; i < ndelta; i++){
    deltabeta[i] = log10(deltabeta[i]);
    if(deltaeff[i] > 0) deltaeff[i] = log10(deltaeff[i]);
  }

  TGraph * deltarayeff = new TGraph (ndelta, deltabeta, deltaeff);

  const double minlogbeta = -3.75;
  const double maxlogbeta = -0.022;

  const int N = 400;
  for(int i = 0; i <= N; i++){
    const double logbeta = minlogbeta + i*(maxlogbeta-minlogbeta)/N;
    const double beta = pow(10., logbeta);

    const double beff = pow(10, base.Eval(logbeta));
    const double crosseff = beta > 0.68? crosstalkeff:
                            beta > 0.64? (beta-0.64)/(0.68-0.64) * (crosstalkeff-1) + 1
                            :1;
    const double deltaeff = pow(10, deltarayeff->Eval(logbeta));

    const double limit = base.Eval(logbeta) == 0? 2.125e-11:
      baselimit / beff / crosseff / deltaeff;

    const double highmass_solidangle = 4*M_PI;
    const double lowmass_solidangle = 1*M_PI;

    printf("%f %g %g\n", logbeta, limit*highmass_solidangle/lowmass_solidangle, limit);
    #if 0
      printf("%f %g %g  %f %g %g %g\n", logbeta, limit, limit/2, beta, beff, crosseff, deltaeff);
    #endif
  }

  printf("-0.022 1.125e-11 1.125e-11\n");

  deltarayeff->Draw("ap*");
}
