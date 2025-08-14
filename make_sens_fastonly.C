void make_sens_fastonly()
{
  const double plottop = 1.5e-11;

  printf("-3.75 %g %g\n", plottop, plottop);

  const int nbase = 43;

  double basebeta[nbase] = {
         1e-5, // -5. [0]
     0.000105, // -3.979 [1]
     0.000125, // -3.903 [2]
     0.000155, // -3.81 [3]
     0.000195, // -3.71 [4]
     0.000255, // -3.593 [5]
     0.000315, // -3.502 [6]
     0.000395, // -3.403 [7]
     0.000505, // -3.297 [8]
     0.000635, // -3.197 [9]
     0.000700, // -3.155 [10]

     0.000795, // -3.1 [11]
     0.001005, // -3.0 [12]
     0.001255, // -2.9 [13]
     0.001585, // -2.8 [14]
     0.001995, // -2.7 [15]
     0.002515, // -2.6 [16]
     0.003165, // -2.5 [17]
     0.003985, // -2.4 [18]
     0.005015, // -2.3 [19]
     0.006305, // -2.2 [20]
     0.007945, // -2.1 [21]
     0.010005, // -2. [22]
     0.012584, // -1.9 [23]
     0.015844, // -1.8 [24]
     0.019954, // -1.7 [25]
     0.025124, // -1.6 [26]
     0.031623, // -1.5 [27]
     0.039813, // -1.4 [28]
     0.050122, // -1.3 [29]
     0.063092, // -1.2 [30]
     0.079431, // -1.1 [31]
          0.1, // -1. [32]
     0.125889, // -0.9 [33]
     0.158487, // -0.8 [34]
     0.199525, // -0.7 [35]
     0.251192, // -0.6 [36]
     0.316229, // -0.5 [37]
     0.398105, // -0.4 [38]
      0.50119, // -0.3 [39]
     0.630953, // -0.2 [40]
     0.794325, // -0.1 [41]

             1 // 0. [42]
  };

  double baselimit4pi[nbase] = {
    plottop,  // [0]
    plottop,  // [1]
    plottop,  // [2]
    plottop,  // [3]
    plottop,  // [4]
    plottop,  // [5]
    plottop,  // [6]
    plottop,  // [7]
    plottop,  // [8]
    plottop,  // [9]
    plottop,  // [10]

    3.808e-13,  // [11]
    1.448e-14,  // [12]
    3.117e-15,  // [13]
    1.384e-15,  // [14]
    7.227e-16,  // [15]
    4.555e-16,  // [16]
    3.119e-16,  // [17]
    2.382e-16,  // [18]
    1.958e-16,  // [19]
    1.673e-16,  // [20]
    1.494e-16,  // [21]
    1.401e-16,  // [22]
    1.350e-16,  // [23]
    1.333e-16,  // [24]
    1.316e-16,  // [25]
    1.301e-16,  // [26]
    1.303e-16,  // [27]
    1.299e-16,  // [28]
    1.286e-16,  // [29]
    1.278e-16,  // [30]
    1.269e-16,  // [31]
    1.271e-16,  // [32]
    1.269e-16,  // [33]
    1.264e-16,  // [34]
    1.266e-16,  // [35]
    1.263e-16,  // [36]
    1.268e-16,  // [37]
    1.272e-16,  // [38]
    1.276e-16,  // [39]
    1.270e-16,  // [40]
    1.254e-16,  // [41]

    1.254e-16,  // [42]
  };

  double baselimit1pi[nbase] = {
    plottop,  // [0]
    plottop,  // [1]
    plottop,  // [2]
    plottop,  // [3]
    plottop,  // [4]
    plottop,  // [5]
    plottop,  // [6]
    plottop,  // [7]
    plottop,  // [8]
    plottop,  // [9]
    plottop,  // [10]

    2.473e-12,  // [11]
    8.019e-14,  // [12]
    1.481e-14,  // [13]
    5.622e-15,  // [14]
    3.000e-15,  // [15]
    1.759e-15,  // [16]
    1.179e-15,  // [17]
    8.929e-16,  // [18]
    7.451e-16,  // [19]
    6.374e-16,  // [20]
    5.674e-16,  // [21]
    5.307e-16,  // [22]
    5.113e-16,  // [23]
    5.030e-16,  // [24]
    4.995e-16,  // [25]
    4.906e-16,  // [26]
    4.959e-16,  // [27]
    4.973e-16,  // [28]
    4.884e-16,  // [29]
    4.838e-16,  // [30]
    4.822e-16,  // [31]
    4.803e-16,  // [32]
    4.769e-16,  // [33]
    4.845e-16,  // [34]
    4.808e-16,  // [35]
    4.778e-16,  // [36]
    4.798e-16,  // [37]
    4.848e-16,  // [38]
    4.871e-16,  // [39]
    4.868e-16,  // [40]
    4.694e-16,  // [41]

    4.694e-16,  // [42]
  };

  for(int i = 0; i < nbase; i++){
    basebeta[i] = log10(basebeta[i]);
    baselimit4pi[i] = log10(baselimit4pi[i]);
    baselimit1pi[i] = log10(baselimit1pi[i]);
  }

  TGraph base4pi(nbase, basebeta, baselimit4pi);
  TGraph base1pi(nbase, basebeta, baselimit1pi);

  const double crosstalkeff = (4519+4487+4596.)/3./6702.;

  const int ndelta = 9;
  double deltabeta[ndelta] = {
    1e-5,
    0.65, // smooth interpolation
    0.7,
    0.75,
    0.8,
    0.85,
    0.9,
    0.995,
    1
  };

  const double eff_for_no_delta_rays = 6702.;
  double deltaeff[ndelta] = {
    1,
    1,
    // Assume for the moment it is ok to base this on the 4pi study and
    // apply also to 1pi
    6306./eff_for_no_delta_rays,
    6104./eff_for_no_delta_rays,
    5795./eff_for_no_delta_rays,
    5311./eff_for_no_delta_rays,
    4456./eff_for_no_delta_rays,
    48./eff_for_no_delta_rays,
    1e-8
  };

  for(int i = 0; i < ndelta; i++){
    deltabeta[i] = log10(deltabeta[i]);
    if(deltaeff[i] > 0) deltaeff[i] = log10(deltaeff[i]);
  }

  TGraph * deltarayeff = new TGraph (ndelta, deltabeta, deltaeff);

  const double maxlogbeta = log10(0.995);

  double final_analysis_mult = 0;
  {
    const double livetime_2025_analysis = 234422296.38;

    const double additional_livetime =
      // A "few months" of as-yet unanalyzed data prior to 2025 analysis cutoff
      3 * 0.9 * 30 * 86400
      // From 2025-02-25 to 2028-12-31, optmistically assuming we get another
      // year of data because of LBNF delays and/or DARPA funding.
      // Assume 90% uptime for good data.
      + (10/12. + 3.0) * 365 * 86400 * 0.9;

     final_analysis_mult =
       (additional_livetime + livetime_2025_analysis)/livetime_2025_analysis;
  }

  // Use a large number of divisions across most of the plot, then handle
  // the last little bit separately.  Why?  Because I had accidentally
  // not handled beta = 0.95 to 0.99 and I want to keep the same points
  // that were working under 0.95.  Nothing fundamental.
  vector<double> logbetas;
  {
    const double minlogbeta = -3.75;
    const double thismaxlogbeta = -0.022;
    const int N = 400;
    for(int i = 0; i <= N; i++)
      logbetas.push_back(minlogbeta + i*(thismaxlogbeta-minlogbeta)/N);
  }
  {
    const int N = 80;
    const double thisminlogbeta = -0.022;
    for(int i = 1 /* sic: don't duplicate border point */; i <= N; i++)
      logbetas.push_back(thisminlogbeta + i*(maxlogbeta-thisminlogbeta)/N);
  }

  for(int i = 0; i < logbetas.size(); i++){
    const double logbeta = logbetas[i];
    const double beta = pow(10., logbeta);

    const double blimit4pi = pow(10, base4pi.Eval(logbeta));
    const double blimit1pi = pow(10, base1pi.Eval(logbeta));
    const double crosseff = beta > 0.68? crosstalkeff:
                            beta > 0.64? (beta-0.64)/(0.68-0.64) * (crosstalkeff-1) + 1
                            :1;
    const double deltaeff = pow(10, deltarayeff->Eval(logbeta));

    const double limit4pi = base4pi.Eval(logbeta) == 0? plottop:
      blimit4pi / crosseff / deltaeff;
    const double limit1pi = base1pi.Eval(logbeta) == 0? plottop:
      blimit1pi / crosseff / deltaeff;

    //#define FINAL
    #ifdef FINAL
      printf("%.4f %#10g %#10g\n", logbeta,
             limit1pi/final_analysis_mult, limit4pi/final_analysis_mult);
    #else
      printf("%.4f %#10g %#10g\n", logbeta,
             limit1pi, limit4pi);
    #endif
  }

  // The first element of baselimit4pi is a dummy value that makes the
  // fill area connect properly on the plot
  printf("%f %g %g\n", maxlogbeta, plottop, plottop);

  deltarayeff->Draw("ap*");
}
