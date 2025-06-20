void make_sens_slowonly()
{
  const double x = 7.65 + 7.58, y = 7.65 + 7.49, z = 59.62;

  const double Aproj = (x*y + y*z + z*x)/2 * 100 * 100;

  // Limit for 100% efficiency for heavy monopoles
  const double livetime = 236.977204e6;
  const double baselimit = log(10)/ ( Aproj * livetime * 4*M_PI);

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
    -2.25,
    -2.2,
    -2.15,
    -2.1,
    -2.05,
    -2.0,

    -1.99,
    -1.98,
    -1.97
  };

  double eff[n] = {
    0.001,
    0.053,
    0.17,
    0.321,
    0.418,
    0.446,
    0.462,
    0.467,
    0.471,
    0.473,
    0.475,
    0.479,
    0.473,
    0.469,

    0.414,
    0.24585361,
    0.146,
    0.092021737,
    0.058,
    0.025258662,
    0.011,

    0.005,
    0.0005,
    0.00005,
  };


  // Product of the ratio of projected area for the top 1pi solid angle
  // to the naive quarter and the efficiency in that top 1pi 
  double corr_1pi[n] = {
     3.1311, //    -3.7   
     1.0726, //    -3.6   
     1.0878, //    -3.5   
      1.041, //    -3.4   
     1.0585, //    -3.3   
     1.0025, //    -3.2   
     1.0167, //    -3.1   
    0.98575, //      -3   
     1.0314, //    -2.9   
    0.98624, //    -2.8   
     1.0345, //    -2.7   
     1.0314, //    -2.6   
     1.0081, //    -2.5   
     1.0041, //    -2.4   

    0.97383, //    -2.3   
    0.52733, // XXX -2.25 -2.2218   
    0.51923, //    -2.2   
    0.50181, // XXX -2.15 -2.1549   
    0.50038, //    -2.1   
    0.50874, // XXX -2.05 -2.0457   
    0.30541, //      -2   
    // I don't think the above numbers can be right, because they imply negative efficiency for near-horizon monopoles

    0.3, // -1.99
    0.3, // -1.98
    0.3  // -1.97
  };

  double logeff4pi[n] = { 0 }, logeff1pi[n] = { 0 }, eff1pi[n] = { 0 };

  const double highmass_solidangle = 4*M_PI;
  const double lowmass_solidangle = 1*M_PI;

  for(int i = 0; i < n; i++) logeff4pi[i] = log10(eff[i]);
  for(int i = 0; i < n; i++) eff1pi[i] = eff[i] / corr_1pi[i] * lowmass_solidangle/highmass_solidangle;
  for(int i = 0; i < n; i++) logeff1pi[i] = log10(eff1pi[i]);

  TSpline3 spline1pi("whatever", logbeta, eff1pi, n);
  TGraph graph1pi(n, logbeta, eff1pi);

  TSpline3 spline4pi("whatever", logbeta, eff, n);
  TGraph graph4pi(n, logbeta, eff);

  const double minlogbeta = logbeta[0];
  const double maxlogbeta = logbeta[n-1];

  const int npoints = 200;

  printf("%f 1.125e-11 1.125e-11\n", logbeta[0]);

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
    const double the_eff1pi = spline1pi.Eval(lbeta);
    const double the_effg1pi = graph1pi.Eval(lbeta);

    const double the_eff4pi = spline4pi.Eval(lbeta);
    const double the_effg4pi = graph4pi.Eval(lbeta);

    const double limit1pi = baselimit / (lbeta < -2.5? the_eff1pi: the_effg1pi);

    const double limit4pi = baselimit / (lbeta < -2.5? the_eff4pi: the_effg4pi);

    #define FINAL
    #ifdef FINAL
      printf("%f %g %g\n", lbeta,
             limit1pi/final_analysis_mult, limit4pi/final_analysis_mult);
    #else
      printf("%f %g %g\n", lbeta,
             limit1pi, limit4pi);
    #endif
  }

  printf("%f 1.125e-11 1.125e-11\n", logbeta[n-1]);
}
