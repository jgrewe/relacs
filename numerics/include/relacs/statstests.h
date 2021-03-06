/*
  statstests.h
  

  RELACS - Relaxed ELectrophysiological data Acquisition, Control, and Stimulation
  Copyright (C) 2002-2015 Jan Benda <jan.benda@uni-tuebingen.de>

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  RELACS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _RELACS_STATSTESTS_H_
#define _RELACS_STATSTESTS_H_

#include <relacs/array.h>
#include <relacs/sampledata.h>

namespace relacs {


  /*! The integral over the standard normal distribution upto x. 
      \f[ \alpha = \int_{-\infty}^x \frac{1}{\sqrt{2\pi}}e^{-z^2/2} dz \f]
  */
double alphaNormal( double x );


  /*! The integral over the binominal distribution upto inklusively \a k.
      \f[ \alpha = \sum_{j=0}^k {n \choose j} p^j (1-p)^{n-j} \f]
      The following symmetry holds: 
      alphaBinomial( k, n, p ) = 1.0-alphaBinomial( n-k-1, n, 1-p ).
      For \f$ p \approx 0.5 \f$ and large \a n, alphaBinomial() 
      can be approximated by a normal distribution with mean \f$ np \f$
      and variance \f$ np(1-p) \f$. */
double alphaBinomial( int k, int n, double p );


  /*! Return the number of elements in \a data
      with value greater than \a median. */
int positiveSign( const ArrayD &data, double median );
  /*! The sign test for the median.
      Tests whether the data have a median less than (\a tail < 0).
      equal to (\a tail = 0), or greater than (\a tail > 0) \a median.
      Returns in \a p the significane level. */
void signTest( const ArrayD &data, double median,
	       int tail, int &n, double &p );

  /*! Returns the rank sum of the positive differences \a ydata - \a xdata
      computed from the matched pairs in \a xdata and \a ydata. 
      In \a n the number of samples contributing to Wilcoxon's W is returned. */
double rankSumWilcoxon( const ArrayD &xdata, const ArrayD &ydata, int &n );
  /*! The one-tailed significance for an observed sum of ranks
      \a w and sample size \a n for the Wilcoxon test.
      \f[ \alpha = \int_0^w p(w',n) dw' \f]
      Use this only for \a n < 20, since the execution time is n*2^n!
      For larger \a n, \a z returned by zWilcoxon( w, n )
      is standard normal distributed and
      alphaNormal( zWilcoxon( w, n ) ) is therefore a good approximation
      for alphaWilcoxon(). 
      Note, however, that for a perfect separation of the pairs (\a w = 0)
      the significance is given by \f$ \alpha_{min} = 2^{-n} \f$, i.e.
      a sample size of at least 5 is needed for an 0.05 significance level.
      The following symmetry holds: 
      alphaWilcoxon( w, n ) = 1.0-alphaWilcoxon( n*(n+1)/2-w-1, n ) */
double alphaWilcoxon( double w, int n );
  /*! For a given \a w and sample size \a n from a Wilcoxon test
      returns the corresponding z 
      \f[ z = \frac{W - n(n+1)/4}{\sqrt{n(n+1)*(2*n+1)/24}} \f]
      which is standard normal distributed. */
double zWilcoxon( double w, int n );
  /*! Perform Wilcoxon test on the matched pairs in \a xdata and \a ydata.
      If \a tail > 0 tests whether ydata > xdata (one tailed ),
      if \a tail < 0 tests whether ydata < xdata (one tailed ),
      If \a tail = 0 tests whether ydata = xdata (two tailed ).
      Returns in \a w the rank sum of the specified sign
      and in \a p the significance level. */
void wilcoxonTest( const ArrayD &xdata, const ArrayD &ydata, int tail,
		   double &w, double &p );

  /*! Returns the significance level of Pearson's correlation coefficient \a r
      optained from \a n pairs of data.
      If n<=2, then 1.0 is returned. */
double pearsonTest( double r, int n );


  /*! The Kolmogorov-Smirnov test for comparing a set of data values with
      a theoretically known distribution.
      \note \a data must be a sorted array of data values.
      \param[in] data the observed data values (not their distribution or cumulative!).
      \param[in] density the computed probability distribution 
      (does not need to be normalized) to which \a data are compared
      (see SampleData::cumulative()).
      From this function the cumulative is computed with linear interpolation.
      \param[out] d K-S statistics D
      \param[out] p significance level of the disproof of the null hypothesis
      that the distributions are the same. */
void KSTest( const ArrayD &data, const SampleDataD &density, double &d, double &p );

  /*! The runs test (or Wald–Wolfowitz test) checks a randomness hypothesis 
      for a two-valued data sequence.
      \param[in] data the data from which the runs are determined.
      Each run is a series of consecutive positive or negative data values.
      \param[out] z the Z statistics of the runs test. 
      \param[out] p the p-value for the runs being random. */
void runsTest( const ArrayD &data, double &z, double &p );

  /*! Compute the complementary normalized incomplete Gamma Function 
      \f[ P(a,x) = 1/\Gamma(a) \int_0^x t^{a-1} \exp(-t) dt \f] for a > 0, x >= 0.
      If a = 0, 1-exp(-x) is returned. */
double gammaP( double a, double x );
  /*! Compute the normalized incomplete Gamma Function 
      \f[ Q(a,x) = 1/\Gamma(a) \int_x^\infty t^{a-1} \exp(-t) dt \f] for a > 0, x >= 0
      If a = 0, exp(-x) is returned. */
double gammaQ( double a, double x );
  /*! Compute the normalized incomplete beta function 
      \f[ B_x(a,b)/B(a,b) \f] for \f$ a > 0\f$, \f$b > 0\f$, and \f$ 0 \le x \le 1\f$,
      where
      \f[ B(a,b) = \int_0^1 t^{a-1} (1-t)^{b-1} dt = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)} \f]
      is the beta function and
      \f[ B(x;a,b) = \int_0^x t^{a-1} (1-t)^{b-1} dt \f]
      is the incomplete beta function.
      If a = 0, \f$ 1-(1-x)^b\f$ is returned. */
double incBeta( double a, double b, double x );


}; /* namespace relacs */

#endif /* ! _RELACS_STATSTESTS_H_ */
