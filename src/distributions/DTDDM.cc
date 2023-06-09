/* ===================================
    Thurstonian drift-diffusion model
   ===================================
*/

#include <config.h>
#include "DTDDM.h"

#include <util/dim.h>
#include <rng/RNG.h>
#include <util/nainf.h>
#include <cmath>
#include <JRmath.h>
#include <stdio.h>

using std::vector;
using std::log;
using std::string;

static inline double XDRIFT (vector<double const*> const &par) { return *par[0]; }
static inline double YDRIFT (vector<double const*> const &par) { return *par[1]; }
static inline double BOUND  (vector<double const*> const &par) { return *par[2]; }
static inline double TZERO  (vector<double const*> const &par) { return *par[3]; }

const double inv2pi = 0.159154943091895;
const double log2pi = 1.837877066409345;

const double j0_over_J1_of_j0[] =
	{
		4.632259e+00, -1.622289e+01,  3.187937e+01, -5.072504e+01,  7.228843e+01,
	       -9.626154e+01,  1.224225e+02, -1.506013e+02,  1.806628e+02, -2.124954e+02,
		2.460056e+02, -2.811133e+02,  3.177488e+02, -3.558508e+02,  3.953649e+02,
	       -4.362423e+02,  4.784390e+02, -5.219150e+02,  5.666336e+02, -6.125612e+02,
		6.596669e+02, -7.079218e+02,  7.572992e+02, -8.077741e+02,  8.593232e+02,
	       -9.119245e+02,  9.655574e+02, -1.020202e+03,  1.075841e+03, -1.132456e+03,
		1.190031e+03, -1.248550e+03,  1.307997e+03, -1.368360e+03,  1.429623e+03,
	       -1.491774e+03,  1.554801e+03, -1.618691e+03,  1.683434e+03, -1.749017e+03,
		1.815430e+03, -1.882663e+03,  1.950707e+03, -2.019550e+03,  2.089186e+03,
	       -2.159604e+03,  2.230795e+03, -2.302752e+03,  2.375467e+03, -2.448931e+03,
	};

const double j0_squared[] =
	{
	    5.783185962947,	   30.471262343662,	   74.887006790695,	  139.040284426460,	  222.932303617634,
	  326.563352932328,	  449.933528518036,	  593.042869655955,	  755.891394783933,	  938.479113475694,
	 1140.806031099645,	 1362.872150854105,	 1604.677474740231,	 1866.222004061853,	 2147.505739697844,
	 2448.528682258052,	 2769.290832176359,	 3109.792189768249,	 3470.032755267558,	 3850.012528850570,
	 4249.731510652210,	 4669.189700777160,	 5108.387099307649,	 5567.323706308981,	 6045.999521833563,
	 6544.414545923833,	 7062.568778614458,	 7600.462219933975,	 8158.094869906055,	 8735.466728550457,
	 9332.577795883786,	 9949.428071920078,	10586.017556671255,	11242.346250147502,	11918.414152357564,
	12614.221263308977,	13329.767583008270,	14065.053111461117,	14820.077848672467,	15594.841794646662,
	16389.344949387512,	17203.587312898377,	18037.568885182227,	18891.289666241712,	19764.749656079162,
	20657.948854696686,	21570.887262096137,	22503.564878279212,	23455.981703247398,	24428.137737002071
	};

#define smax 50  /* Maximum number of steps in infinite sum */

#define DEBUG FALSE

namespace jags {
	namespace tddm {

		DTDDM::DTDDM()
			: VectorDist("dtddm", 4)
		{}

		unsigned int DTDDM::length(vector<unsigned int> const &len) const
		{
			return 2;
		}

		bool DTDDM::checkParameterLength(vector<unsigned int> const &len) const
		{
			if (DEBUG) printf("checkParameterLength() has been called\n");

			if (len[0] != 1)  return false;
			if (len[1] != 1)  return false;
			if (len[2] != 1)  return false;
			if (len[3] != 1)  return false;
			return true;
		}

		bool DTDDM::checkParameterValue(vector<double const *> const &par,
				vector<unsigned int> const &len) const
		{
			if (DEBUG) printf("checkParameterValue() has been called\n");

			double bound = BOUND(par);
			double tzero = TZERO(par);

			if (tzero < 0)  return false;
			if (bound < 0)  return false;

			return true;
		}

		double DTDDM::logDensity(double const *x, unsigned int length,
				PDFType type,
				vector<double const *> const &par,
				vector<unsigned int> const &len,
				double const *lower, double const *upper) const
		{
			double c = x[0];
			double t = x[1];

			double bound = BOUND(par);
			double tzero = TZERO(par);
			double mu1 = XDRIFT(par);
			double mu2 = YDRIFT(par);
			double inva2 = 1.0 / (bound*bound);

			if (DEBUG) printf("c: %f | t: %f\n", c, t);
			if (DEBUG) printf("xdrift: %f | ydrift: %f | bound: %f | tzero %f\n", mu1, mu2, bound, tzero);

			double exponand, sum = 0.0, logPDF;

			for (int i=0; i<smax; i++) {
				exponand = j0_squared[i] * (t-tzero) * inva2 * -0.5;
				sum += exp(exponand) * j0_over_J1_of_j0[i];
				if (DEBUG) printf("exponand = %f\n", exponand);
				if (DEBUG) printf("sum = %f\n", sum);
			}

			logPDF = log(sum) + log(inva2);
			logPDF += bound*(mu1*cos(c)+mu2*sin(c));
			logPDF -= ((mu1*mu1 + mu2*mu2)*(t-tzero))*0.5;
			if (DEBUG) printf("logPDF = %f\n", logPDF);

			return isnan(logPDF) ? JAGS_NEGINF : logPDF;
		}

		void DTDDM::randomSample(double *x, unsigned int length,
				vector<double const *> const &par,
				vector<unsigned int> const &len,
				double const *lower, double const *upper,
				RNG *rng) const
		{
                        x[0] = 0.0;
                        x[1] = 0.5;
			return;
		}

		void DTDDM::support(double *lower, double *upper, unsigned int length,
			   vector<double const *> const &par,
			   vector<unsigned int> const &len) const
		{
			lower[0] = JAGS_NEGINF;
			upper[0] = JAGS_POSINF;
			lower[1] = 0.0;
			upper[1] = JAGS_POSINF;
		}

		void DTDDM::typicalValue(double *x, unsigned int length,
				vector<double const *> const &par,
				vector<unsigned int> const &len,
				double const *lower, double const *upper) const
		{
			x[0] = 0.0;
			x[1] = 0.5;
		}

		bool DTDDM::isSupportFixed(vector<bool> const &fixmask) const
		{
			return true;
		}

		unsigned int DTDDM::df(vector<unsigned int> const &len) const
		{
			return 1;
		}

	} //namespace tddm

} //namespace jags
