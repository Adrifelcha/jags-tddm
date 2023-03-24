/*  This class contains functions for the tddm (Thurstonian drift-diffusion model) */
#ifndef DTDDM_H_
#define DTDDM_H_

#include <distribution/VectorDist.h>

namespace jags {

	namespace tddm {

	/**
	 * Implements the Thurstonian drift-diffusion model
	 *
	 * <pre>
         * X[i,1] <- choice (in radians)
         * X[i,2] <- response time (in seconds)
         * X[i,1:2] ~ dtddm(xdrift, ydrift, bound, tzero)
	 * </pre>
         *
	 * @short Thurstonian drift diffusion model
	 */

		class DTDDM : public VectorDist
		{
		public:
			DTDDM();

			unsigned int length(std::vector<unsigned int> const &len) const;

			double logDensity(double const *x, unsigned int length,
					PDFType type,
					std::vector<double const *> const &par,
					std::vector<unsigned int> const &len,
					double const *lower, double const *upper) const;

			void randomSample(double *x, unsigned int length,
					std::vector<double const *> const &par,
					std::vector<unsigned int> const &len,
					double const *lower, double const *upper,
					RNG *rng) const;

			void typicalValue(double *x, unsigned int length,
					std::vector<double const *> const &par,
					std::vector<unsigned int> const &len,
					double const *lower, double const *upper) const;

			bool checkParameterValue(std::vector<double const *> const &par,
					std::vector<unsigned int> const &len) const;

			bool checkParameterLength(std::vector<unsigned int> const &len) const;

			void support(double *lower, double *upper, unsigned int length,
					std::vector<double const *> const &par,
					std::vector<unsigned int> const &len) const;

			bool isSupportFixed(std::vector<bool> const &fixmask) const;

			unsigned int df(std::vector<unsigned int> const &len) const;

		};
	} //namespace tddm
} //namespace jags

#endif /* DTDDM_H_ */
