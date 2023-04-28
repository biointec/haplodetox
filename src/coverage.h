/******************************************************************************
 *   Copyright (C) 2014 - 2023 Jan Fostier (jan.fostier@ugent.be)             *
 *   This file is part of HaploDetox                                          *
 *                                                                            *
 *   This program is free software; you can redistribute it and/or modify     *
 *   it under the terms of the GNU Affero General Public License as published *
 *   by the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                      *
 *                                                                            *
 *   This program is distributed in the hope that it will be useful,          *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *   GNU General Public License for more details.                             *
 *                                                                            *
 *   You should have received a copy of the GNU Affero General Public License *
 *   along with this program; if not, see <https://www.gnu.org/licenses/>.    *
 ******************************************************************************/

#ifndef COVERAGE_H
#define COVERAGE_H

#include <cmath>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include <numeric>
#include "util.h"
#include "arc.h"
#include "ssnode.h"

// ============================================================================
// CLASS PROTOTYPE
// ============================================================================

class Factor;
class DBGraph;

// ============================================================================
// MULTIPLICITY CLASS (Soft assignment)
// ============================================================================

class Multiplicity {

private:
//        std::vector<int> firstMult;     // strain multiplicities correponding to P[0]
        std::vector<int> card;          // number of possible values for each strain
        std::vector<double> P;          // probabilities of multiplicities (log-space)


        /**
         * Helper function to find a strain-multiplicity assignment based on an index in P
         * iterates over smaller indices in firstMult faster
         * e.g.
         * firstMult = [0,0]
         * card = [3,2]
         * 0 0
         * 1 0
         * 2 0
         * 0 1
         * 1 1
         * 2 1
         */
        std::vector<int> index2mult(int idx) const {
                assert(idx < P.size());

                std::vector<int> multiplier(card.size(), 1);
                for (size_t i = 1; i < card.size(); i++)
                        multiplier[i] = multiplier[i-1] * card[i-1];

                std::vector<int> assignment(card.size());
                for (size_t i = assignment.size(); i-- > 0; ) {
                        assignment[i] = idx / multiplier[i];
                        idx -= assignment[i] * multiplier[i];
                }

                return assignment;
        }

        /**
         * Helper function to find a strain-multiplicity assignment based on an index in P,
         * but store in passed along vector
         */
        void index2mult(int idx, std::vector<int>& assignment) const {
                assert(idx < P.size());

                std::vector<int> multiplier(card.size(), 1);
                for (size_t i = 1; i < card.size(); i++)
                        multiplier[i] = multiplier[i-1] * card[i-1];

                assignment.resize(card.size());
                for (size_t i = assignment.size(); i-- > 0; ) {
                        assignment[i] = idx / multiplier[i];
                        idx -= assignment[i] * multiplier[i];
                }
        }

public:
	static int maxM;	// maximal multiplicity, all multiplicities >= maxM get the same wildcard prb
        /**
         * Default constructor
         */
        Multiplicity() : card({1}) {
                P.push_back(0.0);       // log(1) = 0
        }

        /**
         * Constructor with single strain assignment
	 * All multiplicities except the target get zero probability
	 * Target gets probability 1.
         * @param multiplicity Target multiplicity
         */
        Multiplicity(const std::vector<int>& multiplicity) : card(std::vector<int>(multiplicity.size(), maxM+1)) {
		int size = std::accumulate(card.begin(), card.end(), 1, std::multiplies<int>());
		P = std::vector<double>(size,-INFINITY);
		int idx = 0;
		for (int i = 0; i < card.size(); i++){
			idx += std::min(multiplicity[i],maxM) * pow(maxM+1,i);
		}
		P[idx] = 0;
        }

        /**
         * Constructor from factor
         * @param card cardinality for the strain multiplicities
         * @param P Probabilities (log-space)
         */
        Multiplicity(const std::vector<int>& card, const std::vector<double>& P) :
                     card(card), P(P) {
                        int length = std::accumulate(card.begin(), card.end(), 1, std::multiplies<int>());
                        assert(P.size() == length); // make sure the length of P is correct
        }

        /**
         * Get the expected multiplicity (highest probability)
         * @return Expected multiplicity
         */
        std::vector<int> getExpMult() const;

        /**
         * Get the log(odds ratio) of the multiplicity being equal to the
         * expected multiplicity and the multiplicity being different from
         * the expected multiplicity. This provides a measure of confidence
         * in how certain we are that the assignment is correct
         * @return log( P(M = expMult) / P(M != expMult) )
         */
        double getExpMultLogOR() const;

        /**
         * Get the log probability of the multiplicity being equal to the
         * expected multiplicity
         * @return log( P(M = expMult) )
         */
        double getExpMultLProb() const;

        /**
         * Implicit conversion to integral type
         * @return Expected multiplicity
        operator int() const {
                return getExpMult();
        }*/


        /**
         * Get the probability corresponding to a particular multiplicity
         * @param mult Target multiplicity
         * @return Probability between [0..1]
         */
        double operator[](std::vector<int> mult) const;

        /**
         * Normalize the probabilities
         */
        void normalize();

        double getEntropy();


        /**
         * Operator << overloading
         * @param out Output stream
         * @param m Multiplicity to display
         * @return Output stream
         */
        friend std::ostream &operator<<(std::ostream &out, const Multiplicity &m);

        /**
         * Overload equality comparison operator (==)
         * @param m1 Multiplicity object left
         * @param m2 Multiplicity object right
         * @return equality based on expected multiplicity!
         */
        friend bool operator==(const Multiplicity &m1, const Multiplicity  &m2);

        /**
         * Overload inequality comparison operator (!=)
         * @param m1 Multiplicity object left
         * @param m2 Multiplicity object right
         * @return inequality based on expected multiplicity!
         */
        friend bool operator!=(const Multiplicity &m1, const Multiplicity  &m2);
};

// ============================================================================
// COVERAGE MODEL
// ============================================================================

class CovModel {

private:
	double errLambda;       // mean of the error component
        double errODF;          // error overdispersion factor

        double sosLambda;       // mean of sum-of-strains coverage
        double sosODF;          // overdispersion factor (also for strains)

        // relative fractions of the strains (sum to 1.0)
        std::vector<double> fractions;
        size_t maxMult;         // maximum multiplicity

        // log(weight) of the components (size = (maxMult + 1) ^ numStrains)
        std::vector<double> logw;

	double maxMultPVal = 0.0; // in log space TODO maybe check good value

        /**
         * A single iteration of EM to fit the coverage model to a histogram
         * @param data Histogram data as <x, y> pairs
         * @param xMin Minimum x value to use (lower x-values are infered)
         * @param xMax Maximum x value to use (higher x-values are ignored)
         * @param mu Packed mean vector
         * @param ODF Packed ODF vector
         * @param w Packed weight vector
         * @param nStrains Number of strains
         * @param maxEps Maximum relative tolerance for EM convergence
         * @return The relative difference (epsilon) between mu components
         */
        static double
        haploNBMixtureEMIt(const std::map<uint, double>& data, uint xMin,
                           uint xMax, std::vector<double>& mu,
                           std::vector<double>& ODF, std::vector<double>& w,
                           int nStrains, double maxEps, bool genusLevel);

        /**
         * Fit the coverage model to a histogram using expectation-maximization
         * @param data Histogram data as <x, y> pairs
         * @param xMin Minimum x value to use (lower x-values are infered)
         * @param xMax Maximum x value to use (higher x-values are ignored)
         * @param mu Packed mean vector
         * @param ODF Packed ODF vector
         * @param w Packed weight vector
         * @param nStrains Number of strains
         * @param maxEps Maximum relative tolerance for EM convergence
         * @param maxIteration Maximum number of iterations in EM
         * @return <numIterations, epsilon> tuple
         */
        static std::tuple<int, double>
        haploNBMixtureEM(const std::map<uint, double>& data, uint xMin,
                         uint xMax, std::vector<double>& mu,
                         std::vector<double>& ODF, std::vector<double>& w,
                         int nStrains, double maxEps,
                         int maxIterations, bool genusLevel);

public:
        /**
         * Constructor
         * @param errLambda Mean of the error term
         * @param errODF Error overdispersion factor
         * @param strainLambda Mean of each strain
         * @param strainODF Overdispersion factor of the strains
         * @param w Weight of the components
         * @param maxMult Maximum multiplicity
         */
        CovModel(double errLambda, double errODF,
                 const std::vector<double>& strainLambda, double strainODF,
                 const std::vector<double>& w, size_t maxMult);

        /**
         * Constructor from file
         * @param filename Input filename
         */
        CovModel(const std::string& filename);
        
        /**
         * Constructor of empty model
         */
        CovModel(): errLambda(1.0), errODF(1.0), sosLambda(100.0), sosODF(1.5){};

        /**
         * Fit the coverage model to a histogram using expectation-maximization
         * Only the error (00..0), sum-of-strain (11..1), strain (00..1..00)
         * and not-strain (11..0..11) components are learned from the histogram
         * @param histogram Histogram data as <x, y> pairs
         * @param xMin Minimum x value to use (lower x-values are infered)
         * @param xMax Maximum x value to use (higher x-values are ignored)
         * @param maxEps Maximum relative tolerance for EM convergence
         * @param maxIteration Maximum number of iterations in EM
         * @return <numIterations, epsilon> tuple
         */
        std::tuple<int, double> fitEM(const std::map<uint, double>& histogram,
                                      uint xMin, uint xMax, double maxEps,
                                      int maxIterations, bool genusLevel);

        /**
         * Get the error distribution mean
         * @return Error distribution mean
         */
        double getErrLambda() const {
                return errLambda;
        }

        /**
         * Set a new value for the error distribution mean
         * @param newErrLambda Target value for the error distribution mean
         */
        void setErrLambda(double newErrLambda) {
                errLambda = newErrLambda;
        }

        /**
         * Get the error distribution overdispersion factor
         * @return Error distribution overdispersion factor
         */
        double getErrorODF() const {
                return errODF;
        }

        /**
         * Set the error distribution overdispersion factor
         * @param newErrODF Target value for the error distribution ODF
         */
        void setErrorODF(double newErrODF) {
                errODF = newErrODF;
        }

        /**
         * Get the sum-of-strains distribution mean
         * @return The sum-of-strains distribution mean
         */
        double getSosLambda() const {
                return sosLambda;
        }

        /**
         * Set a new value for the sum-of-strains distribution mean
         * @param newSosLambda Target value for the sos distribution mean
         */
        void setSosLambda(double newSosLambda) {
                sosLambda = newSosLambda;
        }

        /**
         * Get the overdispersion factor for the strains and sos
         * @return The overdispersion factor for the strains and sos
         */
        double getODF() const {
                return sosODF;
        }

        /**
         * Set the overdispersion factor for the strains and sos
         * @param newODF The overdispersion factor for the strains and sos
         */
        void setODF(double newODF) {
                sosODF = newODF;
        }

        /**
         * Get the strain distribution mean
         * @param strain Strain index
         * @return The strain distribution mean
         */
        double getStrainLambda(size_t i) const {
                return fractions[i] * sosLambda;
        }

        /**
         * Get the number of strains
         * @return number of strains
         */
        size_t getNumStrains() const {
                return fractions.size();
        }

        /**
         * Get the maximum multiplicity
         * @return Maximum multiplicity
         */
        size_t getMaxMult() const {
                return maxMult;
        }

        /**
         * Get the logarithm of the weight for a given multiplicity
         * @param mult Multiplicity combination
         * @return The logarithm of the weight for the multiplicity
         */
        double getLWeight(const std::vector<int>& mult) const {
                size_t idx = 0;
                for (size_t i = 0, offset = 1; i < mult.size(); i++) {
                        idx += std::min<size_t>(mult[i], maxMult) * offset;
                        offset *= maxMult + 1;
                }
                return logw[idx];
        }

        /**
         * Get the weight for a given multiplicity
         * @param mult Multiplicity combination
         * @return The weight for the multiplicity
         */
        double getWeight(const std::vector<int>& mult) const {
                return exp(getLWeight(mult));
        }

        /**
         * Set the logarithm of the weight for a given multiplicity
         * @param mult Multiplicity combination
         * @param logWeight Target log(weight)
         */
        void setLWeight(const std::vector<int>& mult, double logWeight) {
                size_t idx = 0;
                for (size_t i = 0, offset = 1; i < mult.size(); i++) {
                        idx += std::min<size_t>(mult[i], maxMult) * offset;
                        offset *= maxMult + 1;
                }
                logw[idx] = logWeight;
        }

        /**
         * Set the weight for a given multiplicity
         * @param mult Multiplicity combination
         * @param weight Target weight
         */
        void setWeight(const std::vector<int>& mult, double weight) {
                setLWeight(mult, log(weight));
        }

        // ================

        double getMaxMultP() const {
		return maxMultPVal;
	}

	void setMaxMultP(double newP) {
		maxMultPVal = newP;
	}



        /**
         * Get the vector with the fractions of the strains in the dataset
         * @return the fractions vectior
         */
        std::vector<double> getFractions() const {
                return fractions;
        }

        /**
         * Get the fraction of the strain at position strainID
         * @param strainID the position in the fractions vector
         * @return the appropriate fraction
         */
        double getFraction(int strainID) const {
                return fractions[strainID];
        }

        /**
         * Set the weights (uniform accross multiplicities)
         * @param w Desired weight
         */
        void setWeightUniform(double w) {
                for (auto& it : logw)
                        it = log(w);
        }
        
        void setZeroWeight(double w) {
                logw[0] = log(w);
        }

        /**
         * Get the number of components in the model
         * @return The number of components
         */
        size_t getC() const {
                return maxMult + 1;
        }

        /**
         * Get the size of w
         * @return size of w
         */
        size_t getK() const {
                return logw.size();
        }

        /**
         * Get the hard assignment given the observed coverage
         * TODO adapt --> for now it gets a multiplicity assignment (x,x,x) based on coverage vs. lambda
         * @param obsCov Observed coverage
         * @return The most likely multiplicity
         */
        std::vector<int> getExpMult(double obsCov) const;

        /**
         * Get the soft assignment given the observed coverage (unnormalised)
         * @param obsCov Observed coverage
         * @return Multiplicity (soft assignment)
         */
        Multiplicity getMultSoft(double obsCov) const;

        /**
         * Get the logprob of observing obsCov given the multiplicity
         * Note: NOT normalised: the value might be bigger than 1 depending on w
         * @param obsCov Observed coverage
         * @param mult Multiplicity (0 = error term)
         * @return log[ P(obsCov | multiplicity) ]
         */
        double getLogProb(double obsCov, const std::vector<int>& mult) const;

        /**
         * Get the probability of observing obsCov given the multiplicity
         * Note: NOT normalised: the value might be bigger than 1 depending on w
         * @param obsCov Observed coverage
         * @param mult Multiplicity (0 = error term)
         * @return P(obsCov | multiplicity)
         */
        double getProb(double obsCov, const std::vector<int>& mult) const {
                return exp(getLogProb(obsCov, mult));
        }

        /**
         * Operator << overloading
         * @param out Output stream
         * @param c Coverage model to display
         * @return Output stream
         */
        friend std::ostream &operator<<(std::ostream &out, const CovModel &c);

        /**
         * Write out parameters of coverage model
         * @param covFilename name of file to write to
         */
        void write(const std::string& covFilename) const;

        /**
         * Read parameters of coverage model from disk
         * @param covFilename name of file to read from
         */
        void read(const std::string& covFilename);

        /**
         * Write histogram and model fit to GNUplot file
         *
         * @param baseFilename filename
         * @param hist Coverage data
         */
        void writeGnuplot(const std::string& baseFilename,
                          const std::map<unsigned int, double>& hist,
                          int maxPlotMult = 1) const;

        /**
         * Write histogram and model fit of initialization to GNUplot file
         *
         * @param baseFilename filename
         * @param hist Coverage data
         * @param initW Weights found for initialization
         * @param lambdaHighM extra lambda for mults >1
         * @param ODFhighM ODF for mults >1
         */
        void writeGnuplotInit(const std::string& baseFilename,
                          const std::map<unsigned int, double>& hist,
                          const std::vector<double>& initW,
                          const double& lambdaHighM,
                          const double& ODFhighM) const;

        /**
         * Get the coverage below which an error is most likely
         * TODO Unsure if this will still be effective,
         * as the separate strains might have very low weights compared to w[0]
         *
         * @param epsilon Accuracy of the estimate
         * @return Coverage cutoff value
         */
        double getCovCutOff(double OR, double epsilon = 1e-3) const
        {
                // use bisection algorithm to find good estimate
                double left = errLambda;

                int minIdx = std::min_element(fractions.begin(),fractions.end()) - fractions.begin();
                // Or should we test all possible strainmult = 1 assignments?
                double right = sosLambda * fractions[minIdx];
                std::vector<int> err(fractions.size(),0);
                std::vector<int> minfrac(fractions.size(),0);
                minfrac[minIdx] = 1;

                do {
                        double mid = 0.5 * (left + right);
                        double Pe = getLogProb(mid, err);
                        double Pt = getLogProb(mid, minfrac);

                        if (Pe > OR * Pt)
                                left = mid;
                        else
                                right = mid;
                } while (right - left > epsilon);

                return 0.5 * (right + left);
        }

        /**
         * Compute number of (independent) parameters in the model
         * e.g. to use in a scoring function
         * @return the number of (independent) parameters
         */
        int getNumParams() const
        {
                //err-param +lambda_i + odf + weights
                int nstrains = fractions.size();
                int nparam = 2 + nstrains + 1;
                if( maxMult == 2) {     /// all weights for m_i combinations with same number 0,1,2 mults are equalised
                        int indep_w = 0; 
                        for(int i = 0; i <= nstrains; i++){ // number of 0's chosen
                                indep_w += nstrains - i; // number of 1's chosen
                        }
                        nparam += indep_w;
                } else {                /// only the weights corresponding with only 0's and x's (x = 1,2,...,maxMult) are equalised
                        nparam += logw.size();
                        int w_dependency = 0;
                        for(int i = 1; i < nstrains; i++){
                                w_dependency += (Util::factRatio(nstrains, std::max(nstrains - i,i)) / Util::factorial(std::min(nstrains-i,i)) - 1);
                        }
                        w_dependency *= maxMult;
                        nparam -= w_dependency;
                }
                return nparam;
        }

        /**
         * Compare the CovModel's estimation to the histogram of observations
         * This function assumes the weights were computed as (weighted) counts
         * of the number of observations assigned a certain multiplicity.
         * If the weights are normalised, this function should be adapted to normalise the histogram first!
         * @param hist observed coverage counts of (a subset of) nodes/arcs
         * @return the MSE = (sum_n (y_obs - y_est)^2)/(n - n_param)
         */
        double getMSE(const std::map<unsigned int, double>& hist) const;

        /**
         * Calculate the log likelihood based on this CovModel and all observations
         * Note: While the parameters in CovModel might have been estimated in an EM
         *       that used the CRF, this LL does not take the CRF into account, only
         *       the resulting mixture model
         * @param obsCov: list of the observed coverages
         * @return logL = sum_n log(sum_c (p_c P(cov|param_c))
         */
        double getLL(const std::vector<double>& obsCov) const;

        
        double getCompleteLL(NodeMap<Multiplicity> mult, const DBGraph& dBG);
        double getCompleteLL(EdgeMap<Multiplicity> mult, const DBGraph& dBG);
        
        double getAIC(double ll) const
        {
                double aic = -2*ll;
                // gewichten + lambdas (gezien som, moeten we 1 niet meetellen) + error params + odf
                int nparam = getNumParams(); // TODO subtract dependencies in w
                aic += 2*nparam;
                return aic;
        }

        double getBIC(double ll, int nObs) const
        {
                double bic = -2*ll;
                // gewichten + lambdas (gezien som, moeten we 1 niet meetellen) + error params + odf
                int nparam = getNumParams(); // TODO subtract dependencies in w
                bic += log(nObs) * nparam;
                return bic;
        }

        /**
         * Calculate the entropy of the multiplicity assignments (tau) of this CovModel
         * @param obsCov : observed coverages
         * @return EN(tau) = sum_n(sum_c(tau_nc log(tau_nc)))
         */
        double getEntropy(const std::vector<double>& obsCov) const;

        /**
         * Calculate the entropy in the multiplicity assignments made using the CRF
         * @param mult: the Multiplicity assignments tau made with the CRF
         * @return EN(tau) = sum_n(sum_c(tau_nc log(tau_nc)))
         */
        double getCRFassignmentEntropy(NodeMap<Multiplicity> mult);
        
        double getCRFassignmentEntropy(EdgeMap<Multiplicity> mult);

        double getCLC(double ll, double entropy) const
        {
                double clc = -2 * ll + 2 * entropy;
                return clc;

        }

        /*
         * Biernacki, Celeux & Govaert 2000 (see also Usami2014)
         */
        double getICL(double ll, double entropy, double nObs) const
        {
                return getBIC(ll, nObs) + 2 * entropy;
        }

        /**
         * Criterion defined by Ramaswamy, DeSarbo, Reibstein and Robinson (1993)
         * Larger is better here!
         * @param entropy the entropy of the predictions
         * @param nObs the number of observations or predictions
         * @return entropy-based information criterion
         */
        double entropyCriterion(double entropy, double nObs) const
        {
                return 1.0 - (entropy/(nObs * log(logw.size())));
        }

};

#endif
