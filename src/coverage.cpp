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

#include <functional>
#include <numeric>

#include "dbgraph.h"
#include "kmernpp.h"
#include "settings.h"
#include "coverage.h"
#include "pgm/factor.h"

using namespace std;


// ============================================================================
// OVERLOAD vector OS
// ============================================================================

std::ostream& operator<<(std::ostream& os, const std::vector<int> &input)
{
	for (auto const& i: input) {
		os << " " << i;
	}
	return os;
}

// ============================================================================
// MULTIPLICITY CLASS
// ============================================================================

int Multiplicity::maxM = 3;

vector<int> Multiplicity::getExpMult() const
{
        // sanity check
        assert(!P.empty());

        // return the multiplicity with the highest probability
        int maxOffset = 0; double maxProb = P[0];
        for (int i = 1; i < (int)P.size(); i++) {
                if (P[i] > maxProb) {
                        maxProb = P[i];
                        maxOffset = i;
                }
        }

        vector<int> offset2mult = index2mult(maxOffset);

        return offset2mult;
}

double Multiplicity::getExpMultLogOR() const
{
        // if we only have one option, the odds ratio is infinite
        if (P.size() <= 1)
                return numeric_limits<double>::max();

        // get the most likely option
        int maxOffset = 0; double maxProb = P[0];
        for (int i = 1; i < (int)P.size(); i++) {
                if (P[i] > maxProb) {
                        maxProb = P[i];
                        maxOffset = i;
                }
        }

        // compute the sum of the probabilities of the other options (log-space)
        double logSum = numeric_limits<double>::lowest();       // log(0) = -oo
        for (int i = 0; i < (int)P.size(); i++) {
                if (i == maxOffset)
                        continue;
                logSum = Util::log_sum_exp(logSum, P[i]);
        }

        return maxProb - logSum;
}

double Multiplicity::getExpMultLProb() const
{
        // if we only have one option, the log of the probability is zero
        if (P.size() <= 1)
                return 0.0;

        // return the multiplicity with the highest probability
        int maxOffset = 0; double maxProb = P[0];
        for (int i = 1; i < (int)P.size(); i++) {
                if (P[i] > maxProb) {
                        maxProb = P[i];
                        maxOffset = i;
                }
        }

        return maxProb;
}

double Multiplicity::operator[](vector<int> mult) const
{
        assert(mult.size() == card.size());
        int strain;
        int idx = 0;
	int cardProd = 1;
        for (strain=0 ; strain<mult.size(); strain++){
            if (mult[strain]<0){
                return 0.0; //smaller than firstMult
            }
            if (mult[strain] >= card[strain]){
                return 0.0; //bigger than maxMult
            }

            idx += mult[strain] * cardProd;
	    cardProd *= card[strain];
        }
        double retValue = exp(P[idx]);
        if (retValue < 0.0 | retValue > 1.0){
            return 0.0;
        }
        return retValue;
}

void Multiplicity::normalize()
{
        if (P.empty())
                return;

        double logSum = P.front();
        for (size_t i = 1; i < P.size(); i++) {
                double m = max<double>(logSum, P[i]);
                logSum = m + log( exp(P[i]-m) + exp(logSum-m) );
        }

        for (size_t i = 0; i < P.size(); i++)
                P[i] = P[i] - logSum;
}

double Multiplicity::getEntropy()
{
        double entropy = 0.0;
        normalize();
        for (int i = 0; i < P.size(); i++){
                double tau = P[i];
                entropy += (exp(tau)*tau);
        }
        return -entropy;
}

std::ostream &operator<<(std::ostream &out, const Multiplicity &m)
{
        vector<int> multContainer(m.card.size());
        for (int i = 0; i < (int)m.P.size(); i++) {
                m.index2mult(i,multContainer);
                out << "( ";
                for (int mult : multContainer)
                        out << mult << " ";
                out << ")\t" << m.P[i] << "\n";
        }

        return out;
}

bool operator==(const Multiplicity &m1, const Multiplicity  &m2)
{
        if (m1.card.size() != m2.card.size())
                return false;

        vector<int> e1 = m1.getExpMult();
        vector<int> e2 = m2.getExpMult();

        return (e1 == e2);
}

bool operator!=(const Multiplicity &m1, const Multiplicity  &m2)
{
        return !(m1 == m2);
}

// ============================================================================
// COVERAGE MODEL
// ============================================================================

CovModel::CovModel(double errLambda, double errODF,
                   const std::vector<double>& strainLambda, double strainODF,
                   const std::vector<double>& w, size_t maxMult) :
        errLambda(errLambda), errODF(errODF),
        sosODF(strainODF), maxMult(maxMult)
{
        assert(maxMult >= 2);    // we want at least mult = 1 and mult >= 2
        // make sure we have a correctly sized weight vector
        assert(round(pow(maxMult + 1, strainLambda.size())) == w.size());

        // compute sum-of-strains lambda and the strain fractions
        sosLambda = accumulate(strainLambda.begin(), strainLambda.end(), 0.0);
        for (auto sl : strainLambda)
                fractions.push_back(sl / sosLambda);

        // store the logarithm of the weights
        logw.resize(w.size());
        const double pseudoCount = 1.0;
        for (size_t i = 0; i < logw.size(); i++)
                logw[i] = log(max<double>(w[i], pseudoCount));
}

tuple<int, double> CovModel::fitEM(const std::map<uint, double>& histogram,
                                   uint xMin, uint xMax, double maxEps,
                                   int maxIterations, bool genusLevel)
{
        const size_t& numStrains = fractions.size();

        // create a "packed" lambda vector that contains: a) error lambda,
        // b) strain lambdas, c) non-strain lambdas, d) sos-lambdas
        vector<double> pLambda = { errLambda };
        for (size_t i = 0; i < numStrains; i++)
                pLambda.push_back(fractions[i] * sosLambda);
        if (numStrains > 2)
                for (size_t i = 0; i < numStrains; i++)
                        pLambda.push_back((1.0 - fractions[i]) * sosLambda);
        for (size_t i = 1; i <= maxMult; i++)
                pLambda.push_back(i * sosLambda);

        // create a "packed" ODF vector that contains error ODF + sos ODF
        vector<double> pODF(pLambda.size(), sosODF);
        pODF[0] = errODF;

        // create a "packed" weight vector
        vector<double> pWeight = { exp(logw[0]) };
        for (size_t i = 0; i < numStrains; i++) {
                vector<int> a(numStrains, 0); a[i] = 1;
                pWeight.push_back( getWeight(a) );
        }
        if (numStrains > 2)
                for (size_t i = 0; i < fractions.size(); i++) {
                        vector<int> a(numStrains, 1); a[i] = 0;
                        pWeight.push_back( getWeight(a) );
                }
        for (size_t i = 1; i <= maxMult; i++) {
                vector<int> a(numStrains, i);
                pWeight.push_back( getWeight(a) );
        }

        // fit the spectrum using Expectation - Maximization
        tuple<int, double> result = haploNBMixtureEM(histogram, xMin, xMax,
                pLambda, pODF, pWeight, fractions.size(), maxEps, maxIterations, genusLevel);

        // unpack the pLambda vector
        errLambda = pLambda[0];
        sosLambda = accumulate(pLambda.begin() + 1,
                               pLambda.begin() + 1 + numStrains, 0.0);
        for (size_t i = 0; i < numStrains; i++)
                fractions[i] = pLambda[i+1] / sosLambda;

        // fill in the error weight
        size_t counter = 0;
        logw[0] = log(pWeight[counter++]);

        // fill in the strain weights
        Assignment ass(vector<int>(numStrains, maxMult + 1));
        for (int i = 0; i < numStrains; i++) {
                vector<int> a(numStrains, 0); a[i] = 1;
                size_t index = ass.assignmentToIndex(a);
                logw[index] = log(pWeight[counter++]);
        }

        // fill in the not-strain weights
        for (int i = 0; (i < numStrains) && (numStrains > 2); i++) {
                vector<int> a(numStrains, 1); a[i] = 0;
                size_t index = ass.assignmentToIndex(a);
                logw[index] = log(pWeight[counter++]);
        }

        // fill in sum components
        for (int i = 1; i <= maxMult; i++) {
                vector<int> a(numStrains, i);
                size_t index = ass.assignmentToIndex(a);
                logw[index] = log(pWeight[counter++]);
        }

        // set the ODF components
        errODF = pODF[0];
        sosODF = pODF[1];

        return result;
}

double CovModel::haploNBMixtureEMIt(const std::map<uint, double>& data, uint xMin,
                             uint xMax, std::vector<double>& mu,
                             std::vector<double>& ODF, std::vector<double>& w,
                             int nStrains, double maxEps, bool genusLevel)
{
        // number of components = error + strain-only + not-strain + sos
        size_t numComp = 2 + ((nStrains > 2) ? 2 : 1) * nStrains;
        // number of NB = number of components + higher-order for sos
        size_t numNB = mu.size();

        vector<double> muPrev = mu;

        // compute, for each x, the relative contribution of the NB
        map<uint, vector<double>> fNB;
        for (const auto& element : data) {
                const uint& x = element.first;
                if (x > xMax)
                        break;

                // compute for each NB the value P(x)
                vector<double> P(numNB);
                for (size_t i = 0; i < numNB; i++)
                        P[i] = w[i]*Util::negbinomialPDF(x, mu[i], ODF[i]*mu[i]);

                // compute relative fractions
                fNB[x] = vector<double>(numNB, 0.0);
                double sum = accumulate(P.begin(), P.end(), 0.0);
                if (sum > 0.0)          // avoid numeric underflow
                        for (size_t i = 0; i < numNB; i++)
                                fNB[x][i] = P[i] / sum;
                else
                        fNB[x].back() = 1.0;  // last component = 1
        }

        // We fit a negative binomial (NB) distribution to the error histogram.
        // The error histogram is always truncated: a) coverage zero is never
        // observed and b) low-coverage k-mers might have been removed by
        // a preprocessing tool like BCALM. We therefore use EM to fit the NB
        // parameters and infer the missing values from the spectrum.
        map<uint, double> errData;
        for (const auto& element : data) {
                uint x = element.first;
                if (x > xMax)
                        break;
                errData[x] = fNB[x][0] * data.at(x);
        }

        // We truncate the node histogram to the provided -abundance-min value.
        // This is relevant when using q-mer counts.
        for (uint k = 0; k < xMin; k++)
                errData.erase(k);

        Util::fitTruncNegBinomEM(errData, mu[0], ODF[0], w[0],
                                 maxEps, 10.0 / maxEps);

        // for each non-error NB
        for (size_t i = 1; i < numNB; i++) {
                // get the right fraction of the data
                map<uint, double> cData;
                for (const auto& element : data) {
                        const uint& x = element.first;
                        if (x > xMax)
                                break;
                        cData[x] = fNB[x][i] * data.at(x);
                }

                // data the data weight and mean
                w[i] = 0.0;
                double xySum = 0.0;
                for (const auto& element : cData) {
                        w[i] += element.second;
                        xySum += element.first * element.second;
                }
                mu[i] = xySum / w[i];
        }

        // take into account the relationships between weights
        double sum = 0.0;
        if(genusLevel) {
                for (size_t i = 1; i < nStrains + 1; i++)
                        sum += w[i];
                for (size_t i = 1; i < nStrains + 1; i++)
                        w[i] = sum / (nStrains);
        } else {
                for (size_t i = 1; i < numComp - 1; i++)
                        sum += w[i];
                for (size_t i = 1; i < numComp - 1; i++)
                        w[i] = sum / (numComp - 2);
        }

        // take into account the relationships between mean values
        int rows = numComp - 1;
        vector<double> A(rows * nStrains, 0.0);
        vector<double> rhs(rows);

        // create matrix A and right-hand side vector
        for (int i = 0; i < nStrains; i++) {    // strain i unique
                A[i*rows+i] = sqrt(w[i+1]);
                rhs[i] = sqrt(w[i+1])*mu[i+1];
        }
        if (nStrains > 2) {                     // sum \ strain i
                for (int i = nStrains; i < 2*nStrains; i++) {
                        for (int j = 0; j < nStrains; j++) {
                                if (i == j+nStrains)
                                        continue;
                                A[j*rows+i] = sqrt(w[i+1]);
                        }
                        rhs[i] = sqrt(w[i+1])*mu[i+1];
                }
        }
        for (int j = 0; j < nStrains; j++)       // sum peak
                A[j*rows+rows-1] = sqrt(w[rows]);
        rhs[rows-1] = sqrt(w[rows])*mu[rows];

        /*cout << "rows: " << rows << endl;
        for (int i = 0; i < rows; i++) {
                for (int j = 0; j < nStrains; j++)
                        cout << A[j*rows+i] << "\t";
                cout << rhs[i] << "\n";
        }*/

        // least-squares solution
        char N = 'N'; int nrhs = 1, info;
        int lwork = 2* rows * nStrains;
        vector<double> work(lwork);
        dgels_(&N, &rows, &nStrains, &nrhs, A.data(), &rows, rhs.data(),
               &rows, work.data(), &lwork, &info);

        // plug the result back into the mu vector
        sum = accumulate(rhs.begin(), rhs.begin() + nStrains, 0.0);
        for (int i = 0; i < nStrains; i++)      // strain i unique
                mu[i+1] = max(rhs[i], mu[0]);
        if (nStrains > 2)                       // sum \ strain i
                for (int i = 0; i < nStrains; i++)
                        mu[i+nStrains+1] = max(sum-rhs[i], mu[0]);
        mu[rows] = max(mu[0], sum);                         // sum

        // compute the variance
        for (size_t i = 1; i < numComp; i++) {
                // get the right fraction of the data
                map<uint, double> cData;
                for (const auto& element : data) {
                        const uint& x = element.first;
                        if (x > xMax)
                                break;
                        cData[x] = fNB[x][i] * data.at(x);
                }

                // data the data weight and mean
                double total = 0.0, count = 0.0;
                for (const auto& element : cData) {
                        const double x = element.first;
                        const double y = element.second;

                        count += y;
                        total += (x - mu[i]) * (x - mu[i]) * y;
                }

                double var = total/count;
                ODF[i] = var / mu[i];
                if (ODF[i] < 1.0)
                        ODF[i] = 1.0;   // ODF cannot be smaller than 1
        }

        // average the ODF over all strains
        double total = 0.0, count = 0.0;
        for (size_t i = 1; i < numComp; i++) {
                total += sqrt(w[i]) * ODF[i];
                count += sqrt(w[i]);
        }
        double avgODF = total / count;
        for (size_t i = 1; i < numNB; i++)
                ODF[i] = avgODF;

        //for (int i = 0; i < mu.size(); i++)
        //        cout << mu[i] << "\t" << ODF[i] << "\t" << w[i] << endl;

        return Util::relDiff(mu, muPrev);
}

tuple<int, double> CovModel::haploNBMixtureEM(const std::map<uint, double>& data, uint xMin,
                           uint xMax, std::vector<double>& mu,
                           std::vector<double>& ODF, std::vector<double>& w,
                           int nStrains, double maxEps,
                           int maxIterations, bool genusLevel)
{
        double eps = 1.0;

        // main EM loop
        for (int itCount = 1; itCount <= maxIterations; itCount++) {
                eps = haploNBMixtureEMIt(data, xMin, xMax, mu, ODF, w,
                                                nStrains, maxEps, genusLevel);
                if (eps < maxEps)       // check for convergence
                        return tuple<int, double> (itCount, eps);

        }

        return tuple<int, double> (maxIterations, eps);
}

CovModel::CovModel(const std::string& filename)
{
        read(filename);
}

std::vector<int> CovModel::getExpMult(double obsCov) const
{
        // TODO find a smart way to do as little computations as possible
        //      Try to compute only for all strains equal multiplicity + for +1 on one or more strains
        //      Or maybe even only all strains equal multiplicity? (we have a margin afterwards in the CRF anyway?)

        vector<int> card(fractions.size(),Multiplicity::maxM+1);

        int numEntries = accumulate(card.begin(), card.end(), 1, multiplies<int>());
        double P_max = 0.0;
        vector<int> expM(fractions.size(), 0);

        //fill P appropriately
        Assignment ass(card);
        if (obsCov > Multiplicity::maxM*sosLambda){
                for (int i = 0; i < numEntries; i++) {
                        double p = (std::find(ass.begin(), ass.end(), Multiplicity::maxM) == ass.end()) ? getLogProb(obsCov, ass) : maxMultPVal;
                        if(p > P_max){
                                P_max = p;
                                expM = ass;
                        }
                        ass.increment();
                }
        }else{
                for (int i = 0; i < numEntries; i++) {
                        double p = getLogProb(obsCov, ass);
                        if(p > P_max){
                                P_max = p;
                                expM = ass;
                        }
                        ass.increment();
                }
        }
        return expM;
}

Multiplicity CovModel::getMultSoft(double obsCov) const
{
        vector<int> card(fractions.size(),Multiplicity::maxM+1);

        int numEntries = accumulate(card.begin(), card.end(), 1, multiplies<int>());
        vector<double> P(numEntries,0.0);

        //fill P appropriately
	Assignment ass(card);
	if (obsCov > Multiplicity::maxM*sosLambda){
		for (int i = 0; i < numEntries; i++) {
			P[i] = (std::find(ass.begin(), ass.end(), Multiplicity::maxM) == ass.end()) ? getLogProb(obsCov, ass) : maxMultPVal;
			ass.increment();
		}
	}else{
		for (int i = 0; i < numEntries; i++) {
			P[i] = getLogProb(obsCov, ass);
			ass.increment();
		}
	}
        return Multiplicity(card, P);
}

double CovModel::getLogProb(double obsCov, const vector<int>& mult) const
{
        // make sure there is a strain multiplicity assignment for each strain
        assert(mult.size() == fractions.size());

        // if all strains have mult 0 we are dealing with an error-node/arc
        if (accumulate(mult.begin(), mult.end(), 0, plus<int>()) == 0)
                return logw[0] + Util::logNegbinomialPDF(obsCov, errLambda, errODF * errLambda);

        double myLambda = 0;
        for (int i = 0; i < mult.size(); i++)
                myLambda += mult[i] * fractions[i];
        myLambda *= sosLambda;
        double w = getLWeight(mult);
        return w + Util::logNegbinomialPDF(obsCov, myLambda, sosODF * myLambda);
}

void CovModel::write(const std::string& covFilename) const
{
        ofstream covFile(covFilename.c_str());

        covFile << errLambda << "\t" << errODF << "\t"
                << sosLambda << "\t" << sosODF << "\n" << fractions.size();
        for (double frac: fractions)
                covFile << "\t" << frac;

        covFile << "\n" << logw.size() << "\t" << maxMult << "\n";

        for (size_t i = 0; i < logw.size(); i++) {
                covFile << exp(logw[i]);
                covFile << (i % maxMult == 0 ? "\n" : "\t");
        }

        covFile.close();
}

void CovModel::read(const std::string& covFilename)
{
        ifstream ifs(covFilename.c_str());
        size_t numStrains;
        ifs >> errLambda >> errODF >> sosLambda >> sosODF >> numStrains;
        fractions.resize(numStrains);
        for (size_t f = 0 ; f < numStrains ; f++)
                ifs >> fractions[f];
        size_t w_size;
        ifs >> w_size >> maxMult;
        logw.resize(w_size);
        for (size_t i = 0; i < w_size; i++) {
                double w;
                ifs >> w;
                logw[i] = log(w);
        }
        ifs.close();

        Assignment mult(vector<int>(numStrains, Multiplicity::maxM+1));

        double maxLProb = std::numeric_limits<double>::lowest();
        for (int i = 0; i < getK(); i++) {
                if (find(mult.begin(), mult.end(), Multiplicity::maxM) == mult.end()){
                        mult.increment();
                        continue;
                }
                double multP = getLogProb(((double)Multiplicity::maxM) * sosLambda, mult);
                if (multP > maxLProb)
                        maxLProb = multP;
                mult.increment();
        }
        setMaxMultP(maxLProb);
}

void CovModel::writeGnuplot(const string& baseFilename,
                            const map<unsigned int, double>& hist,
                            int maxPlotMult) const
{
        // put everything in a .dat file
        ofstream ofs((baseFilename + ".dat").c_str());
        for (int c = 0; c <= (maxMult + 0.5) * sosLambda; c++) {
                ofs << c << "\t";
                Assignment ass(vector<int>(getNumStrains(), maxMult + 1));
                for (int m = 0; m < logw.size(); m++, ass.increment())
                        ofs << getProb(c, ass) << "\t";
                ofs << (hist.find(c) == hist.end() ? 0 : hist.at(c));
                ofs << "\n";
        }
        ofs.close();

        // write a sensible .gnuplot file
        //int maxPlotMult = 1;  // change this for extra components
        int xMax = (maxPlotMult + 0.5) * sosLambda;
        double yweight = getWeight( vector<int>(getNumStrains(), 1));
        for (int i = 0; i < getNumStrains(); i++){
                vector<int> temp(getNumStrains(), 0);
                temp[i] = 1;
                double w = getWeight(temp);
                if (w > yweight)
                        yweight = w;
        }
        int yMax = 3.0 * ((getProb(sosLambda, vector<int>(getNumStrains(), 1))/getWeight( vector<int>(getNumStrains(), 1))) * yweight);

        ofs.open((baseFilename + ".gnuplot").c_str());
        ofs << "set terminal pdf\n"
            << "set output \"" << baseFilename << "spectrum.pdf\"\n"
            << "set title \"Coverage histogram with fitted model\"\n"
            << "set xlabel \"coverage\"\n"
            << "set ylabel \"number of observations\"\n"
            << "set xrange [0:" << xMax << "]\n"
            << "set yrange [0:" << yMax << "]\n"
            << (maxPlotMult == 1 ? "set key top left\n" : "set key top right\n")
            << "plot \"" << baseFilename << ".dat\" using 1:" << 2+logw.size()
            << " title \'k-mer data\' with boxes,\\\n"
            << "\t\"" << baseFilename << ".dat\" using 1:2"
               " title \'seq. errors\' with lines,\\\n";

        Assignment ass(vector<int>(getNumStrains(), maxMult + 1));
        ass.increment();        // skip sequencing errors (already printed)
        for (int m = 1; m < logw.size(); m++, ass.increment()) {
                // only print a limited number of multplicities
                if (*max_element(ass.begin(), ass.end()) > maxPlotMult)
                        continue;

                ofs << "\t\"" << baseFilename << ".dat\" using 1:"
                    << 2+m << " title \'mult = (";
                for (int m : ass)
                        ofs << " " << m;
                ofs << " )\' with lines,\\\n";
        }
        ofs.close();
}

void CovModel::writeGnuplotInit(const string& baseFilename,
                            const map<unsigned int, double>& hist,
                            const vector<double>& initW,
                            const double& lambdaHighM,
                            const double& ODFhighM) const
{
    double maxProb = 0.0;
    double sumInitW = accumulate(initW.begin(),initW.end(),0.0);
    // Do put everything in .dat
    ofstream ofs((baseFilename + ".dat").c_str());
    // Using wInit instead of weights
    for (int c = 1; c < 2.0 * sosLambda; c++) {
        ofs << c << "\t";
        Assignment ass(vector<int>(getNumStrains(),2)); //only consider mult 0 or 1 for init
        ass.increment(); //don't consider (0,0,...)
        double mixt = 0.0;
        for (int m = 0; m < initW.size()-1; m++, ass.increment()){
            double mLambda = 0.0;
            for (int strain = 0; strain < getNumStrains(); strain++) {
                mLambda += ass[strain] * fractions[strain]; //ass[strain] is 0 or 1
            }
            mLambda *= sosLambda;
            double prob = log(initW[m]) + Util::logNegbinomialPDF(double(c), mLambda, sosODF * mLambda);
            ofs << exp(prob) << "\t";
            if (exp(prob)>maxProb)
                maxProb = exp(prob);
            mixt += exp(prob);
        }

        //highM (extra bubble)
        double prob = log(initW[initW.size()-1]) + Util::logNegbinomialPDF(double(c), lambdaHighM, ODFhighM * lambdaHighM);
        ofs << exp(prob) << "\t";
        if (exp(prob)>maxProb)
            maxProb = exp(prob);
        mixt += exp(prob);

        ofs << (hist.find(c) == hist.end() ? 0 : hist.at(c)) << "\t";
        ofs << mixt; //mixture
        ofs << "\n";
    }

    ofs.close();

    // But simplify gnuplot?
    ofs.open((baseFilename + ".gnuplot").c_str());
    ofs << "set terminal pdf\n"
        << "set output \"" << baseFilename << "spectrum.pdf\"\n"
        << "set title \"Coverage histogram with fitted model\"\n"
        << "set xlabel \"coverage\"\n"
        << "set ylabel \"number of observations\"\n"
        << "set xrange [0:" << int(2.0 * sosLambda) << "]\n"
        << "set yrange [0:" << int(maxProb*1.25) << "]\n"
        << "plot \t\"" << baseFilename << ".dat\" using 1:"  << initW.size()+2 << " title \'k-mers\' with boxes, \\\n";
    Assignment ass(vector<int>(getNumStrains(),2)); // Show only strain_mult 0,1 combinations

    int maxMultip = 1;
    while(accumulate(ass.begin(), ass.end(), 1, multiplies<int>()) != maxMultip) {
        ass.increment(); // skip sequencing errors (allready printed)
        ofs << "\t\"" << baseFilename << ".dat\" using 1:";
        int column = 1;
        for (size_t i = 0; i < getNumStrains(); i++) {
            column += std::min<size_t>(ass[i], maxMult) * pow(2,i);
        }
        ofs << column << " title \'mult = (";
        for ( int m : ass)
            ofs << " " << m;
        ofs << " )\' with lines,\\\n";
    }
    ofs << "\t\"" << baseFilename << ".dat\" using 1:" << initW.size()+1 << " title \'highM\' with lines, \\\n";
    ofs << "\t\"" << baseFilename << ".dat\" using 1:" << initW.size()+3 << " title \'Mixture\' with lines, \\\n"; //TODO
    //ofs << "\t\"" << baseFilename << ".dat\" using 1:" << initW.size()+2 << " title \'k-mers\' with boxes,\n";
    ofs.close();
}

std::ostream &operator<<(std::ostream& out, const CovModel& c)
{
        out << fixed;
        out.precision(2);
        out << "[error: " << c.errLambda << " (" << c.errODF << "); "
            << "sum-of-strains: " << c.sosLambda << " (" << c.sosODF << ")]\n";
        out << "\tstrain fractions: [";
        for (size_t i = 0; i < c.getNumStrains(); i++)
                out << " " << c.fractions[i] << "%";
        out << " ]\n\tweights: [";
        out.precision(1);
        Assignment ass(vector<int>(c.getNumStrains(), c.maxMult + 1));
        for (size_t i = 0; i < c.logw.size(); i++, ass.increment())
                out << ((i%(c.maxMult+1) == 0) ? "\n\t\t" : " ")
                    << exp(c.logw[i]) << " (" << ass << " )";
        out << "]";

        return out;
}

double CovModel::getMSE(const map<unsigned int, double>& hist) const
{
        double mse = 0.0;
        int numobs = 0;
        for (const auto& element : hist) {
                unsigned int x = element.first;
                double y = element.second;

                double est = 0.0;
                Assignment ass(vector<int>(getNumStrains(), maxMult + 1));
                for (int m = 0; m < logw.size(); m++, ass.increment())
                        est += exp(getLogProb(x, ass));

                double diff = y-est;
                mse += diff * diff;
                numobs ++;
        }
        mse /= (numobs - getNumParams());
        return mse;

}

double CovModel::getLL(const vector<double>& obsCov) const
{
        double logL = 0.0;
        double sum_w = 0.0;
        for (size_t i = 0; i < logw.size(); i++)
                sum_w += exp(logw[i]);
        sum_w = log(sum_w);

        for (size_t i = 0; i < obsCov.size(); i++){
                Assignment mult(vector<int>(getNumStrains(),getC()));
                double toLog = DOUBLE_SMALL;
                for (int m = 0; m < logw.size(); m++, mult.increment()){
                        toLog += exp(getLogProb(obsCov[i],mult));
                }
                logL += (log(toLog) - sum_w);
        }
        return logL;
}

double CovModel::getCompleteLL(NodeMap<Multiplicity> mult, const DBGraph& dBG){
        double sum_w = 0.0;
        for (size_t i = 0; i < logw.size(); i++)
                sum_w += exp(logw[i]);
        sum_w = log(sum_w);
        
        double logL = 0.0;
        for (auto& it: mult){
                vector<int> ass = it.second.getExpMult();
                logL += (getLogProb(dBG.getSSNode(it.first.getNodeID()).getAvgCov(), ass) - sum_w); 
        }
        return logL;
}

double CovModel::getCompleteLL(EdgeMap<Multiplicity> mult, const DBGraph& dBG){
        double sum_w = 0.0;
        for (size_t i = 0; i < logw.size(); i++)
                sum_w += exp(logw[i]);
        sum_w = log(sum_w);
        
        double logL = 0.0;
        for (auto& it: mult){
                vector<int> ass = it.second.getExpMult();
                logL += (getLogProb(dBG.getArc(dBG.getArcID(it.first)).getCov(), ass) - sum_w); 
        }
        return logL;
}


double CovModel::getEntropy(const std::vector<double>& obsCov) const
{
        double entropy = 0.0;
        for (size_t i = 0; i < obsCov.size(); i++){
                Assignment mult(vector<int>(getNumStrains(),getC()));
                vector<double> unNormProb(logw.size());
                for (int m = 0; m < logw.size(); m++, mult.increment())
                        unNormProb[m] = getLogProb(obsCov[i],mult);
                //normalise
                double logSum = unNormProb.front();
                for (size_t i = 1; i < unNormProb.size(); i++) {
                        double m = max<double>(logSum, unNormProb[i]);
                        logSum = m + log( exp(unNormProb[i]-m) + exp(logSum-m) );
                }

                for (size_t i = 0; i < unNormProb.size(); i++){
                        double logTau = unNormProb[i] - logSum;
                        entropy += exp(logTau)*logTau;
                }
        }
        return -entropy;
}

double CovModel::getCRFassignmentEntropy(NodeMap<Multiplicity> mult)
{
        double entropy = 0.0;
        for (auto& it: mult){
                it.second.normalize();
                Assignment ass(vector<int>(getNumStrains(),getC()));
                for (int i = 0; i < logw.size(); i++, ass.increment()){
                        double tau = max(DOUBLE_SMALL,it.second[ass]);
                        entropy += (log(tau)*tau);
                }
        }
        return -entropy;
}

double CovModel::getCRFassignmentEntropy(EdgeMap<Multiplicity> mult)
{
        double entropy = 0.0;
        for (auto& it: mult){
                it.second.normalize();
                Assignment ass(vector<int>(getNumStrains(),getC()));
                for (int i = 0; i < logw.size(); i++, ass.increment()){
                        double tau = max(DOUBLE_SMALL,it.second[ass]);
                        entropy += (log(tau)*tau);
                }
        }
        return -entropy;
}


// ============================================================================
// DBGRAPH
// ============================================================================

void DBGraph::rescaleKmerToNodeCts(const vector<NodeRep>& nodes, CovModel& nodeModel){
        vector<double> wN(nodeModel.getK(), 0.0);
        map<unsigned int, double> errorHist;
        vector<int> zeroMult = vector<int>(nodeModel.getNumStrains(),0);
        
        double oldZero = nodeModel.getWeight(zeroMult);
        double oldOne = nodeModel.getWeight(vector<int>(nodeModel.getNumStrains(),1));
        
        for (size_t i = 0; i < nodes.size(); i++) {
                SSNode node = getSSNode(nodes[i]);
                double cov = node.getAvgCov(); //get average kmer coverage
                vector<int> expMult = nodeModel.getExpMult(cov);
                
                if (cov > (nodeModel.getMaxMult()+1)*nodeModel.getSosLambda()) //error node
                        continue;
                
                vector<double> P(nodeModel.getK());
                Assignment mult(vector<int>(nodeModel.getNumStrains(), nodeModel.getC()));
                for (int idx=0; idx < nodeModel.getK(); idx ++){
                        P[idx] = nodeModel.getProb(cov,mult);
                        mult.increment();
                }
                
                double sum = accumulate(P.begin(), P.end(), 0.0);
                double p_zero = P[0]/sum;
                wN[0] += p_zero;
                /*if (p_zero > DOUBLE_SMALL) {
                        double f = node.getAvgCov() - floor(node.getAvgCov());
                        errorHist[node.getAvgCov()] += (1.0-f) * p_zero;// * node.getMarginalLength();
                        errorHist[node.getAvgCov() + 1] += f * p_zero;// * node.getMarginalLength();
                }*/
                /*if (accumulate(expMult.begin(), expMult.end(), 0) == 0 || cov > (nodeModel.getMaxMult()+1)*nodeModel.getSosLambda()) //error node
                        continue;*/
                //counts for weight matrix
                for (int idx=1; idx < nodeModel.getK(); idx ++){ //don't compute wNode[0], already computed
                        wN[idx] += P[idx]/sum;
                }
        }
        
        /*// We truncate the node histogram to the provided -abundance-min value.
        // This is relevant when using q-mer counts.
        unsigned int smallSupport = (getAbundanceMin() < 0)?
        errorHist.begin()->first : getAbundanceMin();
        for (unsigned int k = 0; k < smallSupport; k++)
                errorHist.erase(k);
        
        double nodeErrorLambda = nodeModel.getErrLambda();
        double nodeErrorODF = nodeModel.getErrorODF();
        double nodeErrorWeight = 2.0; //nodeModel.getWeight(zeroMult);
        
        // The below EM procedure might not convergence, but the obtained NB
        // parameters should be sensible in all cases: only small support
        // values are inferred, therefore, the mean should never be larger than
        // the mean of the error histogram (and thus be close to zero).
        int maxIter = 10.0 / settings.getEMConvEps();
        int nIter = Util::fitTruncNegBinomEM(errorHist, nodeErrorLambda,
                                             nodeErrorODF, nodeErrorWeight,
                                             settings.getEMConvEps(), maxIter);
        
        cout << "\tEM algorithm to fit edge error model converged after " << nIter << " iterations" << endl;
        
        wN[0] = nodeErrorWeight;
        
        if (nIter > maxIter)
                cout << "\tWARNING: EM algorithm to fit node error model did not converge\n";*/
        
        Assignment mult(vector<int>(nodeModel.getNumStrains(), nodeModel.getC()));
        for (int idx=0; idx < nodeModel.getK(); idx ++){
                nodeModel.setWeight(mult,max(1.0,wN[idx]));
                mult.increment();
        }
        
        wN[0] = nodeModel.getWeight(vector<int>(nodeModel.getNumStrains(),1))/oldOne * oldZero;
        nodeModel.setWeight(zeroMult, max(1.0,wN[0]));
        
        //save All relevant weights
        double pWeight = 0.0;
        int totCase = 0;
        int numStrains = nodeModel.getNumStrains();
        for (size_t i = 0; i < numStrains; i++) {
                vector<int> a(numStrains, 0); a[i] = 1;
                pWeight += ( nodeModel.getWeight(a) );
                totCase ++;
        }
        if( ! settings.useGenusLevelModel()){
                if (numStrains > 2)
                        for (size_t i = 0; i < numStrains; i++) {
                                vector<int> a(numStrains, 1); a[i] = 0;
                                pWeight += ( nodeModel.getWeight(a) );
                                totCase ++;
                        }
        }
        pWeight /= totCase;
        vector<double> allStrains(nodeModel.getC());
        for (int i = 0; i < nodeModel.getC(); i++)
                allStrains[i] = nodeModel.getWeight(vector<int>(numStrains, i));
        
        
        //make sure new W has same layout as init estimation
        Assignment ass(vector<int>(nodeModel.getNumStrains(), nodeModel.getC()));
        for (int idx=1; idx < nodeModel.getK(); idx ++){
                ass.increment();
                nodeModel.setWeight(ass,1.0);
        }
        
        
        for (size_t i = 0; i < numStrains; i++) {
                vector<int> a(numStrains, 0); a[i] = 1;
                nodeModel.setWeight(a,pWeight);
        }
        if (numStrains > 2)
                for (size_t i = 0; i < numStrains; i++) {
                        vector<int> a(numStrains, 1); a[i] = 0;
                        settings.useGenusLevelModel() ? nodeModel.setWeight(a, max(1.0,wN[ass.assignmentToIndex(a)])): nodeModel.setWeight(a, pWeight);
                }
                
        for (int i = 0; i < nodeModel.getC(); i++)
                nodeModel.setWeight(vector<int>(numStrains, i), allStrains[i]);
        
                
        
}

double DBGraph::estimateInitKmerCov(double errLambda, double p) const
{
        // sanity checks
        assert(errLambda > 0.0);
        assert(p > 0.0);
        assert(p < 1.0);

        // Given a Poisson distribution for the error model, find a cutoff
        // value for the coverage for which the probability of observing
        // a coverage is less than p under this error model
        double cutoff = ceil(errLambda);
        for ( ; cutoff < 10.0 * errLambda; cutoff++)
                if (Util::poissonPDF((unsigned int)cutoff, errLambda) < p)
                        break;

        double totCoverage = 0, totSize = 0;
        for ( NodeID id = 1; id <= numNodes; id++ ) {
                SSNode node = getSSNode(id);
                if (!node.isValid())
                        continue;
                if (node.getAvgCov() < cutoff)
                        continue;
                totCoverage += node.getCov();
                totSize += node.getMarginalLength();
        }

        if (totSize > 0)
                return (double)totCoverage / (double)totSize;
        else
                return 2.0 * errLambda;      // pathological case
}

CovModel DBGraph::getInitCovModel(const Settings& settings, const int nStrains)
{
        // get initial error model
        double initErrCov = settings.getInitErrCov();
        double initErrorODF = settings.getInitErrODF();

        // get initial sum-of-strains (sos) model
        double initSosCov = settings.getInitCov();
        if (initSosCov < 0.0)           // if not user-provided, estimate
                initSosCov = estimateInitKmerCov(initErrCov);
        double initStrainODF = settings.getInitODF();

        // initialize the strain coverages
        vector<double> strainFrac = settings.getInitFractions();
        if( strainFrac.empty() || strainFrac.size() != nStrains){
                strainFrac.resize(nStrains); 
                int sum = 0;
                for (int i = 0; i < nStrains; i++) {
                        strainFrac[i] = i+1;
                        sum += (i+1);
                }
                for (int i = 0; i < nStrains; i++)
                        strainFrac[i] /= sum;
        }

//        size_t nStrains = strainFrac.size();
        if(settings.useGenusLevelModel())
                initSosCov *= nStrains; // average of largest k-mers is average of nStrains peaks, we want the sum of the peaks for initSosCov
        vector<double> initStrainCov;   // lambda of the strains
        for (size_t i = 0; i < nStrains; i++)
                initStrainCov.push_back(strainFrac[i] * initSosCov);

        // initialize the weights
        size_t maxMult = settings.getCRFMaxMult();
        const double pseudoCount = 1.0;
        vector<double> initW(Util::pow(maxMult + 1, nStrains), pseudoCount);

        CovModel covModel(initErrCov, initErrorODF, initStrainCov,
                          initStrainODF, initW, maxMult);

        covModel.setWeight(vector<int>(nStrains, 0), 1e3);      // error weight
        
        if(settings.useGenusLevelModel()){
                for (size_t i = 0; i < nStrains; i++) {
                        // weight of strain i
                        vector<int> a(nStrains, 0);
                        a[i] = 1;
                        covModel.setWeight(a, 1e4);
                        
                        if (nStrains <= 2)
                                continue;
                        
                        // weight of not-strain i
                        a = vector<int>(nStrains, 1);
                        a[i] = 0;
                        covModel.setWeight(a, 1e2);
                }
                
                for (size_t i = 0; i < maxMult; i++) {  // weight of the sum-of-strains
                        vector<int> a(nStrains, i+1);
                        covModel.setWeight(a, max(10 / pow(10, i), pseudoCount));
                } 
        } else { 
                for (size_t i = 0; i < nStrains; i++) {
                        // weight of strain i
                        vector<int> a(nStrains, 0);
                        a[i] = 1;
                        covModel.setWeight(a, 1e2);

                        if (nStrains <= 2)
                                continue;

                        // weight of not-strain i
                        a = vector<int>(nStrains, 1);
                        a[i] = 0;
                        covModel.setWeight(a, 1e2);
                }

                for (size_t i = 0; i < maxMult; i++) {  // weight of the sum-of-strains
                        vector<int> a(nStrains, i+1);
                        covModel.setWeight(a, max(1e3 / pow(10, i), pseudoCount));
                }
        }

        return covModel;
}

int DBGraph::getInitialLambdaEstimates(double cutoff, int nStrains, vector<double>& initLambdas, double& lambdaHighM, vector<double>& initW, double& initODF, double& ODFhighM, int maxIter) const
{

        assert(initW.size()==(pow(2,nStrains)));
        vector<int> combinations(nStrains,2);
        double maxCutoff = 2*accumulate(initLambdas.begin(), initLambdas.end(), 0.0);
        int highMidx = initW.size()-1;

        //EM algorithm: find lambdas in mixture model (assume only mult 0 or 1 exist for each strain)

        map<int, vector<double>> nodeProbs;
        vector<double> zeroProb(initW.size(),0.0);
        for (NodeID id = 1; id <= numNodes; id++){
            SSNode node = getSSNode(id);
            double obsCov = node.getAvgCov();
            if (!node.isValid())
                continue;
            if (obsCov < cutoff)
                continue; //cov too low
            if (obsCov > maxCutoff)
                continue; //cov too high
            nodeProbs.insert({id,zeroProb});
        }

        int iter = 0;
        while(iter < maxIter){
            iter ++;
            double nodesChanged = 0.0;

            //normalize weights
            double sumInitW = accumulate(initW.begin(), initW.end(),0.0);
            for (int i = 0; i < initW.size(); i++){
                initW[i] /= sumInitW;
            }



            //E-step
            for (auto it = nodeProbs.begin(); it != nodeProbs.end(); it++ ) {
                SSNode node = getSSNode(it->first);
                double obsCov = node.getAvgCov();
                double prevMaxIdx = distance(it->second.begin(), max_element(it->second.begin(), it->second.end()));

                Assignment ass(combinations);
                ass.increment(); //don't consider (0,0,...,0)

                vector<int> logProb(initW.size());
                double sumProb = 0.0; //for normalization

                //compute probabilities
                for (int i = 0; i < initW.size()-1; i++) {
                    double expCov = 0;
                    for (int strain = 0; strain < ass.size(); strain++)
                        expCov += ass[strain] * initLambdas[strain];

                    logProb[i] = log(initW[i]) + Util::logNegbinomialPDF(obsCov, expCov, initODF * expCov);
                    sumProb += exp(logProb[i]);
                    ass.increment();
                }
                // compute probability for highM
                logProb[highMidx] = log(initW[highMidx]) + Util::logNegbinomialPDF(obsCov, lambdaHighM, ODFhighM * lambdaHighM);
                sumProb += exp(logProb[highMidx]);

                //normalize prob
                for (int i = 0; i < initW.size(); i++){
                    it->second[i] = logProb[i] - log(sumProb);
                }

                double maxIdx = distance(it->second.begin(), max_element(it->second.begin(), it->second.end()));
                if (prevMaxIdx != maxIdx)
                    nodesChanged++;
            }


            //Mstep

            //initialize W
            vector<double> newW(initW.size(),0.0);

            //initialize mtm & mtx
            vector<double> mtm(nStrains*nStrains, 0.0);
            vector<double> mtX(nStrains, 0.0);

            // for each (1,0,0), (0,1,0), (0,0,1) add PC
            for (int strain = 0; strain < nStrains; strain++){
                mtm[strain * (nStrains+1)] = 1.0;
                mtX[strain] = initLambdas[strain];
                newW[pow(2,strain) - 1] = 1.0;
            }

            //initialize parameters highM (>1)
            int nrHighM = 0;
            double avgCovHighM = 0.0;

            for (auto it = nodeProbs.begin(); it != nodeProbs.end(); it++ ) {
                SSNode node = getSSNode(it->first);
                double maxProb = std::numeric_limits<double>::lowest();
                int idx = 0;
                vector<int> expMult(nStrains,0);
                Assignment ass(combinations);
                ass.increment(); //don't consider (0,0,...,0)

                for (int i = 0; i < initW.size()-1; i++){
                    //fill newW :soft assignments
                    newW[i] += exp(it->second[i])*node.getMarginalLength();

                    //get expMult
                    if (it->second[i] > maxProb){
                        maxProb = it->second[i];
                        idx = i;
                        expMult = ass;
                    }

                    ass.increment();
                }

                //parameters for highM (>1)
                newW[highMidx] += exp(it->second[highMidx])*node.getMarginalLength();
                if (it->second[highMidx] > maxProb){ //expected highM
                    //hard assignments to compute lambdaHighM
                    nrHighM += node.getMarginalLength();
                    avgCovHighM += node.getCov();
                }else{ //expected expMult
                    // fill mtm & mtX for OLS
                    for (int strain = 0; strain < nStrains; strain++) {
                        mtX[strain] += expMult[strain] * node.getMarginalLength() * node.getCov();
                        for (size_t j = strain; j < nStrains; j++)
                            mtm[strain * nStrains + j] += expMult[strain] * node.getMarginalLength() * expMult[j] * node.getMarginalLength();
                    }
                }
            }
            //mtm is symmetric
            for (size_t j = 0; j < nStrains; j++){
                for (size_t i = j+1; i < nStrains; i++) {
                    mtm[i*nStrains + j] = mtm[j*nStrains + i];
                }
            }

            //compute new lambda estimates
            Util::choleskyLinSolve(mtm,mtX);

            //check for negative lambdas
            for (int strain = 0; strain < nStrains; strain++){
                if (mtX[strain] < 0){ //TODO: better message
                    //assign previous smallest lambda
                    mtX[strain] = *min_element(initLambdas.begin(), initLambdas.end());
                    cerr << "Warning: encountered negative lambda in initLambdaEstimates \n"
                            "set to smallest lambda in previous EM iteration \n";
                }
            }

            //new lambdaHighM
            lambdaHighM = avgCovHighM/nrHighM;

            //check if lambdaHighM > sumLambdas
            if (lambdaHighM <= accumulate(mtX.begin(), mtX.end(), 0.0)){
                lambdaHighM = 1.05 * accumulate(mtX.begin(), mtX.end(), 0.0);
                cerr << "Warning: encountered highM smaller than sum of lambdas\n"; //TODO: better message
            }

            cout << "OLS Lambda est" << "\t";
            for(int strain = 0; strain < nStrains; strain++){
                cout << mtX[strain] << "\t";
            }
            cout << "\t" << "lambdaHighM: " << lambdaHighM ;
            cout << endl;

            initLambdas = mtX;
            initW = newW;

            //compute ODF
            double nom = 0.0;
            double nomHighM = 0.0;
            for (auto it = nodeProbs.begin(); it != nodeProbs.end(); it++ ) {
                SSNode node = getSSNode(it->first);

                //ODF using soft assignments
                Assignment ass(combinations);
                ass.increment(); //don't consider (0,0,...,0)
                for (int i = 0; i < initW.size()-1; i++){ //compute ODFhighM separately
                    double expCov = 0;
                    for (size_t strain = 0; strain < nStrains; strain++) {
                        expCov += ass[strain] * initLambdas[strain];
                    }
                    double delta = node.getAvgCov() - expCov;
                    nom += exp(it->second[i]) * delta * delta * node.getMarginalLength() / expCov;
                    ass.increment();
                }
                //compute ODFhighM
                double deltaHighM = node.getAvgCov() - lambdaHighM;
                nomHighM += exp(it->second[highMidx]) * deltaHighM * deltaHighM * node.getMarginalLength()/ lambdaHighM;
            }

            initODF = nom/(accumulate(initW.begin(), initW.end(), 0.0) - initW[highMidx]); //don't consider last element
            initODF = max(1.0, initODF);
            cout << "ODF nodes : " << initODF << endl;

            ODFhighM = nomHighM/initW[highMidx];
            ODFhighM = max(1.0, ODFhighM);
            cout << "ODF highM : " << ODFhighM << endl;


            double sumLambdas = accumulate(initLambdas.begin(), initLambdas.end(), 0.0);
            double size = nodeProbs.size();
            if (2*sumLambdas > maxCutoff){ //add extra nodes
                for (NodeID id = 1; id <= numNodes; id++){
                    SSNode node = getSSNode(id);
                    double obsCov = node.getAvgCov();
                    if (!node.isValid())
                        continue;
                    if (obsCov < maxCutoff)
                        continue; //cov too low or already in map
                    if (obsCov > 2*sumLambdas)
                        continue; //cov too high
                    nodeProbs.insert({id,zeroProb});
                }
            }else{ //remove nodes
                for (NodeID id = 1; id <= numNodes; id++){
                    SSNode node = getSSNode(id);
                    double obsCov = node.getAvgCov();
                    if (!node.isValid())
                        continue;
                    if (obsCov < 2*sumLambdas)
                        continue; //keep in map
                    if (obsCov > maxCutoff)
                        continue; //cov too high
                    nodeProbs.erase(id);
                }
            }
            maxCutoff = 2*sumLambdas;

            //check convergence
            double change = nodesChanged/nodeProbs.size();
            if (nodesChanged/nodeProbs.size() < 0.001)
                return iter;

        }

        cout << "Warning: initial lambda estimates did not converge after " << iter << " iterations." << endl;
        return iter;
}

void DBGraph::covCount(const FastQRecord& rr, const KmerNPPTable& table)
{
        const string& read = rr.getRead();
        const string& qual = rr.getQual();

        NodePosPair prevNpp; NodePosition prevOff = 0;
        for (KmerIt it(read); it.isValid(); it++) {
                NodePosPair npp = table.find(it.getKmer());
                if (!npp.isValid())
                        continue;

                // increase the node coverage
                NodeID id = npp.getNodeID();
                SSNode n = getSSNode(id);
                double score = Util::phred2prob(qual, it.getOffset(),
                                                it.getOffset() + Kmer::getK());
                n.incCov(score);

                // if current kmer succeeds a valid previous kmer
                if ((prevOff+1 == it.getOffset()) && (crossesArc(prevNpp, npp))) {
                        score = Util::phred2prob(qual, prevOff,
                                                 prevOff + Kmer::getK() + 1);
                        getSSNode(prevNpp.getNodeID()).rightArc(id)->incCov(score);
                        // palindromic arcs exist only once! Don't add coverage!
                        if (prevNpp.getNodeID() != - id)
                                n.leftArc(prevNpp.getNodeID())->incCov(score);
                }

                prevNpp = npp;
                prevOff = it.getOffset();
        }
}

void DBGraph::covCountThread(FastQReader& inputs,
                             const KmerNPPTable& table)
{
        // local storage of reads
        vector<FastQRecord> readBuf;

        size_t chunkID;
        while (inputs.getNextChunk(readBuf, chunkID))
                for (const auto& readRecord : readBuf)
                        covCount(readRecord, table);
}

void DBGraph::getCovFromReads(LibraryContainer &inputs, const KmerNPPTable& table)
{
        // reset all coverages (might have been loaded from BCALM)
        for (size_t i = 1; i <= numNodes; i++)
                getSSNode(i).setCov(0);

        for (size_t i = 1; i <= numArcs; i++)
                arcs[i].setCov(0);

        // initialize Phred score lookup table
        if (settings.useQual()) {
                cout << "Using Phred quality scores (ASCII base="
                     << settings.getPhredBase() << ") to weigh k-mers\n";
                Util::enablePhred(settings.getPhredBase());
        } else {
                cout << "Not using Phred quality scores\n";
        }

        const unsigned int& numThreads = settings.getNumThreads();
        cout << "Number of threads: " << numThreads << endl;

        const size_t ws = settings.getThreadIOWorkSize();
        for (size_t i = 0; i < inputs.size(); i++) {
                string fn1, fn2;
                tie(fn1, fn2) = inputs.getFilename(i);
                
                FastQReader myReader(fn1, fn2);
                myReader.startReaderThread(ws, ws * settings.getNumThreads());

                // start worker threads
                vector<thread> workerThreads(numThreads);
                for (size_t i = 0; i < workerThreads.size(); i++)
                        workerThreads[i] = thread(&DBGraph::covCountThread,
                                                  this, ref(myReader), cref(table));
                // wait for worker threads to finish
                for_each(workerThreads.begin(), workerThreads.end(), mem_fn(&thread::join));

                myReader.joinReaderThread();
        }

        size_t numZero = 0;
        for (ArcID id = 1; id <= numArcs; id++)
                if (arcs[id].getCov() == 0)
                        numZero++;

        if (numZero > 0)
                cerr << "WARNING: found " << numZero << " arcs with coverage 0\n";
}
