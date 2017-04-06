/*
 * Copyright (C) 2016 Brian L. Browning
 *
 * This file is part of Beagle
 *
 * Beagle is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Beagle is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package vcf;

/**
 * <p>Class {@code R2Estimator} estimates the correlation between the
 * estimated allele dose and true allele dose for a set of genotypes.
 * </p>
 * <p>Instances of class {@code R2Estimator} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class R2Estimator {

    private static final double MIN_R2_DEN = 1e-8;

    private int nGenotypes = 0;
    private double sumCall = 0;
    private double sumSquareCall = 0;
    private double sumExpected = 0;
    private double sumExpectedSquare = 0;
    private double sumSquareExpected= 0;
    private double sumCallExpected = 0;

    /**
     * Constructs a new {@code R2Estimator} instance.
     */
    public R2Estimator() {
    }

    /**
     * Clears all genotype data and sets the number of genotype with
     * allele dose data to 0.
     */
    public void clear() {
        this.nGenotypes = 0;
        this.sumCall = 0.0;
        this.sumSquareCall = 0.0;
        this.sumExpected = 0.0;
        this.sumExpectedSquare = 0.0;
        this.sumSquareExpected = 0.0;
        this.sumCallExpected = 0.0;
    }

    /**
     * Returns the current number of genotypes with allele dose data.
     * @return the current number of genotypes with allele dose data
     */
    public int nGenotypes() {
        return nGenotypes;
    }

    /**
     * Adds the specified allele dose probabilities for a genotype to the stored
     * allele dose data.
     * @param doseProbs an array of length 3 whose {@code j}-th element
     * is the probability that the genotype contains {@code j} non-reference
     * alleles.
     * @throws IllegalArgumentException if {@code doseProbs.length != 3}
     * @throws IllegalArgumentException if any element of {@code doseProbs} is
     * less than 0
     * @throws IllegalArgumentException if the sum of the elements in
     * {@code doseProbs} differs from 1.0 by more than {@code 1e-5}
     * @throws NullPointerException if {@code doseProbs == null}
     */
    public void addSampleData(double[] doseProbs) {
        if (doseProbs.length != 3) {
            throw new IllegalArgumentException(String.valueOf(doseProbs));
        }
        ++nGenotypes;
        int call = checkDataAndReturnMaxIndex(doseProbs);
        double exp = (doseProbs[1] + 2*doseProbs[2]);
        double expSquare = (doseProbs[1] + 4*doseProbs[2]);
        sumCall += call;
        sumSquareCall += call*call;
        sumExpected += exp;
        sumExpectedSquare += expSquare;
        sumSquareExpected += (exp*exp);
        sumCallExpected += (call*exp);
    }

    private static int checkDataAndReturnMaxIndex(double[] probs) {
        int maxIndex = 0;
        double sum = 0.0;
        for (int j=0; j<probs.length; ++j) {
            if (probs[j] < 0) {
                throw new IllegalArgumentException(String.valueOf(probs[j]));
            }
            sum += probs[j];
            if (probs[j] > probs[maxIndex]) {
                maxIndex = j;
            }
        }
        if (Math.abs(sum - 1.0) > 1e-5) {
            throw new IllegalArgumentException(String.valueOf(sum));
        }
        return maxIndex;
    }

    /**
     * Returns the estimated squared correlation between the most probable
     * ALT allele dose and the true ALT allele dose for the current
     * genotype data.
     * Returns 0 if the marker is monomorphic or if the most probable ALT
     * allele dose is monomorphic.
     *
     * @return the estimated squared correlation between the most likely
     * allele dose and the true allele dose
     */
    public double allelicR2() {
        double f = 1.0/nGenotypes;
        double cov = sumCallExpected - (sumCall * sumExpected * f);
        double varBest = sumSquareCall - (sumCall * sumCall * f);
        double varExp = sumExpectedSquare - (sumExpected * sumExpected * f);
        double den = varBest*varExp;
        return (den < MIN_R2_DEN) ? 0.0f : (cov*cov/den);
    }

    /**
     * Returns the estimated squared correlation between the expected
     * ALT allele dose and the true ALT allele dose for the current
     * genotype data. Returns 0 if the marker is monomorphic.
     *
     * @return the estimated squared correlation between the expected ALT
     * allele dose and the true ALT allele dose
     */
    public double doseR2() {
        double f = 1.0/nGenotypes;
        double num = sumSquareExpected - (sumExpected * sumExpected * f);
        double den = sumExpectedSquare - (sumExpected * sumExpected * f);
        if (num < 0.0) {
            num = 0.0;
        }
        return (den < MIN_R2_DEN) ? 0.0f : (num / den);
    }
}
