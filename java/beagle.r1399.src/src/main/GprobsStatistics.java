/*
 * Copyright (C) 2014 Brian L. Browning
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
package main;

import vcf.Marker;
import blbutil.Const;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.Random;

/**
 * Class {@code GprobsStatistics} has methods for computing statistics
 * from posterior genotype probabilities.  Multi-allelic markers are
 * transformed to diallelic markers by combining all ALT alleles.
 *
 * The squared correlation statistics computed by this class can be derived
 * using the arguments found in Appendix 1 of
 * "Browning BL and Browning SR, Am J Hum Genet 2009;84(2):210-23".
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class GprobsStatistics {

    private static final Random random = new Random(123);
    private static final DecimalFormat df = new DecimalFormat("0.####");

    private final Marker marker;
    private final int nSamples;
    private final float[] alleleFreq;

    private float sumCall = 0;             // sum_i max_j{P(j)_i}
    private float sumSquareCall = 0;       // sum_i (max_j{P(i)_i})^2
    private float sumExpected = 0;         // sum_i E[X_i]
    private float sumExpectedSquare = 0;   // sum_i E[(X_i)^2]
    private float sumSquareExpected= 0;    // sum_i (E[X_i])^2
    private float sumCallExpected = 0;     // sum_i E[X_i] * max_j{P(j)_i}

    /**
     * Constructs a new {@code GprobsStatistics} instance from the
     * specified scaled genotype probabilities.
     * @param gv scaled sample posterior genotype probabilities.
     * @param marker a marker index
     * @throws IndexOutOfBoundsException if
     * {@code marker<0 || marker>=gv.nMarkers()}
     */
    public GprobsStatistics(GenotypeValues gv, int marker) {
        int nAlleles = gv.marker(marker).nAlleles();
        this.marker = gv.marker(marker);
        this.nSamples = gv.nSamples();
        this.alleleFreq = new float[nAlleles];
        float[] gtProbs = new float[3];
        for (int j=0; j<this.nSamples; ++j) {
            setProbs(gv, marker, j, gtProbs, alleleFreq);
            int call = maxIndex(gtProbs, random);
            float exp = (gtProbs[1] + 2*gtProbs[2]);
            float expSquare = (gtProbs[1] + 4*gtProbs[2]);
            sumCall += call;
            sumSquareCall += call*call;
            sumExpected += exp;
            sumExpectedSquare += expSquare;
            sumSquareExpected += (exp*exp);
            sumCallExpected += (call*exp);
        }
        divideBySum(alleleFreq);
    }

    private static void setProbs(GenotypeValues gv, int marker, int sample,
            float[] gtProbs, float[] alleleFreq) {
        Arrays.fill(gtProbs, 0.0f);
        int nGt = gv.marker(marker).nGenotypes();
        float sum = 0.0f;
        for (int gt=0; gt<nGt; ++gt) {
            sum += gv.value(marker, sample, gt);
        }

        int gt = -1;
        for (int a2=0; a2<alleleFreq.length; ++a2) {
            for (int a1=0; a1<a2; ++a1) {
                float gprob = gv.value(marker, sample, ++gt)/sum;
                alleleFreq[a1] += gprob;
                alleleFreq[a2] += gprob;
                gtProbs[ (a1==0) ? 1 : 2 ] += gprob;
            }
            float gprob = gv.value(marker, sample, ++gt)/sum;
            alleleFreq[a2] += 2*gprob;
            gtProbs[ (a2==0) ? 0 : 2] += gprob;
        }
    }

    private static int maxIndex(float[] fa, Random random) {
        int pivot = random.nextInt(fa.length);
        int maxIndex = pivot;
        for (int j=pivot+1; j<fa.length; ++j) {
            if (fa[j]>fa[maxIndex]) {
                maxIndex = j;
            }
        }
        for (int j=0; j<pivot; ++j) {
            if (fa[j]>fa[maxIndex]) {
                maxIndex = j;
            }
        }
        return maxIndex;
    }

    private static void divideBySum(float[] fa) {
        float sum = 0.0f;
        for (float f : fa) {
            sum += f;
        }
        for (int j=0; j<fa.length; ++j) {
            fa[j] /= sum;
        }
    }

    /**
     * Returns the marker.
     * @return the marker.
     */
    public Marker marker() {
        return marker;
    }

    /**
     * Returns an array of length {@code this.marker().nAlleles()} whose
     * {@code j}-th element is the sample frequency of allele {@code j}.
     * @return the allele frequencies.
     */
    public float[] alleleFreq() {
        return alleleFreq.clone();
    }

    /**
     * Returns the estimated squared correlation between the most likely
     * allele dose and the true allele dose.
     * Returns 0.0f if the marker is monomorphic or if most likely allele
     * dose is monomorphic.
     *
     * @return the estimated squared correlation between the most likely
     * allele dose and the true allele dose.
     */
    public float allelicR2() {
        float f = 1.0f /  nSamples;
        float cov = sumCallExpected - (sumCall * sumExpected * f);
        float varBest = sumSquareCall - (sumCall * sumCall * f);
        float varExp = sumExpectedSquare - (sumExpected * sumExpected * f);
        float den = varBest * varExp;
        return (den==0.0f) ? 0.0f : Math.abs( (cov*cov) / den );
    }

    /**
     * Returns the estimated squared correlation between the estimated
     * allele dose and the true allele dose.  Returns 0.0f if the
     * marker is monomorphic.
     *
     * @return the estimated squared correlation between the estimated
     * allele dose and the true allele dose.
     */
    public float doseR2() {
        float f = 1.0f / (float) nSamples;
        float num = sumSquareExpected - (sumExpected * sumExpected * f);
        float den = sumExpectedSquare - (sumExpected * sumExpected * f);
        return (den==0.0f) ? 0.0f : Math.abs(num / den);
    }

    /**
     * Returns the estimated squared correlation between the estimated
     * allele dose and the true allele dose when the variance of
     * the true allele dose is calculated from the estimated
     * allele frequency. Returns 0.0f if the marker is monomorphic.
     *
     * @return the estimated squared correlation between the estimated
     * allele dose and the true allele dose
     */
    public float hweDoseR2() {
        float f = 1.0f / nSamples;
        float num = (sumSquareExpected - (sumExpected*sumExpected*f))/nSamples;
        float altFreq = sumExpected / (2.0f * nSamples);
        float den = 2.0f * altFreq * (1.0f - altFreq);
        return (den==0.0f) ? 0.0f : Math.abs(num / den);
    }

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}.
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(80);
        sb.append(marker);
        sb.append(Const.tab);
        for (int j=0; j<alleleFreq.length; ++j) {
            sb.append( (j==0) ? "AF=" : Const.comma);
            sb.append(alleleFreq[j]);
        }
        sb.append(Const.tab);
        sb.append("AR2=");
        sb.append(format(allelicR2()));
        sb.append(Const.tab);
        sb.append("DR2=");
        sb.append(format(doseR2()));
        sb.append(Const.tab);
        sb.append("HDR2=");
        sb.append(format(hweDoseR2()));
        return sb.toString();
    }

    private static String format(float d) {
        if (Double.isNaN(d)) {
            return "NaN";
        }
        else if (d==Double.POSITIVE_INFINITY) {
            return "Infinity";
        }
        else if (d==Double.NEGATIVE_INFINITY) {
            return "-Infinity";
        }
        else {
            return df.format(d);
        }
    }
}
