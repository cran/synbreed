/*
 * Copyright (C) 2015 Brian L. Browning
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
package sample;

import java.util.Arrays;
import main.HapAlleleProbs;
import main.LowMemHapAlleleProbs;
import vcf.Markers;

/**
 * <p>Class {@code LSHapBaum} implements the Baum hidden Markov model
 * forward and backward algorithms for imputing missing alleles on a
 * target haplotype.
 * </p>
 * <p>Instances of class {@code LSHapBaum} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class LSHapBaum {

    private final ImputationData impData;
    private final boolean lowMem;
    private final int n;    // number of reference haplotypes
    private final Markers refMarkers;
    private final float[] alleleProbs;
    private final float[][] fwdVal;
    private final float[] bwdVal;
    private final int[] fwdValueIndex2Marker;

    private final RefHapSegs refHapSegs;
    private final float[][] fwdHapProbs;
    private final float[][] bwdHapProbs;

    private int windowIndex = -9999;
    private int arrayIndex = -9999;

    /**
     * Creates a {@code LSHapBaum} instance from the specified data.
     *
     * @param impData the input data for genotype imputation
     * @param lowMem {@code true} if a low-memory checkpoint algorithm
     * should be used, and {@code false} otherwise
     *
     * @throws NullPointerException if {@code impData == null}
     */
    public LSHapBaum(ImputationData impData, boolean lowMem) {
        this.impData = impData;
        this.lowMem = lowMem;
        this.n = impData.refHapPairs().nHaps();
        this.refMarkers = impData.refHapPairs().markers();
        this.alleleProbs = new float[refMarkers.sumAlleles()];

        int nClusters = impData.nClusters();
        int size = lowMem ? (int) Math.ceil(Math.sqrt(1 + 8*nClusters)/2.0) + 1
                : nClusters;
        this.fwdValueIndex2Marker = new int[size];
        this.fwdVal = new float[size][n];
        this.bwdVal = new float[n];

        this.refHapSegs = impData.refHapSegs();
        this.fwdHapProbs = new float[impData.nClusters()][];
        this.bwdHapProbs = new float[impData.nClusters()][];
        for (int j=0; j < nClusters; ++j) {
            this.fwdHapProbs[j] = new float[refHapSegs.nSeq(j+1)];
            this.bwdHapProbs[j] = new float[refHapSegs.nSeq(j)];
        }
    }

    /**
     * <p>Estimates and returns allele probabilities for the specified target
     * haplotype. Estimated allele probabilities are conditional on the hidden
     * Markov model (HMM) and the input data represented by
     * {@code this.imputationData()}.
     * </p>
     *
     * @param hap a target data haplotype index
     * @return allele probabilities for the specified target haplotype
     *
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap >= this.imputationData().targetHapPairs().nHaps()}
     */
    public HapAlleleProbs randomHapSample(int hap) {
        Arrays.fill(alleleProbs, 0f);
        int nMarkers = impData.nClusters();
        windowIndex = 0;
        arrayIndex = -1;
        setForwardValues(0, nMarkers, hap);
        Arrays.fill(bwdVal, 1.0f/n);
        setStateProbs(nMarkers-1, currentIndex());
        for (int m=nMarkers-2; m>=0; --m) {
            setBwdValue(m, hap);
            setStateProbs(m, previousIndex(hap));
        }
        setAlleleProbs(alleleProbs);
        return new LowMemHapAlleleProbs(refMarkers, impData.targetSamples(),
                hap, alleleProbs);
    }

    /**
     * Returns the input data for genotype imputation.
     * @return the input data for genotype imputation
     */
    public ImputationData imputationData() {
        return impData;
    }

    private void setForwardValues(int start, int end, int hap) {
        float lastSum = 1.0f;
        for (int m=start; m<end; ++m) {
            float probRec = impData.pRecomb(m);
            float probNoRec = 1.0f - probRec;
            float noErrProb = impData.noErrProb(m);
            float errProb = impData.errProb(m);
            float shift = probRec/n;
            float scale = probNoRec/lastSum;
            int prev = currentIndex();
            int next = nextIndex();
            float sum = 0.0f;
            fwdValueIndex2Marker[next] = m;
            int a = impData.targetAllele(m, hap);
            for (int h=0; h<n; ++h) {
                float em = (a == impData.refAllele(m, h)) ? noErrProb : errProb;
                fwdVal[next][h] = m==0 ? em : em*(scale*fwdVal[prev][h] + shift);
                sum += fwdVal[next][h];
            }
            lastSum = sum;
        }
    }

    private void setBwdValue(int m, int hap) {
        int mP1 = m + 1;
        float probRec = impData.pRecomb(mP1);
        float probNoRec = 1.0f - probRec;
        float noErrProb = impData.noErrProb(mP1);
        float errProb = impData.errProb(mP1);
        float sum = 0f;
        int al = impData.targetAllele(mP1, hap);
        for (int h=0; h<n; ++h) {
            float em = (al==impData.refAllele(mP1, h)? noErrProb : errProb);
            bwdVal[h] *= em;
            sum += bwdVal[h];
        }
        float scale = probNoRec/sum;
        float shift = probRec/n;
        for (int h=0; h<n; ++h) {
            bwdVal[h] = scale*bwdVal[h] + shift;
        }
    }

    private void setStateProbs(int m, int fwdIndex) {
        Arrays.fill(fwdHapProbs[m], 0f);
        Arrays.fill(bwdHapProbs[m], 0f);
        for (int h=0; h<n; ++h) {
            float stateProbs = fwdVal[fwdIndex][h]*bwdVal[h];
            fwdHapProbs[m][refHapSegs.seq(m+1, h)] += stateProbs;
            bwdHapProbs[m][refHapSegs.seq(m, h)] += stateProbs;
        }
        float sum = sum(fwdHapProbs[m]);
        scale(fwdHapProbs[m], sum);
        scale(bwdHapProbs[m], sum);
    }

    private static float sum(float[] fa) {
        float sum = 0f;
        for (float f : fa) {
            sum += f;
        }
        return sum;
    }

    private static void scale(float[] fa, float divisor) {
        for (int j=0; j<fa.length; ++j) {
            fa[j] /= divisor;
        }
    }

    private static float threshold(int nSeq) {
        return Math.min(0.005f, 1.0f/nSeq);
    }

    private void setAlleleProbs(float[] alleleProbs) {
        setFirstAlleleProbs(alleleProbs);
        int nSegsM1 = refHapSegs.nSegs() - 1;
        for (int j=1; j<nSegsM1; ++j) {
            setAlleleProbs(alleleProbs, j);
        }
        setLastAlleleProbs(alleleProbs);
    }

    private void setFirstAlleleProbs(float[] alleleProbs) {
        int segment = 0;
        int nSeq = refHapSegs.nSeq(segment);
        int endRefMarker = refHapSegs.segStart(segment + 1);
        float threshold = threshold(nSeq);
        for (int seq=0; seq<nSeq; ++seq) {
            if (bwdHapProbs[segment][seq] >= threshold) {
                for (int m=0; m<endRefMarker; ++m) {
                    int start = refMarkers.sumAlleles(m);
                    int allele = refHapSegs.allele(segment, m, seq);
                    alleleProbs[start + allele] += bwdHapProbs[segment][seq];
                }
            }
        }
    }

    private void setAlleleProbs(float[] alleleProbs, int segment) {
        assert segment > 0;
        int clustStart = refHapSegs.segStart(segment);
        int clustEnd = refHapSegs.segEnd(segment - 1);
        int nextClustStart = refHapSegs.segStart(segment + 1);
        int nSeq = refHapSegs.nSeq(segment);
        float threshold = threshold(nSeq);
        for (int seq=0; seq<nSeq; ++seq) {
            boolean useFwd = fwdHapProbs[segment-1][seq] >= threshold;
            boolean useBwd = bwdHapProbs[segment][seq] >= threshold;
            if (useFwd) {
                for (int m=clustStart; m<clustEnd; ++m) {
                    int start = refMarkers.sumAlleles(m);
                    int allele = refHapSegs.allele(segment, m - clustStart, seq);
                    alleleProbs[start + allele] += fwdHapProbs[segment-1][seq];
                }
            }
            if (useFwd || useBwd) {
                for (int m=clustEnd; m<nextClustStart; ++m) {
                    int start = refMarkers.sumAlleles(m);
                    int allele = refHapSegs.allele(segment, m - clustStart, seq);
                    double wt = impData.weight(m);
                    alleleProbs[start + allele] += wt*fwdHapProbs[segment-1][seq];
                    alleleProbs[start + allele] += (1-wt)*bwdHapProbs[segment][seq];
                }
            }
        }
    }

    private void setLastAlleleProbs(float[] alleleProbs) {
        int segment = refHapSegs.nSegs() - 1;
        int cluster = segment - 1;
        int refMarkerStart = refHapSegs.segStart(segment);
        int refMarkerEnd = refHapSegs.segEnd(segment);
        int nSeq = refHapSegs.nSeq(segment);
        float threshold = threshold(nSeq);
        for (int seq=0; seq<nSeq; ++seq) {
            if (fwdHapProbs[cluster][seq] >= threshold) {
                for (int m=refMarkerStart; m<refMarkerEnd; ++m) {
                    int start = refMarkers.sumAlleles(m);
                    int allele = refHapSegs.allele(segment, m - refMarkerStart, seq);
                    alleleProbs[start + allele] += fwdHapProbs[cluster][seq];
                }
            }
        }
    }

    private int nextIndex() {
        ++arrayIndex;
        if (arrayIndex == fwdVal.length) {
            ++windowIndex;
            arrayIndex = windowIndex;
        }
        return arrayIndex;
    }

    private int currentIndex() {
        return arrayIndex;
    }

    private int previousIndex(int hap) {
        if (arrayIndex == windowIndex) {
            --windowIndex;
            arrayIndex = windowIndex;
            int start = fwdValueIndex2Marker[arrayIndex] + 1;
            int end = start + ( fwdVal.length - (arrayIndex + 1) );
            setForwardValues(start, end, hap);
            return arrayIndex;
        }
        else {
            return --arrayIndex;
        }
    }
}
