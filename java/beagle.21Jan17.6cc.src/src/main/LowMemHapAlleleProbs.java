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
package main;

import beagleutil.Samples;
import vcf.Marker;
import vcf.Markers;

/**
 * <p>Class {@code LowMemHapAlleleProbs} stores allele probabilities for
 * a haplotype.
 * </p>
 * Instances of class {@code LowMemHapAlleleProbs} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class LowMemHapAlleleProbs implements HapAlleleProbs {

    private static final int N_BINS = 256;
    private static final int SHIFT = 128;
    private static final float INCREMENT = 1f/N_BINS;

    private final Markers markers;
    private final Samples samples;
    private final int hap;
    private final byte[] alleleBin;

    /**
     * Constructs a new {@code LowMemHapAlleleProbs} instance.  The
     * {@code alleleProbs} array lists the probability of each allele for
     * each marker, sorted first by marker index, and then by allele index.
     * @param markers the markers
     * @param samples the samples
     * @param hap the haplotype index
     * @param alleleProbs the allele probabilities
     * @throws IllegalArgumentException if
     * {@code alleleProbs.length != markers.sumAlleles()}
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap >= 2*samples.nSamples()}
     * @throws NullPointerException if
     * {@code markers == null || samples == null}
     */
    public LowMemHapAlleleProbs(Markers markers, Samples samples, int hap,
            float[] alleleProbs) {
        if (alleleProbs.length != markers.sumAlleles()) {
            throw new IllegalArgumentException("inconsistent data");
        }
        if (hap < 0 || hap >= 2*samples.nSamples()) {
            throw new IndexOutOfBoundsException(String.valueOf(hap));
        }
        byte[] nonRefProbs = new byte[alleleProbs.length - markers.nMarkers()];
        int index = 0;
        for (int m=0, n=markers.nMarkers(); m<n; ++m) {
            float sum = 0f;
            int start = markers.sumAlleles(m);
            int end = markers.sumAlleles(m+1);
            for (int j=start; j<end; ++j) {
                float p = alleleProbs[j];
                if (p < 0 || p > 1.01f || Float.isNaN(p)) {
                    throw new IllegalArgumentException(String.valueOf(p));
                }
                sum += alleleProbs[j];
            }
            if (sum > 1.01f) {
                throw new IllegalArgumentException(String.valueOf(sum));
            }
            for (int j=start; j<end - 1; ++j) {
                nonRefProbs[index++] = convertToByte(alleleProbs[j] / sum);
            }
        }
        assert index == nonRefProbs.length;

        this.markers = markers;
        this.samples = samples;
        this.hap = hap;
        this.alleleBin = nonRefProbs;
    }

    private static byte convertToByte(float f) {
        if (f >= 1f) {
            f = 0.99999f;
        }
        int bin = ((int) Math.floor(f*N_BINS)) - SHIFT;
        return (byte) bin;
    }

    private static float convertToFloat(byte b) {
        return (b + 128.5f)*INCREMENT;
    }

    @Override
    public float allele(int marker, int allele) {
        int nAlleles = markers.marker(marker).nAlleles();
        if (allele < 0 || allele >= nAlleles) {
            throw new IllegalArgumentException(String.valueOf(allele));
        }
        int start = markers.sumAlleles(marker) - marker;
        if (nAlleles == 2) {
            float f = convertToFloat(alleleBin[start]);
            return (allele == 0) ? f : 1f - f;
        }
        else if (allele == nAlleles - 1) {
            return lastAlleleProb(marker);
        }
        else {
            return convertToFloat(alleleBin[start + allele]);
        }
    }

    private float lastAlleleProb(int marker) {
        int nAlleles = markers.marker(marker).nAlleles();
        int start = markers.sumAlleles(marker) - marker;
        int end = start + nAlleles - 1;
        float sum = 1f;
        for (int j = start; j < end; ++j) {
            sum -= convertToFloat(alleleBin[j]);
        }
        return (sum < 0f) ? 0f : sum;
    }

    @Override
    public int alleleWithMaxProb(int marker) {
        int nAlleles = markers.marker(marker).nAlleles();
        int start = markers.sumAlleles(marker) - marker;
        if (nAlleles == 2) {
            return alleleBin[start] >= 0 ? 0 : 1;
        }
        else {
            int bestIndex = start;
            int end = start + nAlleles - 1;
            float sumProb = 0f;
            for (int j = start; j<end; ++j) {
                sumProb += convertToFloat(alleleBin[j]);
                if (alleleBin[j] > alleleBin[bestIndex]) {
                    bestIndex = j;
                }
            }
            if ( sumProb < 0.5f) {
                return nAlleles - 1;
            }
            else {
                return bestIndex - start;
            }
        }
    }

    @Override
    public int nMarkers() {
        return markers.nMarkers();
    }

    @Override
    public Markers markers() {
        return markers;
    }

    @Override
    public Marker marker(int marker) {
        return markers.marker(marker);
    }

    @Override
    public int hapIndex() {
        return hap;
    }

    @Override
    public Samples samples() {
        return samples;
    }

    /**
     * Returns a string representation of {@code this}. The exact details
     * of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(100);
        sb.append('[');
        sb.append(this.getClass().toString());
        sb.append(": hap=");
        sb.append(hap);
        sb.append(']');
        return sb.toString();
    }
}
