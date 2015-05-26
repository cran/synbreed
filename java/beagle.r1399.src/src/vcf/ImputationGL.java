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
package vcf;

import beagleutil.Samples;
import haplotype.SampleHapPairs;

/**
 * <p>Class {@code ImputationGL} represents genotype emission
 * probabilities for a list of samples that have missing genotypes
 * at a subset of markers, and phased, non-missing genotypes
 * at all other markers.
 * </p>
 * <p>Instances of class {@code ImputationGL} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class ImputationGL implements GL {

    private final Markers refMarkers;
    private final int[] nonRefIndex;
    private final SampleHapPairs haps;

    private static int[] nonRefIndex(Markers refMarkers, Markers nonRefMarkers) {
        int n = refMarkers.nMarkers();
        int[] nonRefIndex = new int[n];
        int index = 0;
        for (int j=0; j<n; ++j) {
            if (index < nonRefMarkers.nMarkers()
                    && refMarkers.marker(j).equals(nonRefMarkers.marker(index))) {
                nonRefIndex[j] = index++;
            }
            else {
                nonRefIndex[j] = -1;
            }
        }
        if (index != nonRefMarkers.nMarkers()) {
            String s = "inconsistent markers";
            throw new IllegalArgumentException(s);
        }
        return nonRefIndex;
    }

    /**
     * Constructs an {@code ImputationGL} instance for the specified data.
     * @param refMarkers the reference markers.
     * @param haps haplotype pairs for the target samples.  The list of
     * markers returned by the {@code haps.markers()} method must be
     * a sublist of {@code refMarkers}.
     *
     * @throws IllegalArgumentException if the list of markers returned
     * by the {@code haps.markers()} method is not a sublist of
     * {@code refMarkers}.
     * @throws NullPointerException if
     * {@code refMarkers==null || haps==null}.
     */
    public ImputationGL(Markers refMarkers, SampleHapPairs haps) {
        if (refMarkers==null) {
            throw new NullPointerException("refMarkers==null");
        }
        this.refMarkers = refMarkers;
        this.nonRefIndex = nonRefIndex(refMarkers, haps.markers());
        this.haps = haps;
    }

    @Override
    public boolean isRefData() {
        return haps.nMarkers()==refMarkers.nMarkers();
    }

    @Override
    public float gl(int marker, int sample, byte allele1, byte allele2) {
        int nAlleles = refMarkers.marker(marker).nAlleles();
        if (allele1 < 0 || allele1 >= nAlleles) {
            String s = "marker=" + marker + " allele1: " + allele1;
            throw new IllegalArgumentException(s);
        }
        if (allele2 < 0 || allele2 >= nAlleles) {
            String s = "marker=" + marker + " allele2: " + allele2;
            throw new IllegalArgumentException(s);
        }
        if (nonRefIndex[marker] != -1) {
            byte a1 = haps.allele1(nonRefIndex[marker], sample);
            byte a2 = haps.allele2(nonRefIndex[marker], sample);
            return (a1==allele1 && a2==allele2) ? 1.0f : 0.0f;
        }
        else {
            return 1.0f;
        }
    }

    @Override
    public byte allele1(int marker, int sample) {
        if (nonRefIndex[marker] != -1) {
            return haps.allele1(nonRefIndex[marker], sample);
        }
        else {
            return -1;
        }
    }

    @Override
    public byte allele2(int marker, int sample) {
        if (nonRefIndex[marker] != -1) {
            return haps.allele2(nonRefIndex[marker], sample);
        }
        else {
            return -1;
        }
    }

    @Override
    public int nMarkers() {
        return refMarkers.nMarkers();
    }

    @Override
    public Marker marker(int marker) {
        return refMarkers.marker(marker);
    }

    @Override
    public Markers markers() {
        return refMarkers;
    }

    @Override
    public int nSamples() {
        return haps.nSamples();
    }

    @Override
    public Samples samples() {
        return haps.samples();
    }

    @Override
    public String toString() {
        StringBuilder sb  = new StringBuilder();
        sb.append("[ImputationGL: nMarkers=");
        sb.append(nMarkers());
        sb.append(" nSamples=");
        sb.append(nSamples());
        sb.append(']');
        return sb.toString();
    }
}
