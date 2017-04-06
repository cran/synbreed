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

import vcf.Markers;
import vcf.Marker;
import beagleutil.Samples;
import haplotype.SampleHapPairs;

/**
 * <p>Class {@code ConstrainedAlleleProbs} is a wrapper for an
 * {@code AlleleProbs} instance that changes the wrapped haplotype allele
 * probabilities for a subset of markers.
 * </p>
 * <p>Instances of class {@code ConstrainedAlleleProbs} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class ConstrainedAlleleProbs implements AlleleProbs {

    private final AlleleProbs alProbs;
    private final SampleHapPairs shp;
    private final int[] indexMap;

    /**
     * Construct a new {@code ConstrainedAlleleProbs} instance.  The alleles
     * in the specified {@code SampleHapPairs} object will
     * have probability 1 in the new {@code ConstrainedAlleleProbs} object.
     * All other allele probabilities in the constructed object
     * are determined by the wrapped {@code AlleleProbs} object.
     * @param shp phased haplotype pairs for a subset of markers
     * @param alProbs the allele probabilities
     * @param indexMap an array of length {@code alProbs.nMarkers()}
     * whose {@code j}-th element is the index of marker
     * {@code alProbs.marker(j)} in {@code shp.markers()}, or is -1 if the
     * marker is not present in {@code shp.markers()}
     *
     * @throws IllegalArgumentException if
     * {@code alProbs.nMarkers() != indexMap.length}
     * @throws IllegalArgumentException if
     * {@code (indexMap[j] != -1
     * && alProbs.marker(j).equals(shp.marker(indexMap[j])) == false)}
     * for any {@code j} satisfying {@code (0 <= j && j < indexMap.length)}
     * @throws NullPointerException if
     * {@code shp == null || alProbs == null || indexMap == null}
     */
    public ConstrainedAlleleProbs(SampleHapPairs shp, AlleleProbs alProbs,
            int[] indexMap) {
        if (alProbs.nMarkers() != indexMap.length) {
            throw new IllegalArgumentException("inconsistent markers");
        }
        for (int j=0; j<indexMap.length; ++j) {
            if (indexMap[j] != -1) {
                if (alProbs.marker(j).equals(shp.marker(indexMap[j]))==false) {
                    throw new IllegalArgumentException("inconsistent markers");
                }
            }
        }
        if (shp.samples().equals(alProbs.samples())==false) {
            throw new IllegalArgumentException("inconsistent sample");
        }
        this.shp = shp;
        this.alProbs = alProbs;
        this.indexMap = indexMap.clone();
    }

    /**
     * Returns {@code true} if the specified marker is not present in the
     * input data and returns {@code false} otherwise.
     * @param marker a marker index
     * @return {@code true} if the specified marker is not present in the
     * input target data
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    public boolean isImputed(int marker) {
        return indexMap[marker] == -1;
    }

    @Override
    public float alProb1(int marker, int sample, int allele) {
        int targetMarker = indexMap[marker];
        if (targetMarker == -1) {
            return alProbs.alProb1(marker, sample, allele);
        }
        else {
            return shp.allele1(targetMarker, sample) == allele ? 1f : 0f;
        }
    }

    @Override
    public float alProb2(int marker, int sample, int allele) {
        int targetMarker = indexMap[marker];
        if (targetMarker == -1) {
            return alProbs.alProb2(marker, sample, allele);
        }
        else {
            return shp.allele2(targetMarker, sample) == allele ? 1f : 0f;
        }
    }

    @Override
    public float gtProb(int marker, int sample, int allele1, int allele2) {
        return alProb1(marker, sample, allele1)*alProb2(marker, sample, allele2);
    }

    @Override
    public int allele1(int marker, int sample) {
        int targetMarker = indexMap[marker];
        if (targetMarker == -1) {
            return alProbs.allele1(marker, sample);
        }
        else {
            return shp.allele1(targetMarker, sample);
        }
    }

    @Override
    public int allele2(int marker, int sample) {
        int targetMarker = indexMap[marker];
        if (targetMarker == -1) {
            return alProbs.allele2(marker, sample);
        }
        else {
            return shp.allele2(targetMarker, sample);
        }
    }

    @Override
    public int nMarkers() {
        return alProbs.nMarkers();
    }

    @Override
    public Markers markers() {
        return alProbs.markers();
    }

    @Override
    public Marker marker(int marker) {
        return alProbs.marker(marker);
    }

    @Override
    public int nSamples() {
        return alProbs.nSamples();
    }

    @Override
    public Samples samples() {
        return alProbs.samples();
    }
}
