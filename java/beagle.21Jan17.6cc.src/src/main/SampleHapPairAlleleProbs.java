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
import haplotype.SampleHapPairs;
import vcf.Marker;
import vcf.Markers;

/**
 * <p>Class {@code SampleHapPairAlleleProbs} is a wrapper for a
 * {@code SampleHapPairs} instance.
 * </p>
 * <p>Instances of class {@code HaplotypeAlleleProbs} are immutable.
 * </p>
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class SampleHapPairAlleleProbs implements AlleleProbs {

    private final SampleHapPairs sampleHapPairs;

    /**
     * Constructs a new {@code SampleHapPairAlleleProbs} instance that wraps
     * the specified {@code SampleHapPairs} object.  The alleles in
     * the specified {@code SampleHapPairs} instance will have
     * probability 1.
     *
     * @param sampleHapPairs the sample haplotype pairs that will
     * be wrapped by {@code this}
     *
     * @throws NullPointerException if {@code sampleHapPairs == null}
     */
    public SampleHapPairAlleleProbs(SampleHapPairs sampleHapPairs) {
        if (sampleHapPairs==null) {
            throw new NullPointerException("sampleHapPairs==null");
        }
        this.sampleHapPairs = sampleHapPairs;
    }

    @Override
    public float alProb1(int marker, int sample, int allele) {
        return allele==sampleHapPairs.allele1(marker, sample) ? 1f : 0f;
    }

    @Override
    public float alProb2(int marker, int sample, int allele) {
        return allele==sampleHapPairs.allele2(marker, sample) ? 1f : 0f;
    }

    @Override
    public float gtProb(int marker, int sample, int allele1, int allele2) {
        return alProb1(marker, sample, allele1)*alProb2(marker, sample, allele2);
    }

    @Override
    public int allele1(int marker, int sample) {
        return sampleHapPairs.allele1(marker, sample);
    }

    @Override
    public int allele2(int marker, int sample) {
        return sampleHapPairs.allele2(marker, sample);
    }

    @Override
    public Marker marker(int marker) {
        return sampleHapPairs.marker(marker);
    }

    @Override
    public Markers markers() {
        return sampleHapPairs.markers();
    }

    @Override
    public int nMarkers() {
        return sampleHapPairs.nMarkers();
    }

    @Override
    public int nSamples() {
        return sampleHapPairs.nSamples();
    }

    @Override
    public Samples samples() {
        return sampleHapPairs.samples();
    }
}
