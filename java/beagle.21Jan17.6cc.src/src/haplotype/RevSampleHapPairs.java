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
package haplotype;

import beagleutil.Samples;
import vcf.Marker;
import vcf.Markers;

/**
 * <p>Class {@code RevSampleHapPairs} is a wrapper for a {@code SampleHapPairs}
 * instance.  The wrapper reverses the order of markers in the wrapped object.
 * </p>
 * <p>Instances of class {@code RevSampleHapPairs} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class RevSampleHapPairs implements SampleHapPairs {

    /*
     * All instances of the {@code SampleHapPairs} interface are required to be
     * immutable.
     */
    private final SampleHapPairs hapPairs;
    private final int lastMarker;

    /**
     * Creates a new {@code RevSampleHapPairs} instance from the specified data.
     * @param hapPairs the sample haplotype pairs that will be wrapped by the
     * new instance
     * @throws NullPointerException if {@code hapPairs == null}
     */
    public RevSampleHapPairs(SampleHapPairs hapPairs) {
        this.hapPairs = hapPairs;
        this.lastMarker = hapPairs.nMarkers() - 1;
    }

    @Override
    public int allele1(int marker, int hapPair) {
        return hapPairs.allele1(lastMarker - marker, hapPair);
    }

    @Override
    public int allele2(int marker, int hapPair) {
        return hapPairs.allele2(lastMarker - marker, hapPair);
    }

    @Override
    public int allele(int marker, int haplotype) {
        return hapPairs.allele(lastMarker - marker, haplotype);
    }

    @Override
    public int nMarkers() {
        return hapPairs.nMarkers();
    }

    @Override
    public Markers markers() {
        return hapPairs.markers().reverse();
    }

    @Override
    public Marker marker(int marker) {
        return hapPairs.marker(lastMarker - marker);
    }

    @Override
    public int nHaps() {
        return hapPairs.nHaps();
    }

    @Override
    public int nHapPairs() {
        return hapPairs.nHapPairs();
    }

    @Override
    public int nSamples() {
        return hapPairs.nSamples();
    }

    @Override
    public Samples samples() {
        return hapPairs.samples();
    }

    @Override
    public Samples samples(int hapPair) {
        return hapPairs.samples(hapPair);
    }

    @Override
    public int sampleIndex(int hapPair) {
        return hapPairs.sampleIndex(hapPair);
    }

    @Override
    public int nAlleles(int marker) {
        return hapPairs.nAlleles(lastMarker - marker);
    }

    @Override
    public boolean storesNonMajorIndices(int marker) {
        return hapPairs.storesNonMajorIndices(lastMarker - marker);
    }

    @Override
    public int majorAllele(int marker) {
        return hapPairs.majorAllele(lastMarker - marker);
    }

    @Override
    public int alleleCount(int marker, int allele) {
        return hapPairs.alleleCount(lastMarker - marker, allele);
    }

    @Override
    public int hapIndex(int marker, int allele, int copy) {
        return hapPairs.hapIndex(lastMarker - marker, allele, copy);
    }
}
