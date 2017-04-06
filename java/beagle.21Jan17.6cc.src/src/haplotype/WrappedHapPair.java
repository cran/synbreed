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
 * Class {@code WrappedHapPair} is a {@code HapPair} instance
 * that wraps a {@code SampleHapPairs} object.

* @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class WrappedHapPair implements HapPair {

    private final SampleHapPairs haps;
    private final int hapPair;

    /**
     * Creates a {@code WrappedHapPair} instance representing
     * the specified haplotype pair.
     * @param sampleHapPairs the {@code SampleHapPairs} object that
     * will be "wrapped" by {@code this}
     * @param hapPair a haplotype pair index
     * @throws IllegalArgumentException if
     * {@code hapPair < 0 || hapPair >= sampleHapPairs.nHapPairs()}
     * @throws NullPointerException if {@code sampleHapPairs == null}
     */
    public WrappedHapPair(SampleHapPairs sampleHapPairs, int hapPair) {
        if (hapPair < 0 || hapPair >= sampleHapPairs.nHapPairs()) {
            throw new IllegalArgumentException("hapPair: " + hapPair);
        }
        this.haps = sampleHapPairs;
        this.hapPair = hapPair;
    }

    @Override
    public int allele1(int marker) {
        return haps.allele1(marker, hapPair);
    }

    @Override
    public int allele2(int marker) {
        return haps.allele2(marker, hapPair);
    }

    @Override
    public Markers markers() {
        return haps.markers();
    }

    @Override
    public Marker marker(int marker) {
        return haps.marker(marker);
    }

    @Override
    public int nMarkers() {
        return haps.nMarkers();
    }

    @Override
    public Samples samples() {
        return haps.samples();
    }

    @Override
    public int sampleIndex() {
        return hapPair;
    }
}
