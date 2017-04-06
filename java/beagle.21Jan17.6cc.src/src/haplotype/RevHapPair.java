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
 * <p>Class {@code RevHapPair} is a wrapper for a {@code HapPair}
 * instance.  The wrapper reverses the order of markers in the wrapped object.
 * </p>
 * <p>Instances of class {@code RevHapPair} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class RevHapPair implements HapPair {

    /*
     * All instances of the {@code HapPair} interface are required to be
     * immutable.
     */
    private final HapPair hapPair;
    private final int lastMarker;

    /**
     * Creates a new {@code RevHapPair} instance from the specified data.
     * @param hapPair the haplotype pair that will be wrapped by the
     * new instance
     * @throws NullPointerException if {@code hapPair == null}
     */
    public RevHapPair(HapPair hapPair) {
        this.hapPair = hapPair;
        this.lastMarker = hapPair.nMarkers() - 1;
    }

    @Override
    public int allele1(int marker) {
        return hapPair.allele1(lastMarker - marker);
    }

    @Override
    public int allele2(int marker) {
        return hapPair.allele2(lastMarker - marker);
    }

    @Override
    public Markers markers() {
        return hapPair.markers().reverse();
    }

    @Override
    public Marker marker(int marker) {
        return hapPair.marker(lastMarker - marker);
    }

    @Override
    public int nMarkers() {
        return hapPair.nMarkers();
    }

    @Override
    public Samples samples() {
        return hapPair.samples();
    }

    @Override
    public int sampleIndex() {
        return hapPair.sampleIndex();
    }
}
