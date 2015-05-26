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

/**
 * <p>Class {@code RevAL} is wrapper for an {@code AL} instance.
 * The wrapper reverses the order of markers in the wrapped object.
 * </p>
 * <p>Instances of class {@code RevAL} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class RevAL implements AL {

    private final int lastMarker;
    private final AL al;

    /**
     * Constructs a new {@code RevAL} instance.
     * @param al allele emission probabilities that will be wrapped by the
     * new instance.
     * @throws NullPointerException if {@code al==null}.
     */
    public RevAL(AL al) {
        this.lastMarker = al.nMarkers() - 1;
        this.al = al;
    }

    @Override
    public float al(int marker, int haplotype, byte allele) {
        int revMarker = lastMarker - marker;
        return al.al(revMarker, haplotype, allele);
    }


    @Override
    public Marker marker(int marker) {
        int revMarker = lastMarker - marker;
        return al.marker(revMarker);
    }

    @Override
    public Markers markers() {
        return al.markers().reverse();
    }

    @Override
    public int nMarkers() {
        return al.nMarkers();
    }

    @Override
    public int nSamples() {
        return al.nSamples();
    }

    @Override
    public Samples samples() {
        return al.samples();
    }

    @Override
    public String toString() {
        return al.toString();
    }
}
