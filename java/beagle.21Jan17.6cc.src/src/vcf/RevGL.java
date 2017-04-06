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
 * <p>Class {@code RevGL} is wrapper for a {@code GL} instance.  The wrapper
 * reverses the order of markers in the wrapped object.
 * </p>
 * <p>Instances of class {@code RevGL} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class RevGL implements GL {

    /*
     * All instances of the {@code GL} interface are required to be immutable.
     */
    private final GL gl;
    private final int lastMarker;

    /**
     * Constructs a new {@code RevGL} instance.
     * @param gl genotype emission probabilities that will be
     * wrapped by the new instance
     * @throws NullPointerException if {@code gl == null}
     */
    public RevGL(GL gl) {
        this.gl = gl;
        this.lastMarker = gl.nMarkers() - 1;
    }

    @Override
    public boolean isRefData() {
        return gl.isRefData();
    }

    @Override
    public float gl(int marker, int sample, int allele1, int allele2) {
        int revMarker = lastMarker - marker;
        return gl.gl(revMarker, sample, allele1, allele2);
    }

    @Override
    public boolean isPhased(int marker, int sample) {
        int revMarker = lastMarker - marker;
        return gl.isPhased(revMarker, sample);
    }

    @Override
    public int allele1(int marker, int sample) {
        int revMarker = lastMarker - marker;
        return gl.allele1(revMarker, sample);
    }

    @Override
    public int allele2(int marker, int sample) {
        int revMarker = lastMarker - marker;
        return gl.allele2(revMarker, sample);
    }

    @Override
    public int allele(int marker, int hap) {
        int revMarker = lastMarker - marker;
        return gl.allele(revMarker, hap);
    }

    @Override
    public Marker marker(int marker) {
        int revMarker = lastMarker - marker;
        return gl.marker(revMarker);
    }

    @Override
    public Markers markers() {
        return gl.markers().reverse();
    }

    @Override
    public int nMarkers() {
        return gl.nMarkers();
    }

    @Override
    public int nHaps() {
        return gl.nHaps();
    }

    @Override
    public int nSamples() {
        return gl.nSamples();
    }

    @Override
    public Samples samples() {
        return gl.samples();
    }

    @Override
    public String toString() {
        return gl.toString();
    }
}
