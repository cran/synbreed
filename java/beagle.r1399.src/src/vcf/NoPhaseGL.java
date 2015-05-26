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
 * <p>Class {@code NoPhaseGL} is a wrapper for a {@code GL}
 * instance.  The wrapper's {@code gl()} method ignores genotype phase
 * information in the wrapped object.
 * </p>
 * <p>Instances of class {@code NoPhaseGL} are immutable.
 * </p>
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class NoPhaseGL implements GL {

    private final GL gl;

    /**
     * Constructs a new {@code NoPhaseGL} instance.
     * @param gl genotype emission probabilities that will be wrapped by
     * the new instance.
     * @throws NullPointerException if {@code gl==null}.
     */
    public NoPhaseGL(GL gl) {
        if (gl==null) {
            throw new NullPointerException("gl==null");
        }
        this.gl = gl;
    }

    @Override
    public float gl(int marker, int sample, byte a1, byte a2) {
        if (a1==a2) {
            return gl.gl(marker, sample, a1, a2);
        }
        else {
            float f1 = gl.gl(marker, sample, a1, a2);
            float f2 = gl.gl(marker, sample, a2, a1);
            return Math.max(f1, f2);
        }
    }

    @Override
    public boolean isRefData() {
        return gl.isRefData();
    }

    @Override
    public byte allele1(int marker, int sample) {
        return gl.allele1(marker, sample);
    }

    @Override
    public byte allele2(int marker, int sample) {
        return gl.allele2(marker, sample);
    }

    @Override
    public int nMarkers() {
        return gl.nMarkers();
    }

    @Override
    public Marker marker(int marker) {
        return gl.marker(marker);
    }

    @Override
    public Markers markers() {
        return gl.markers();
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
