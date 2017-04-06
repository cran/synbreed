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
package main;

import vcf.Markers;
import vcf.Marker;
import beagleutil.Samples;

/**
 * <p>Class {@code RevGenotypeValues} is a wrapper for a {@code GenotypeValues}
 * instance.  The wrapper reverses the order of markers in the wrapped object.
 * </p>
 * Instances of class {@code RevGenotypeValues} are thread-safe.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class RevGenotypeValues implements GenotypeValues {

    /*
     * All instances of the {@code GenotypeValues} interface are required to be
     * thread-safe.
     */
    private final GenotypeValues gv;

    /**
     * Constructs a new {@code RevGenotypeValues} instance from the specified
     * data.
     * @param gv genotype values that will be wrapped by the new instance
     * @throws NullPointerException if {@code gv == null}
     */
    public RevGenotypeValues(GenotypeValues gv) {
        this.gv = gv;
    }

    @Override
    public float value(int marker, int sample, int genotype) {
        int revMarker = gv.nMarkers() - 1 - marker;
        return gv.value(revMarker, sample, genotype);
    }

    @Override
    public void add(int sample, double[] values) {
        if (values.length != gv.markers().sumGenotypes()) {
            throw new IllegalArgumentException("values.length=" + values.length);
        }
        int index = 0;
        for (int m=0, n=gv.nMarkers(); m<n; ++m) {
            int revMarker = gv.nMarkers() - 1 - m;
            int nGt = gv.marker(revMarker).nGenotypes();
            for (int gt=0; gt<nGt; ++gt) {
                gv.add(revMarker, sample, gt, values[index++]);
            }
        }
    }

    @Override
    public void add(int marker, int sample, int genotype, double value) {
        int revMarker = gv.nMarkers() - 1 - marker;
        gv.add(revMarker, sample, genotype, value);
    }

    @Override
    public Samples samples() {
        return gv.samples();
    }

    @Override
    public int nSamples() {
        return gv.nSamples();
    }

    @Override
    public Marker marker(int marker) {
        return gv.markers().reverse().marker(marker);
    }

    @Override
    public Markers markers() {
        return gv.markers().reverse();
    }

    @Override
    public int nMarkers() {
        return gv.nMarkers();
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(100);
        sb.append('[');
        sb.append(this.getClass().toString());
        sb.append(']');
        return sb.toString();
    }
}
