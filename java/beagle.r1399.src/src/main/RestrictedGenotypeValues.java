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
 * <p>Class {@code RestrictedGenotypeValues} is a wrapper for a
 * {@code GenotypeValues} instance.  The wrapper restricts the list of
 * markers to a sublist of the markers in the wrapped object.
 * </p>
 * Instances of class {@code RestrictedtGenotypeValues} are thread-safe.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class RestrictedGenotypeValues implements GenotypeValues {

    /*
     * All instances of {@code GenotypeValues} are required to be
     * thread-safe.
     */
    private final GenotypeValues gv;
    private final Markers restriction;
    private final int[] inclusionMap; // restrictedMarkers -> gv.markers()

    /**
     * Constructs a new {@code RestrictedGenotypeValues} instance.
     * @param gv genotype values whose markers will be restricted.
     * @param restriction the list of restricted markers..
     * @throws IllegalArgumentException if the list of markers
     * represented by {@code restrictedMarkers} is not a sublist
     * of the list of markers represented by {@code gv.markers()}.
     * @throws NullPointerException if
     * {@code gv==null || restrictedMarkers==null}.
     */
    public RestrictedGenotypeValues(GenotypeValues gv, Markers restriction) {
        this.gv = gv;
        this.restriction = restriction;

        int[] map = new int[restriction.nMarkers()];
        int index = 0;
        for (int j=0, n=gv.nMarkers(); j<n; ++j) {
            if (index < map.length
                    && gv.marker(j).equals(restriction.marker(index))) {
                map[index++] = j;
            }
        }
        if (index != map.length) {
            throw new IllegalArgumentException("invalid restricted markers");
        }
        this.inclusionMap = map;
    }

    @Override
    public float value(int marker, int sample, int genotype) {
        return gv.value(inclusionMap[marker], sample, genotype);
    }

    @Override
    public void add(int sample, double[] values) {
        if (values.length != restriction.sumGenotypes()) {
            throw new IllegalArgumentException("values.length=" + values.length);
        }
        int index = 0;
        for (int j=0; j<inclusionMap.length; ++j) {
            int nGt = restriction.marker(j).nGenotypes();
            for (int gt=0; gt<nGt; ++gt) {
                gv.add(inclusionMap[j], sample, gt, values[index++]);
            }
        }
        assert index==values.length;
    }

    @Override
    public void add(int marker, int sample, int genotype, double value) {
        gv.add(inclusionMap[marker], sample, genotype, value);
    }

    @Override
    public Samples samples() {
        return gv.samples();
    }

    @Override
    public int nSamples() {
        return gv.samples().nSamples();
    }

    @Override
    public Marker marker(int marker) {
        return restriction.marker(marker);
    }

    @Override
    public Markers markers() {
        return restriction;
    }

    @Override
    public int nMarkers() {
        return restriction.nMarkers();
    }
}
