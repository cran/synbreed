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
 * <p>Class {@code FixedGenotypeValues} is a wrapper for a {@code GenotypeValues}
 * instance.  The wrapper fixes the values of the genotype values at a
 * subset of markers so that the {@code add()} methods do not alter
 * the genotype values at any marker in the subset.
 * </p>
 * Instances of class {@code FixedGenotypeValues} are thread-safe.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class FixedGenotypeValues implements GenotypeValues {

    /*
     * All instances of {@code GenotypeValues} are required to be
     * thread-safe.
     */
    private final GenotypeValues gv;
    private final boolean[] isFixed;

    /**
     * Constructs a {@code FixedGenotypeValues} instance.
     * @param gv genotype values that will be wrapped by the new instance.
     * @param fixedMarkers a list of markers whose genotype values will be
     * fixed.
     * @throws NullPointerException if {@code gv==null || fixedMarkers==null}.
     */
    public FixedGenotypeValues(GenotypeValues gv,
            Markers fixedMarkers) {
        this.gv = gv;
        this.isFixed = new boolean[gv.nMarkers()];
        for (int j=0; j<isFixed.length; ++j) {
            if (fixedMarkers.contains(gv.marker(j))) {
                this.isFixed[j] = true;
            }
        }
    }

    @Override
    public float value(int marker, int sample, int genotype) {
        return gv.value(marker, sample, genotype);
    }

    @Override
    public void add(int sample, double[] values) {
        int index = 0;
        for (int j=0, n=gv.nMarkers(); j<n; ++j) {
            int nGt = gv.marker(j).nGenotypes();
            for (int gt=0; gt<nGt; ++gt) {
                this.add(j, sample, gt, values[index++]);
            }
        }
        assert index==values.length;
    }

    @Override
    public void add(int marker, int sample, int genotype, double value) {
        if (isFixed[marker]==false) {
            gv.add(marker, sample, genotype, value);
        }
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
        return gv.marker(marker);
    }

    @Override
    public Markers markers() {
        return gv.markers();
    }

    @Override
    public int nMarkers() {
        return gv.nMarkers();
    }
}
