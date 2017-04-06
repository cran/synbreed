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
import java.util.concurrent.atomic.AtomicReferenceArray;

/**
 * <p>Class {@code BasicGenotypeValues} stores values for each possible
 * genotype for each sample at each marker.
 * </p>
 * <p>Instances of class {@code BasicGenotypeValues} are thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class BasicGenotypeValues implements GenotypeValues {

    private final Markers markers;
    private final Samples samples;

    /*
     * Class {@code SampleGenotypeValues} is thread-safe.
     */
    private final AtomicReferenceArray<SampleGenotypeValues> values;

    /**
     * Constructs a new {@code BasicGenotypeValues} instance with initial
     * value 0 for each possible genotype for each sample at each marker.
     * @param markers a list of markers
     * @param samples a list of samples
     * @throws NullPointerException if
     * {@code markers == null || samples == null}
     */
    public BasicGenotypeValues(Markers markers, Samples samples) {
        if (markers==null) {
            throw new NullPointerException("markers");
        }
        if (samples==null) {
            throw new NullPointerException("samples");
        }
        this.markers = markers;
        this.samples = samples;
        this.values = new AtomicReferenceArray<>(samples.nSamples());
        for (int j=0, n=samples.nSamples(); j<n; ++j) {
            this.values.set(j, new SampleGenotypeValues(markers, samples, j));
        }
    }

    @Override
    public float value(int marker, int sample, int genotype) {
        return this.values.get(sample).value(marker, genotype);
    }

    @Override
    public void add(int sample, double[] values) {
        this.values.get(sample).add(values);
    }

    @Override
    public void add(int marker, int sample, int genotype, double value) {
        this.values.get(sample).add(marker, genotype, value);
    }

    @Override
    public Samples samples() {
        return samples;
    }

    @Override
    public int nSamples() {
        return samples.nSamples();
    }

    @Override
    public Marker marker(int marker) {
        return markers.marker(marker);
    }

    @Override
    public Markers markers() {
        return markers;
    }

    @Override
    public int nMarkers() {
        return markers.nMarkers();
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
