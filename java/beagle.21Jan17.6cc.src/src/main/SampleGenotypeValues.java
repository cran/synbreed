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

import beagleutil.Samples;
import vcf.Markers;
import vcf.Marker;

/**
 * <p>Class {@code SampleGenotypeValues} stores a value for each possible
 * genotype at each marker for one sample.
 * </p>
 * <p>Class {@code SampleGenotypeValues} is thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class SampleGenotypeValues {

    private final Markers markers;
    private final Samples samples;
    private final int sampleIndex;
    private final float[] gtValues;

    /**
     * Constructs a {@code SampleGenotypeValues} instance for the
     * specified markers and sample with initial value 0 for each possible
     * genotype at each marker.
     * @param markers the list of markers
     * @param samples the list of samples
     * @param sampleIndex a sample index
     * @throws IllegalArgumentException if
     * {@code sampleIndex < 0  || sampleIndex >= sampes.nSamples()}
     * @throws NullPointerException if
     * {@code markers == null || samples == null}
     */
    public SampleGenotypeValues(Markers markers, Samples samples, int sampleIndex) {
        if (sampleIndex < 0 || sampleIndex >= samples.nSamples()) {
            throw new IllegalArgumentException(String.valueOf(sampleIndex));
        }
        this.markers = markers;
        this.samples = samples;
        this.sampleIndex = sampleIndex;
        this.gtValues = new float[markers.sumGenotypes()];
    }

    /**
     * Returns the specified genotype value.
     *
     * @param marker a marker index
     * @param genotype a genotype index
     * @return the specified genotype value
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code genotype < 0 || genotype >= this.marker(marker).nGenotypes()}
     */
    public synchronized float value(int marker, int genotype) {
        checkGenotype(marker, genotype);
        return gtValues[markers.sumGenotypes(marker) + genotype];
    }

    /**
     * Adds the specified genotype values to {@code this}.  This method is
     * equivalent to
     * <pre>
     * for (int m=0; m&lt;this.nMarkers(); ++m) {
     *     offset = this.markers().sumGenotypes(m);
     *     for (int gt=0; gt&lt;this.marker(m).nGenotypes(); ++gt) {
     *         this.add(marker, gt, values[offset + gt])
     *     }
     * }
     * </pre>
     *
     * @param values an array with {@code this.markers.sumGenotypes()}
     * elements containing the genotype values to be added
     * @throws IllegalArgumentException if
     * {@code values.length != this.markers().sumGenotypes()}
     * @throws NullPointerException if {@code values == null}
     */
    public synchronized void add(double[] values) {
        if (values.length != gtValues.length) {
            String s = "values.length=" + values.length;
            throw new IllegalArgumentException(s);
        }
        for (int j=0; j<values.length; ++j) {
            gtValues[j] += values[j];
        }
    }

    /**
     * Adds the specified genotype value to {@code this}.
     * @param marker a marker index
     * @param genotype a genotype index
     * @param value the value to be added
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code genotype < 0 || genotype >= this.marker(marker).nGenotypes()}
     */
    public synchronized void add(int marker, int genotype, double value) {
        checkGenotype(marker, genotype);
        gtValues[markers.sumGenotypes(marker) + genotype] += value;
    }

    private void checkGenotype(int marker, int genotype) {
        int nGenotypes = markers.marker(marker).nGenotypes();
        if (genotype < 0 || genotype >= nGenotypes) {
            throw new IndexOutOfBoundsException("genotype: " + genotype);
        }
    }

    /**
     * Returns the list of markers.
     * @return the list of markers
     */
    public Markers markers() {
        return markers;
    }

    /**
     * Returns the list of samples.
     * @return the list of samples
     */
    public Samples samples() {
        return samples;
    }

    /**
     * Returns the sample index.
     * @return the sample index
     */
    public int sampleIndex() {
        return sampleIndex;
    }

    /**
     * Returns the specified marker.
     * @param marker a marker index
     * @return the specified marker
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    public Marker marker(int marker) {
        return markers.marker(marker);
    }

    /**
     * Returns the number of markers.
     * @return the number of markers
     */
    public int nMarkers() {
        return markers.nMarkers();
    }

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        return this.getClass().toString();
    }
}
