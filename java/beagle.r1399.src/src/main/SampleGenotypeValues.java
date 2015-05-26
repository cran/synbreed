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
import java.util.Arrays;

/**
 * <p>Class {@code SampleGenotypeValues} stores a value for each possible
 * genotype for one sample.
 * </p>
 * Class {@code SampleGenotypeValues} is thread-safe.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class SampleGenotypeValues {

    private final int sample;
    private final Markers markers;
    private final float[] gtValues;

    /**
     * Constructs a {@code SampleGenotypeValues} instance for the
     * specified markers and sample with initial value 0 for each genotype.
     * @param markers the list of markers.
     * @param sample a sample index.
     * @throws IllegalArgumentException if {@code sample < 0}.
     * @throws NullPointerException if {@code markers==null}.
     */
    public SampleGenotypeValues(Markers markers, int sample) {
        if (sample < 0) {
            throw new IllegalArgumentException("sample<0: " + sample);
        }
        this.sample = sample;
        this.markers = markers;
        this.gtValues = new float[markers.sumGenotypes()];
    }

    /**
     * Returns a new {@code SampleGenotypeValues} instance obtained
     * by restricting the data of {@code this} to the specified
     * marker interval.
     * @param start the starting marker index (inclusive).
     * @param end the ending marker index (exclusive).
     * @return a new {@code SampleGenotypeValues} instance obtained
     * by restricting the data of {@code this} to the specified
     * marker interval.
     *
     * @throws IndexOutOfBoundsException if
     * {@code start<0 || end>=this.nMarkers()}.
     * @throws IllegalArgumentException if {@code start>=end}.
     */
    public synchronized SampleGenotypeValues restrict(int start, int end) {
        return new SampleGenotypeValues(this, start, end);
    }

    private SampleGenotypeValues(SampleGenotypeValues gv, int start, int end) {
        this.sample = gv.sample;
        this.markers = gv.markers.restrict(start, end);
        int from = gv.markers.sumGenotypes(start);
        int to = gv.markers.sumGenotypes(end);
        this.gtValues = Arrays.copyOfRange(gv.gtValues, from, to);
    }

    /**
     * Returns the specified genotype value.
     *
     * @param marker a marker index.
     * @param genotype a genotype index.
     * @return the specified genotype value.
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker<0 || marker>=this.nMarkers()}.
     * @throws IndexOutOfBoundsException if
     * {@code genotype<0 || genotype>=this.marker(marker).nGenotypes()}.
     */
    public synchronized float value(int marker, int genotype) {
        checkGenotype(marker, genotype);
        return gtValues[markers.sumGenotypes(marker) + genotype];
    }

    /**
     * Adds the specified genotype values to {@code this}.  This method is
     * equivalent to
     * <pre>
     * for (m=0; m&lt;this.nMarkers(); ++m) {
     *     offset = this.markers().sumGenotypes(m);
     *     for (gt=0; gt&lt;this.marker(m).nGenotypes(); ++gt) {
     *         this.add(marker, gt, values[offset + gt])
     *     }
     * }
     * </pre>
     *
     * @param values an array of length {@code this.markers.sumGenotypes()}
     * containing the genotype values to be added.     * @throws IllegalArgumentException if
     * {@code values.length!=this.markers().sumGenotypes()}.
     * @throws NullPointerException if {@code values==null}.
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
     * @param marker a marker index.
     * @param genotype a genotype index.
     * @param value the value to be added.
     * @throws IndexOutOfBoundsException if
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker<0 || marker>=this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code genotype<0 || genotype>=this.marker(marker).nGenotypes()}.
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
     * Returns the sample index.
     * @return the sample index.
     */
    public int sample() {
        return sample;
    }

    /**
     * Returns the list of markers.
     * @return the list of markers.
     */
    public Markers markers() {
        return markers;
    }

    /**
     * Returns the specified marker.
     * @param marker a marker index.
     * @return the specified marker.
     * @throws IndexOutOfBoundsException if
     * {@code marker<0 || marker>=this.nMarkers()}.
     */
    public Marker marker(int marker) {
        return markers.marker(marker);
    }

    /**
     * Returns the number of markers.
     * @return the number of markers.
     */
    public int nMarkers() {
        return markers.nMarkers();
    }
}
