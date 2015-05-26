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
 * <p>Interface {@code GL} (Genotype Likelihoods) represents genotype
 * emission probabilities for a set of samples.
 * </p>
 * <p>All instances of {@code GL} are required to be immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface GL {

    /**
     * Returns {@code true} if the observed data for each sample
     * and for each marker are phased genotypes that have no missing alleles,
     * and returns {@code false} otherwise.
     * @return {@code true} if the observed data for each sample
     * and for each marker are phased genotypes that have no missing alleles,
     * and returns {@code false} otherwise.
     */
    boolean isRefData();

    /**
     * Returns the probability of the observed data if the specified pair
     * of ordered alleles is the true ordered genotype at the specified marker
     * in the specified sample.
     * @param marker the marker index.
     * @param sample the sample index.
     * @param allele1 the first allele index.
     * @param allele2 the second allele index
     * @return the probability of the observed data if the specified pair
     * of ordered alleles is the true ordered genotype at the specified marker
     * in the specified sample.
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker<0 || marker>=this.nMarkers()}.
     * @throws IndexOutOfBoundsException if
     * {@code samples<0 || samples>=this.nSamples()}.
     * @throws IndexOutOfBoundsException if
     * {@code allele1<0 || allele1>=this.marker(marker).nAlleles()}.
     * @throws IndexOutOfBoundsException if
     * {@code allele2<0 || allele2>=this.marker(marker).nAlleles()}.
     */
    float gl(int marker, int sample, byte allele1, byte allele2);

    /**
     * Returns the first allele for the specified marker and sample
     * if the observed data is a genotype, and returns -1 otherwise.
     * Alleles are arbitrarily ordered if the genotype is unphased.
     * @param marker the marker index.
     * @param sample the sample index.
     * @return the first allele for the specified marker and sample
     * if the observed data for the sample is a genotype, and returns
     * -1 otherwise.
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker<0 || marker>=this.nMarkers()}.
     * @throws IndexOutOfBoundsException if
     * {@code sample<0 || samples>=this.nSamples()}.
     */
    byte allele1(int marker, int sample);

    /**
     * Returns the second allele for the specified marker and sample
     * if the observed data is a genotype, and returns -1 otherwise.
     * Alleles are arbitrarily ordered if the genotype is unphased.
     * @param marker the marker index.
     * @param sample the sample index.
     * @return the second allele for the specified marker and sample
     * if the observed data for the sample is a genotype, and returns
     * -1 otherwise.
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker<0 || marker>=this.nMarkers()}.
     * @throws IndexOutOfBoundsException if
     * {@code sample<0 || samples>=this.nSamples()}.
     */
    byte allele2(int marker, int sample);

    /**
     * Returns the number of markers.
     * @return the number of markers.
     */
    int nMarkers();

    /**
     * Returns the specified marker.
     * @param marker the marker index.
     * @return the specified marker.
     * @throws IndexOutOfBoundsException if
     * {@code marker<0 || marker>=this.nMarkers()}.
     */
    Marker marker(int marker);

    /**
     * Returns the list of markers.
     * @return the list of markers.
     */
    Markers markers();

    /**
     * Returns the number of samples.
     * @return the number of samples.
     */
    int nSamples();

   /**
     * Returns the list of samples.
     * @return the list of samples.
     */
    Samples samples();

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}.
     */
    @Override
    String toString();
}
