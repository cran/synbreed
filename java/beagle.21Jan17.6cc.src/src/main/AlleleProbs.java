/*
 * Copyright (C) 2015 Brian L. Browning
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
import vcf.Marker;
import vcf.Markers;

/**
 * <p>Interface {@code AlleleProbs} represents per-haplotype allele
 * probabilities for a list of samples.
 * </p>
 * <p>All instances of {@code AlleleProbs} are required to be immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface AlleleProbs {

    /**
     * Returns the probability that the specified marker allele is
     * present on the first haplotype of the specified sample.
     *
     * @param marker a marker index
     * @param sample a sample index
     * @param allele an allele index
     * @return the probability that the specified marker allele is
     * present on the first haplotype of the specified sample
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code allele < 0 || allele >= this.marker(marker).nAlleles()}
     */
    float alProb1(int marker, int sample, int allele);

    /**
     * Returns the probability that the specified marker allele is
     * present on the second haplotype of the specified sample.
     *
     * @param marker a marker index
     * @param sample a sample index
     * @param allele an allele index
     * @return the probability that the specified marker allele is
     * present on the second haplotype of the specified sample.
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code allele < 0 || allele >= this.marker(marker).nAlleles()}
     */
    float alProb2(int marker, int sample, int allele);

    /**
     * Returns the phased genotype probability, equal to
     * {@code (this.allele1(marker, sample, allele1)
     *        * this.allele2(marker, sample, allele2))}.
     *
     * @param marker a marker index
     * @param sample a sample index
     * @param allele1 allele index of the allele on the first haplotype
     * @param allele2 allele index of the allele on the second haplotype
     * @return the phased genotype probability equal to
     *       {@code (this.allele1(marker, sample, allele1)
     *              * this.allele2(marker, sample, allele2))}
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code allele1 < 0 || allele1 >= this.marker(marker).nAlleles()}
     * @throws IndexOutOfBoundsException if
     * {@code allele2 < 0 || allele2 >= this.marker(marker).nAlleles()}
     */
    float gtProb(int marker, int sample, int allele1, int allele2);

    /**
     * Returns the marker allele with maximum probability for the
     * first haplotype of the specified sample. If more than one allele
     * has maximum probability, one of the alleles with maximum
     * probability will be returned.
     * @param marker a marker index
     * @param sample a sample index
     * @return the marker allele with maximum probability for the
     * first haplotype of the specified sample
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     */
    int allele1(int marker, int sample);

    /**
     * Returns the marker allele with maximum probability for the
     * second haplotype of the specified sample. If more than one allele
     * has maximum probability, one of the alleles with maximum
     * probability will be returned.
     * @param marker a marker index
     * @param sample a sample index
     * @return the marker allele with maximum probability for the
     * second haplotype of the specified sample
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     */
    int allele2(int marker, int sample);

    /**
     * Returns the specified marker.
     * @param marker a marker index
     * @return the specified marker
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    Marker marker(int marker);

    /**
     * Returns the list of markers.
     * @return the list of markers
     */
    Markers markers();

    /**
     * Returns the number of markers.
     * @return the number of markers
     */
    int nMarkers();

    /**
     * Returns the number of samples.
     * @return the number of samples
     */
    int nSamples();

    /**
     * Returns the list of samples.
     * @return the list of samples
     */
    Samples samples();
}
