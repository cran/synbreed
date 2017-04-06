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
 * <p>Interface {@code AL} (Allele Likelihoods) represents allele
 * likelihoods for a set of haplotypes.
 * </p>
 * Instances of class {@code AL} are required to be immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface AL {

    /**
     * Returns the probability of the observed data if the specified allele
     * is the true allele for the specified marker and haplotype.
     * @param marker a marker index
     * @param haplotype a haplotype index
     * @param allele a allele index
     * @return the probability of the observed data if the specified allele
     * is the true allele for the specified marker and haplotype
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code haplotype < 0 || haplotype >= this.nHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code allele < 0 || allele >= this.marker(marker).nAlleles()}
     */
    float al(int marker, int haplotype, int allele);

    /**
     * Returns the allele on the specified haplotype if the allele
     * emission probabilities are determined by a called allele, and
     * returns -1 otherwise.
     * @param marker a marker index
     * @param haplotype a haplotype index
     * @return the allele on the specified haplotype if the allele
     * emission probabilities are determined by a called allele, and
     * returns -1 otherwise
     *
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap  >= this.nHaps()}
     */
    int allele(int marker, int haplotype);

    /**
     * Returns the number of markers.
     * @return the number of markers
     */
    int nMarkers();

    /**
     * Returns the specified marker.
     * @param marker the marker index
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
     * Returns the number of samples.
     * @return the number of samples
     */
    int nSamples();

    /**
     * Returns the list of samples.
     * @return the list of samples
     */
    Samples samples();

    /**
     * Returns the number of haplotypes.
     * @return the number of haplotypes
     */
    int nHaps();

    /**
     * Returns the allelic error probability.
     * @return the allelic error probability
     */
    float errProb();

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    String toString();
}
