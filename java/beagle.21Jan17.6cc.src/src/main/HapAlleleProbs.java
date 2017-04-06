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
import vcf.Markers;
import vcf.Marker;

/**
 * <p>Interface {@code HapAlleleProbs} stores allele probabilities for
 * a haplotype.
 * </p>
 * All instances of {@code HapAlleleProbs} are required to be immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface HapAlleleProbs {

    /**
     * Returns the specified allele probability.
     *
     * @param marker a marker index
     * @param allele an allele index
     * @return the specified allele probability
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code allele < 0 || allele >= this.marker(marker).nAlleles()}
     */
    float allele(int marker, int allele);

    /**
     * Returns the allele with maximum posterior probability.  If more than
     * one allele has maximum posterior probability, one of the
     * alleles with maximum posterior probability will be returned.
     * @param marker a marker index
     * @return the allele with maximum posterior probability
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    int alleleWithMaxProb(int marker);

    /**
     * Returns the number of markers.
     * @return the number of markers
     */
    int nMarkers();

    /**
     * Returns the list of markers.
     * @return the list of markers
     */
    Markers markers();

    /**
     * Returns the specified marker.
     * @param marker a marker index
     * @return the specified marker
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    Marker marker(int marker);

    /**
     * Returns the haplotype index. The two haplotypes for sample {@code k}
     * have indices {@code 2*k} and {@code 2*k + 1}.
     * @return the haplotype index
     */
    int hapIndex();

    /**
     * Returns the list of samples.
     * @return the list of samples
     */
    Samples samples();
}
