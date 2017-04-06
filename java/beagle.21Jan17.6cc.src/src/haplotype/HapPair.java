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
package haplotype;

import beagleutil.Samples;
import vcf.Marker;
import vcf.Markers;

/**
 * <p>Interface {@code HapPair} represents a pair of haplotypes for a sample.
 * The pair of haplotypes are guaranteed to have non-missing alleles at each
 * marker.
 * </p>
 * All instances of {@code HapPair} are required to be immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface HapPair {

    /**
     * Returns the first allele for the specified marker.
     * @param marker a marker index
     * @return the first allele for the specified marker
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    int allele1(int marker);

    /**
     * Returns the second allele for the specified marker.
     * @param marker a marker index
     * @return the second allele for the specified marker
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    int allele2(int marker);

    /**
     * Returns the markers.
     * @return the markers
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
     * Returns the number of markers.
     * @return the number of markers
     */
    int nMarkers();

    /**
     * Returns the list of samples containing the sample associated with
     * this haplotype pair.
     * @return the list of samples containing the sample associated with
     * this haplotype pair
     */
    Samples samples();

    /**
     * Returns the index of the sample associated with this haplotype pair
     * in the list of samples returned by {@code this.samples()}.
     * @return the index of the sample associated with this haplotype pair
     * in the list of samples returned by {@code this.samples()}
     */
    int sampleIndex();
}
