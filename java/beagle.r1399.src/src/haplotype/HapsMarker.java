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

import vcf.Marker;

/**
 * <p>Interface {@code HapsMarker} represents data at marker from a
 * {@code HapsInterface} object.
 * </p>
 * All instances of {@code HapsMarkers} are required to be
 * immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface HapsMarker {

     /**
     * Returns the allele for the specified haplotype.
     * @param haplotype a haplotype index.
     * @return the allele for the specified haplotype.
     *
     * @throws IndexOutOfBoundsException if
     * {@code haplotype<0 || haplotype>=this.nHaps()}.
     */
    byte allele(int haplotype);

    /**
     * Returns the first allele for the specified haplotype pair.
     * @param hapPair a haplotype pair index.
     * @return the first allele for the specified haplotype pair.
     *
     * @throws IndexOutOfBoundsException if
     * {@code hapPair<0 || hapPair>=this.nHapPairs()}.
     */
    byte allele1(int hapPair);

    /**
     * Returns the second allele for the specified haplotype pair.
     * @param hapPair a haplotype pair index.
     * @return the second allele for the specified haplotype pair.
     *
     * @throws IndexOutOfBoundsException if
     * {@code hapPair<0 || hapPair>=this.nHapPairs()}.
     */
    byte allele2(int hapPair);

    /**
     * Returns the marker.
     * @return the marker.
     */
    Marker marker();

    /**
     * Returns the number of haplotypes.  The returned value is equal to
     * {@code 2*this.nHapPairs()}.
     * @return the number of haplotypes.
     */
    int nHaps();

    /**
     * Returns the number of haplotype pairs.  The returned value is
     * equal to {@code this.nHaps()/2}.
     * @return the number of haplotype pairs.
     */
    int nHapPairs();

    /**
     * Returns the sample identifier index corresponding to the specified
     * haplotype pair.
     * @param hapPair a haplotype pair index.
     * @return the sample identifier index corresponding to the specified
     * haplotype pair.
     * @throws IndexOutOfBoundsException if
     * {@code hapPair<0 || hapPair>=this.nHapPairs()}
     * @see beagleutil.SampleIds
     */
    int idIndex(int hapPair);
}
