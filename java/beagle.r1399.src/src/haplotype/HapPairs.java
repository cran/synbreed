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
import vcf.Markers;

/**
 * <p>Class {@code HaplotypesInterface} stores an list of haplotype pairs.
 * </p>
 * All instances of {@code HaplotypesInterface} are required to
 * be immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface HapPairs {

     /**
     * Returns the allele for the specified marker and haplotype.
     * @param marker a marker index.
     * @param haplotype a haplotype index.
     * @return the allele for the specified marker and haplotype.
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker<0 || marker>=this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code haplotype<0 || haplotype>=this.nHaps()}
     */
    byte allele(int marker, int haplotype);

    /**
     * Returns the first allele for the specified marker and haplotype pair.
     * @param marker a marker index.
     * @param hapPair a haplotype pair index.
     * @return the first allele for the specified marker and haplotype pair.
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker<0 || marker>=this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code hapPair< 0 || hapPair>=this.nHapPairs()}
     */
    byte allele1(int marker, int hapPair);

    /**
     * Returns the second allele for the specified marker and haplotype pair.
     * @param marker a marker index.
     * @param hapPair a haplotype pair index.
     * @return the second allele for the specified marker and haplotype pair.
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker<0 || marker>=this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code hapPair< 0 || hapPair>=this.nHapPairs()}
     */
    byte allele2(int marker, int hapPair);

    /**
     * Returns the number of markers. The returned value is equal to
     * {@code this.markers().nMarkers()}.
     * @return the number of markers.
     */
    int nMarkers();

    /**
     * Returns the markers in each haplotype.
     * @return the markers in each haplotype
     */
    Markers markers();

    /**
     * Returns the specified marker.
     * @param marker a marker index.
     * @return the specified marker.
     * @throws IndexOutOfBoundsException if
     * {@code marker<0 || marker>=this.nMarkers()}
     */
    Marker marker(int marker);

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
