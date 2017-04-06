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

/**
 * <p>Interface {@code SampleHapPairs} represents a list of samples and a
 * haplotype pair for each sample. Each haplotype pair is guaranteed
 * to have two non-missing alleles at each marker.
 * </p>
 * <p>All instances of {@code SampleHapPairs} are required to be immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface SampleHapPairs extends HapPairs {

    /**
     * Returns the samples.  The {@code k}-th sample corresponds to
     * the {@code k}-th haplotype pair.
     * @return the samples
     */
    Samples samples();

    /**
     * Returns the number of samples.
     * @return the number of samples
     */
    int nSamples();

    /**
     * Returns the number of marker alleles.
     * @param marker a marker index
     * @return the number of marker alleles.
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    int nAlleles(int marker);

    /**
     * Returns {@code true} if this object stores the indices of haplotypes
     * that carry non-major alleles, and returns {@code false} otherwise.
     * @param marker a marker index
     * @return {@code true} if this object stores the indices of haplotypes
     * that carry non-major alleles
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    boolean storesNonMajorIndices(int marker);

    /**
     * Returns the index of the major allele.
     * @param marker a marker index
     * @return the index of the major allele
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws UnsupportedOperationException if
     * {@code storesNonMajorIndices(marker) == false}
     */
    int majorAllele(int marker);

    /**
     * Returns the number of haplotypes that carry the specified allele.
     * @param marker a marker index
     * @param allele an allele index
     * @return the number of haplotypes that carry the specified allele
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IllegalArgumentException if
     * {@code allele == this.majorAllele()}
     * @throws IndexOutOfBoundsException if
     * {@code allele < 0 ||  allele >= this.nAlleles()}
     * @throws UnsupportedOperationException if
     * {@code storesNonMajorIndices(marker) == false}
     */
    int alleleCount(int marker, int allele);

    /**
     * Returns index of the haplotype that carries the specified copy of the
     * specified allele.
     * @param marker a marker index
     * @param allele an allele index
     * @param copy a copy index.
     * @return index of the haplotype that carries the specified allele.
     * @throws IllegalArgumentException if
     * {@code allele == this.majorAllele()}
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code allele < 0 ||  allele >= this.nAlleles()}
     * @throws IndexOutOfBoundsException if
     * {@code copy < 0 ||  copy >= this.alleleCount(allele)}
     * @throws UnsupportedOperationException if
     * {@code storesNonMajorIndices(marker) == false}
     */
    int hapIndex(int marker, int allele, int copy);
}
