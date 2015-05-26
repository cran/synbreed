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
 * <p>Interface {@code VcfEmission} represents genotype emission
 * probabilities for a set of samples at a single marker.
 * </p>
 * <p>All instances of {@code VcfEmission} are required to be immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface VcfEmission {

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
     * Returns the marker.
     * @return the marker.
     */
    Marker marker();

    /**
     * Returns {@code true} if the genotype emission probabilities
     * for each sample are determined by a phased called genotype
     * that has no missing alleles, and returns {@code false} otherwise.
     * @return {@code true} if the genotype emission probabilities
     * for each sample are determined by a phased called genotype
     * that has no missing alleles, and {@code false} otherwise.
     */
    boolean isRefData();

    /**
     * Returns {@code true} if the genotype emission probabilities
     * for each sample are determined by a genotype that contains
     * only missing alleles, and returns {@code false} otherwise.
     * @return {@code true}  if the genotype emission probabilities
     * for each sample are determined by a genotype that contains
     * only missing alleles, and {@code false} otherwise.
     */
    boolean isMissingData();

    /**
     * Returns the probability of the observed data if the specified pair
     * of ordered alleles is the true genotype in the specified sample.
     * @param sample the sample index.
     * @param allele1 the first allele index.
     * @param allele2 the second allele index
     * @return the probability of the observed data if the specified pair
     * of ordered alleles is present at the specified marker in the specified
     * sample.
     *
     * @throws IndexOutOfBoundsException if
     * {@code samples< 0 || samples>=this.nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code allele1< 0 || allele1>=this.marker().nAlleles()}
     * @throws IndexOutOfBoundsException if
     * {@code allele2< 0 || allele2>=this.marker().nAlleles()}
     */
    float gl(int sample, byte allele1, byte allele2);

    /**
     * Returns {@code true} if the genotype emission probabilities
     * for the specified sample are determined by a phased called genotype,
     * and returns {@code false} otherwise.
     * @param sample the sample index.
     * @return {@code true} if the genotype emission probabilities
     * for the specified sample are determined by a phased called genotype,
     * and {@code false} otherwise.
     *
     * @throws IndexOutOfBoundsException if
     * {@code sample<0 || sample>=this.nSamples()}
     */
    boolean isPhased(int sample);

    /**
     * Returns the first allele if the genotype emission probabilities
     * for the specified sample are determined by a called genotype, and
     * returns -1 otherwise.  Alleles are arbitrarily ordered
     * if the genotype is unphased.
     * @param sample the sample index.
     * @return the first allele if the genotype emission probabilities
     * for the specified sample are determined by a called genotype, and
     * -1 otherwise.
     *
     * @throws IndexOutOfBoundsException if
     * {@code sample<0 || samples>=this.nSamples()}
     */
    byte allele1(int sample);

    /**
     * Returns the second allele if the genotype emission probabilities
     * for the specified sample are determined by a called genotype, and
     * returns -1 otherwise.  Alleles are arbitrarily ordered
     * if the genotype is unphased.
     * @param sample the sample index.
     * @return the second allele if the genotype emission probabilities
     * for the specified sample are determined by a called genotype, and
     * -1 otherwise.
     *
     * @throws IndexOutOfBoundsException if
     * {@code sample<0 || samples>=this.nSamples()}
     */
    byte allele2(int sample);

    /**
     * Returns a VCF record corresponding to {@code this}.  The returned
     * VCF record will have missing QUAL and INFO fields and will have "PASS"
     * in the filter field.
     * @return a VCF record corresponding to {@code this}.
     */
    @Override
    String toString();
}
