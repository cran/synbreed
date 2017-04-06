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
import blbutil.Const;

/**
 * <p>Interface {@code VcfEmission} represents genotype emission
 * probabilities for a set of samples at a single marker.
 * </p>
 * <p>All instances of {@code VcfEmission} are required to be immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface VcfEmission extends HapsMarker {

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
     * Returns {@code true} if the genotype emission probabilities
     * for each sample are determined by a phased called genotype
     * that has no missing alleles, and returns {@code false} otherwise.
     * @return {@code true} if the genotype emission probabilities
     * for each sample are determined by a phased called genotype
     * that has no missing alleles
     */
    boolean isRefData();

    /**
     * Returns the probability of the observed data if the specified pair
     * of ordered alleles is the true genotype in the specified sample.
     * @param sample the sample index
     * @param allele1 the first allele index
     * @param allele2 the second allele index
     * @return the probability of the observed data if the specified pair
     * of ordered alleles is the true genotype in the specified sample
     *
     * @throws IndexOutOfBoundsException if
     * {@code samples < 0 || samples >= this.nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code allele1 < 0 || allele1 >= this.marker().nAlleles()}
     * @throws IndexOutOfBoundsException if
     * {@code allele2 < 0 || allele2 >= this.marker().nAlleles()}
     */
    float gl(int sample, int allele1, int allele2);

    /**
     * Returns {@code true} if the genotype emission probabilities for
     * the specified sample are determined by a phased, nonmissing genotype,
     * and returns {@code false} otherwise.
     * @param sample the sample index
     * @return {@code true} if the genotype emission probabilities
     * for the specified sample are determined by a phased, nonmissing genotype
     *
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()}
     */
    boolean isPhased(int sample);

    /**
     * Returns the first allele for the specified sample or -1 if the
     * allele is missing. Alleles are arbitrarily ordered
     * if the genotype is unphased.
     * @param sample the sample index
     * @return the first allele for the specified sample or -1 if the
     * allele is missing
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || samples >= this.nSamples()}
     */
    @Override
    int allele1(int sample);

    /**
     * Returns the second allele for the specified sample or -1 if the
     * allele is missing. Alleles are arbitrarily ordered
     * if the genotype is unphased.
     * @param sample the sample index
     * @return the second allele for the specified sample or -1 if the
     * allele is missing
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || samples >= this.nSamples()}
     */
    @Override
    int allele2(int sample);

    /**
     * Returns the allele on the specified haplotype or -1 if the
     * allele is missing. Alleles are arbitrarily ordered if the genotype
     * is unphased.
     * @param hap the haplotype index
     * @return the allele on the specified haplotype or -1 if the
     * allele is missing
     *
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap  >= this.nHaps()}
     */
    @Override
    int allele(int hap);

    /**
     * Returns the number of marker alleles.
     * @return the number of marker alleles.
     */
    int nAlleles();

    /**
     * Returns {@code true} if this instance stores the indices of haplotypes
     * that carry non-major alleles, and returns {@code false} otherwise.
     *
     * @return {@code true} if this instance stores the indices of haplotypes
     * that carry non-major alleles, and returns {@code false} otherwise
     */
    boolean storesNonMajorIndices();

    /**
     * Returns the index of the major allele.
     * @return the index of the major allele
     * @throws UnsupportedOperationException if
     * {@code storesNonMajorIndices() == false}
     */
    int majorAllele();

    /**
     * Returns the number of haplotypes that carry the specified allele.
     * @param allele an allele index
     * @return the number of haplotypes that carry the specified allele
     * @throws IllegalArgumentException if
     * {@code allele == this.majorAllele()}
     * @throws IndexOutOfBoundsException if
     * {@code allele < 0 ||  allele >= this.nAlleles()}
     * @throws UnsupportedOperationException if
     * {@code storesNonMajorIndices() == false}
     */
    int alleleCount(int allele);

    /**
     * Returns index of the haplotype that carries the specified copy of the
     * specified allele.
     * @param allele an allele index
     * @param copy a copy index
     * @return index of the haplotype that carries the specified allele
     * @throws IllegalArgumentException if
     * {@code allele == this.majorAllele()}
     * @throws IndexOutOfBoundsException if
     * {@code allele < 0 ||  allele >= this.nAlleles()}
     * @throws IndexOutOfBoundsException if
     * {@code copy < 0 ||  copy >= this.alleleCount(allele)}
     * @throws UnsupportedOperationException if
     * {@code storesNonMajorIndices() == false}
     */
    int hapIndex(int allele, int copy);

    /**
     * Returns a VCF record corresponding to {@code this}.  The returned
     * VCF record will have missing QUAL and INFO fields, will have "PASS"
     * in the filter field, and will have a GT format field.
     * @return a VCF record corresponding to {@code this}
     */
    default String toVcfRec() {
        StringBuilder sb = new StringBuilder(100);
        sb.append(marker());
        sb.append(Const.tab);
        sb.append(Const.MISSING_DATA_CHAR);  // QUAL
        sb.append(Const.tab);
        sb.append("PASS");                   // FILTER
        sb.append(Const.tab);
        sb.append(Const.MISSING_DATA_CHAR);  // INFO
        sb.append(Const.tab);
        sb.append("GT");                     // FORMAT
        for (int j=0, n=nSamples(); j<n; ++j) {
            int a1 = allele1(j);
            int a2 = allele2(j);
            sb.append(Const.tab);
            sb.append(a1 == -1 ? Const.MISSING_DATA_CHAR : a1);
            sb.append(isPhased(j) ? Const.phasedSep : Const.unphasedSep);
            sb.append(a2 == -1 ? Const.MISSING_DATA_CHAR : a2);
        }
        return sb.toString();
    }
}
