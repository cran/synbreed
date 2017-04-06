/*
 * Copyright (C) 2015 browning
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package vcf;

import beagleutil.Samples;
import blbutil.IntArray;

/**
 * <p>Class {@code SeqCodedRefGT}  represents phased, non-missing
 * genotypes for a list of reference samples at a single marker.
 * Genotype emission probabilities are determined by the sample
 * genotypes.
 * </p>
 * <p>Instances of class {@code SeqCodedRefGT} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class SeqCodedRefGT implements VcfEmission {

    private final Marker marker;
    private final Samples samples;
    private final IntArray hapToSeq;
    private final IntArray seqToAllele;

    /**
     * Creates a new {@code SeqCodedRefGT} instance with phased,
     * non-missing genotypes from the specified marker, samples,
     * and haplotype alleles.
     * @param marker the marker
     * @param samples the samples
     * @param hapToSeq an array whose {@code j}-th element is the index
     * of the distinct allele sequence carried by the {@code j}-th haplotype
     * @param seqToAllele an array whose {@code j}-th element is the marker
     * allele carried by the {@code j}-th distinct allele sequence
     *
     * @throws IllegalArgumentException if
     * {@code hapToSeq.size() != 2*samples.nSamples()}
     * @throws IndexOutOfBoundsException if any element of {@code hapToSeq}
     * is negative or greater than or equal to {@code seqToAllele.size()}
     * @throws IndexOutOfBoundsException if any element of {@code seqToAllele}
     * is negative or greater than or equal to {@code marker.nAlleles()}
     * @throws NullPointerException if any parameter is {@code null}
     */
    public SeqCodedRefGT(Marker marker, Samples samples, IntArray hapToSeq,
        IntArray seqToAllele) {
        checkData(marker, samples, hapToSeq, seqToAllele);
        this.marker = marker;
        this.samples = samples;
        this.hapToSeq = hapToSeq;
        this.seqToAllele = seqToAllele;
    }

    private static void checkData(Marker marker, Samples samples,
            IntArray hapToSeq, IntArray seqToAllele) {
        String err = "inconsistent data";
        int nHaps = hapToSeq.size();
        if (hapToSeq.size() != 2*samples.nSamples()) {
            throw new IllegalArgumentException(err);
        }
        for (int j=0; j<nHaps; ++j) {
            marker.allele(seqToAllele.get(hapToSeq.get(j)));
        }
    }

    @Override
    public int nSamples() {
        return samples.nSamples();
    }

    @Override
    public Samples samples() {
        return samples;
    }

    @Override
    public int nHaps() {
        return hapToSeq.size();
    }

    @Override
    public int nHapPairs() {
        return samples.nSamples();
    }

    @Override
    public Marker marker() {
        return marker;
    }

    @Override
    public boolean isRefData() {
        return true;
    }

    @Override
    public float gl(int sample, int allele1, int allele2) {
        boolean match = allele1 == allele1(sample)
                && allele2 == allele2(sample);
        return match ? 1f : 0f;
    }

    @Override
    public boolean isPhased(int sample) {
        return true;
    }

    @Override
    public int allele1(int sample) {
        return seqToAllele.get(hapToSeq.get(2*sample));
    }

    @Override
    public int allele2(int sample) {
        return seqToAllele.get(hapToSeq.get(2*sample + 1));
    }

    @Override
    public int allele(int hap) {
        return seqToAllele.get(hapToSeq.get(hap));
    }

    @Override
    public int nAlleles() {
        return this.marker().nAlleles();
    }

    @Override
    public boolean storesNonMajorIndices() {
        return false;
    }

    @Override
    public int majorAllele() {
        String s = "this.storesNonMajorIndices()==false";
        throw new UnsupportedOperationException(s);
    }

    @Override
    public int alleleCount(int allele) {
        String s = "this.storesNonMajorIndices()==false";
        throw new UnsupportedOperationException(s);
    }

    @Override
    public int hapIndex(int allele, int copy) {
        String s = "this.storesNonMajorIndices()==false";
        throw new UnsupportedOperationException(s);
    }

    /**
     * Returns the data represented by {@code this} as a VCF
     * record with a GT format field. The returned VCF record
     * will have missing QUAL and INFO fields, will have "PASS"
     * in the filter field, and will have a GT format field.
     * @return the data represented by {@code this} as a VCF
     * record with a GT format field
     */
    @Override
    public String toString() {
        return toVcfRec();
    }
}
