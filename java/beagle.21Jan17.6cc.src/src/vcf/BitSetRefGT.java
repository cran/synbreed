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
import java.util.BitSet;

/**
 * <p>Class {@code BitSetRefGT} represents genotype emission
 * probabilities for a list reference samples with phased, non-missing
 * genotypes at a single marker.
 * The genotype emission probabilities are determined by the called
 * genotypes for the reference samples.
 * </p>
 * <p>Class {@code BitSetRefGT} stores alleles using
 * {@code java.util.BitSet} objects.
 * </p>
 * <p>Instances of class {@code BitSetRefGT} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class BitSetRefGT implements VcfEmission {

    private final int bitsPerAllele;
    private final Marker marker;
    private final Samples samples;
    private final BitSet allele1;
    private final BitSet allele2;

    /**
     * Creates a new {@code BitSetRefGT} instance from a VCF record
     * storing phased, non-missing reference genotypes.
     *
     * @param vcfHeader meta-information lines and header line for the
     * specified VCF record
     * @param vcfRecord a VCF record corresponding to the specified
     * {@code vcfHeader} object
     *
     * @throws IllegalArgumentException if a format error is detected
     * in the VCF record or if any allele is missing or unphased
     * @throws IllegalArgumentException if {@code rec.nSamples() == 0}
     * @throws IllegalArgumentException if the header line
     * or VCF record does not have a "GT" format field
     * @throws NullPointerException if
     * {@code vcfHeader == null || vcfRecord == null}
     */
    public BitSetRefGT(VcfHeader vcfHeader, String vcfRecord) {
        this(new VcfRecGTParser(vcfHeader, vcfRecord));
    }

    /**
     * Creates a new {@code VcfEmission} instance from a VCF record
     * containing phased, non-missing genotypes for a list of reference samples.
     * @param gtp a parser for the VCF record
     * @throws IllegalArgumentException if a format error, a missing genotype,
     * or an unphased genotype is detected in the VCF record
     * @throws NullPointerException if {@code gtp==null}
     */
    public BitSetRefGT(VcfRecGTParser gtp) {
        this.marker = gtp.marker();
        this.samples = gtp.samples();
        this.bitsPerAllele = bitsPerAllele(marker);
        this.allele1 = new BitSet(gtp.nSamples()*bitsPerAllele);
        this.allele2 = new BitSet(gtp.nSamples()*bitsPerAllele);
        storeAlleles(gtp, bitsPerAllele, allele1, allele2);
    }

    private static int bitsPerAllele(Marker marker) {
        int nAllelesM1 = marker.nAlleles() - 1;
        int nStorageBits = Integer.SIZE - Integer.numberOfLeadingZeros(nAllelesM1);
        return nStorageBits;
    }

    private static void storeAlleles(VcfRecGTParser gtp, int bitsPerAllele,
            BitSet allele1, BitSet allele2) {
        int nSamples = gtp.nSamples();
        for (int sample=0; sample<nSamples; ++sample) {
            int a1 = gtp.allele1();
            int a2 = gtp.allele2();
            if (gtp.isPhased()==false || a1 == -1 || a2 == -2) {
                String s = "Unphased or missing reference genotype at marker: "
                        + gtp.marker();
                throw new IllegalArgumentException(s);
            }
            storeAllele(allele1, sample, bitsPerAllele, a1);
            storeAllele(allele2, sample, bitsPerAllele, a2);
            if (sample + 1 < nSamples) {
                gtp.nextSample();
            }
        }
    }

    private static void storeAllele(BitSet alleles, int sample,
            int bitsPerAllele, int allele) {
        int index = sample*bitsPerAllele;
        int mask = 1;
        for (int k=0; k<bitsPerAllele; ++k) {
            if ((allele & mask)==mask) {
                alleles.set(index);
            }
            ++index;
            mask <<= 1;
        }
    }

    @Override
    public boolean isPhased(int sample) {
        if (sample < 0  || sample >= this.nSamples()) {
            throw new IllegalArgumentException(String.valueOf(sample));
        }
        return true;
    }

    /**
     * Returns the samples. The returned samples are the filtered samples
     * after all sample exclusions.
     *
     * @return the samples.
     */
    @Override
    public Samples samples() {
        return samples;
    }

    @Override
    public int nSamples() {
        return samples.nSamples();
    }

    @Override
    public int nHaps() {
        return 2*samples.nSamples();
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
        if (allele1 < 0 || allele1 >= marker.nAlleles()) {
            throw new IndexOutOfBoundsException(String.valueOf(allele1));
        }
        if (allele2 < 0 || allele2 >= marker.nAlleles()) {
            throw new IndexOutOfBoundsException(String.valueOf(allele2));
        }
        boolean matches = (allele1==allele1(sample) && allele2==allele2(sample));
        return matches ? 1.0f : 0.0f;
    }

    @Override
    public int allele1(int sample) {
        return allele(allele1, bitsPerAllele, sample);
    }

    @Override
    public int allele2(int sample) {
        return allele(allele2, bitsPerAllele, sample);
    }

    @Override
    public int allele(int hap) {
        int sample = hap/2;
        return (hap & 1) == 0 ? allele1(sample) : allele2(sample);
    }

    private int allele(BitSet bits, int bitsPerAllele, int sample) {
        if (sample >= samples.nSamples()) {
            throw new IndexOutOfBoundsException(String.valueOf(sample));
        }
        int start = bitsPerAllele*sample;
        int end = start + bitsPerAllele;
        int allele = 0;
        int mask = 1;
        for (int j=start; j<end; ++j) {
            if (bits.get(j)) {
                allele += mask;
            }
            mask <<= 1;
        }
        return allele;
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
