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
 * <p>Class {@code ByteArrayRefGT} represents phased, non-missing
 * genotypes for a list of reference samples at a single marker.
 * Genotype emission probabilities are determined by the sample
 * genotypes.
 * </p>
 * <p>Instances of class {@code ByteArrayRefGT} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class ByteArrayRefGT implements VcfEmission {

    private final Marker marker;
    private final Samples samples;
    private final byte[] allele1;
    private final byte[] allele2;

    /**
     * Creates a new {@code ByteArrayRefGT} instance from the
     * specified {@code VcfHeader} and string VCF record.
     * The VCF record must contain a GT FORMAT field, and
     * all genotypes must have phased, non-missing genotypes.
     *
     * @param vcfHeader meta-information lines and header line for the
     * specified VCF record
     * @param vcfRecord a VCF record corresponding to the specified
     * {@code vcfHeader} object
     *
     * @throws IllegalArgumentException if a format error is detected
     * in the VCF record
     * @throws IllegalArgumentException if any VCF record genotype
     * is unphased or missing
     * @throws NullPointerException if
     * {@code vcfHeader == null || vcfRecord == null}
     */
    public ByteArrayRefGT(VcfHeader vcfHeader, String vcfRecord) {
        this(new VcfRecGTParser(vcfHeader, vcfRecord));
    }

    /**
     * Creates a new {@code VcfEmission} instance from a VCF record
     * containing phased, non-missing genotypes for a list of reference samples.
     * The VCF record must contain a GT FORMAT field, and
     * all genotypes must have phased, non-missing genotypes.
     *
     * @param gtp a parser for the VCF record
     * @throws IllegalArgumentException if a format error is detected
     * in the VCF record
     * @throws IllegalArgumentException if any VCF record genotype
     * is unphased or missing
     * @throws NullPointerException if {@code gtp == null}
     */
    public ByteArrayRefGT(VcfRecGTParser gtp) {
        if (gtp.marker().nAlleles() > Byte.MAX_VALUE) {
            throw new IllegalArgumentException(gtp.marker().toString());
        }
        if (gtp.currentSample() > 0) {
            throw new IllegalArgumentException(
                    String.valueOf(gtp.currentSample()));
        }
        this.marker = gtp.marker();
        this.samples = gtp.samples();
        this.allele1 = new byte[gtp.nSamples()];
        this.allele2 = new byte[gtp.nSamples()];
        storeAlleles(gtp, allele1, allele2);
    }

    private static void storeAlleles(VcfRecGTParser gtp, byte[] allele1,
            byte[] allele2) {
        int nSamples = gtp.nSamples();
        for (int sample=0; sample<nSamples; ++sample) {
            byte a1 = (byte) gtp.allele1();
            byte a2 = (byte) gtp.allele2();
            if (gtp.isPhased()==false || a1 == -1 || a2 == -2) {
                String s = "Unphased or missing reference genotype at marker: "
                        + gtp.marker();
                throw new IllegalArgumentException(s);
            }
            allele1[sample] = a1;
            allele2[sample] = a2;
            if (sample + 1 < nSamples) {
                gtp.nextSample();
            }
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
        return allele1[sample];
    }

    @Override
    public int allele2(int sample) {
        return allele2[sample];
    }

    @Override
    public int allele(int hap) {
        int sample = hap/2;
        return (hap & 1) == 0 ? allele1(sample) : allele2(sample);
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
