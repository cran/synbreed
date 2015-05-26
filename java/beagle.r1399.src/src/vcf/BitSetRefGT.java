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
import java.util.BitSet;

/**
 * <p>Class {@code BitSetRefGT} represent represents genotype emission
 * probabilities for a set of reference samples at a single marker.
 * The genotype emission probabilities are determined by the called
 * genotypes for the samples.  All genotypes are required to be phased and
 * to have no missing alleles.
 * </p>
 * <p>Instances of class {@code BitSetRefGT} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class BitSetRefGT implements VcfEmission {

    /**
     * The VCF FORMAT code for genotype data: "GT".
     */
    private static final String GT_FORMAT = "GT";

    private final int bitsPerAllele;
    private final Marker marker;
    private final Samples samples;
    private final BitSet allele1;
    private final BitSet allele2;

    /**
     * Constructs a new {@code LowMemRefGT} instance representing
     * the specified VCF record's GT format field data.
     *
     * @param rec a VCF record.
     *
     * @throws IllegalArgumentException if {@code rec.nSamples()==0}.
     * @throws IllegalArgumentException if the VCF record does not have a
     * GT format field.
     * @throws IllegalArgumentException if any genotype has a missing allele
     * or if any genotype is unphased.
     * @throws NullPointerException if {@code rec==null}.
     */
    public BitSetRefGT(VcfRecord rec) {
        if (rec.nSamples()==0) {
            throw new IllegalArgumentException("missing sample data: " + rec);
        }
        if (rec.hasFormat(GT_FORMAT)==false) {
            throw new IllegalArgumentException("missing GT FORMAT: " + rec);
        }
        checkAlleles(rec);
        this.bitsPerAllele = bitsPerAllele(rec.marker());
        this.marker = rec.marker();
        this.samples = rec.samples();
        this.allele1 = new BitSet(rec.nSamples()*bitsPerAllele);
        this.allele2 = new BitSet(rec.nSamples()*bitsPerAllele);
        storeAlleles(rec, bitsPerAllele, allele1, allele2);
    }

    private static void checkAlleles(VcfRecord rec) {
        int nAlleles = rec.marker().nAlleles();
        for (int sample=0, n=rec.nSamples(); sample<n; ++sample) {
            byte a1 = rec.gt(sample, 0);
            byte a2 = rec.gt(sample, 1);
            if (rec.isPhased(sample)==false) {
                String sampleId = rec.samples().id(sample);
                String s = "Reference genotype is not phased for sample: "
                        + sampleId + " marker: " + rec.marker();
                throw new IllegalArgumentException(s);
            }
            if (a1 < 0 || a2 < 0) {
                String sampleId = rec.samples().id(sample);
                String s = "Reference genotype has a missing allele for sample: "
                        + sampleId + " marker: " + rec.marker();
                throw new IllegalArgumentException(s);
            }
            if (a1 >= nAlleles || a2 >= nAlleles) {
                String sampleId = rec.samples().id(sample);
                String s = "Invalid allele index (" + Math.max(a1, a2)
                        + ") for sample: " + sampleId + " marker: "
                        + rec.marker();
                throw new IllegalArgumentException(s);
            }
        }
    }

    private static int bitsPerAllele(Marker marker) {
        int nAllelesM1 = marker.nAlleles() - 1;
        int nStorageBits = Integer.SIZE - Integer.numberOfLeadingZeros(nAllelesM1);
        return nStorageBits;
    }

    private static void storeAlleles(VcfRecord rec, int bitsPerAllele,
            BitSet allele1, BitSet allele2) {
        int index1 = 0;
        int index2 = 0;
        for (int sample=0, n=rec.nSamples(); sample<n; ++sample) {
            byte a1 = rec.gt(sample, 0);
            byte a2 = rec.gt(sample, 1);

            int mask = 1;
            for (int k=0; k<bitsPerAllele; ++k) {
                if ((a1 & mask)==mask) {
                    allele1.set(index1);
                }
                ++index1;
                mask <<= 1;
            }

            mask = 1;
            for (int k=0; k<bitsPerAllele; ++k) {
                if ((a2 & mask)==mask) {
                    allele2.set(index2);
                }
                ++index2;
                mask <<= 1;
            }
        }
    }

    @Override
    public Samples samples() {
        return samples;
    }

    @Override
    public int nSamples() {
        return samples.nSamples();
    }

    @Override
    public Marker marker() {
        return marker;
    }


    @Override
    public boolean isMissingData() {
        return false;
    }

    @Override
    public boolean isRefData() {
        return true;
    }

    @Override
    public float gl(int sample, byte allele1, byte allele2) {
        boolean matches = (allele1==allele1(sample) && allele2==allele2(sample));
        return matches ? 1.0f : 0.0f;    }

    @Override
    public boolean isPhased(int sample) {
        return true;
    }

    @Override
    public byte allele1(int sample) {
        return allele(allele1, sample);
    }

    @Override
    public byte allele2(int sample) {
        return allele(allele2, sample);
    }

    private byte allele(BitSet bits, int sample) {
        int start = bitsPerAllele*sample;
        int end = start + bitsPerAllele;
        byte allele = 0;
        byte mask = 1;
        for (int j=start; j<end; ++j) {
            if (bits.get(j)) {
                allele += mask;
            }
            mask <<= 1;
        }
        return allele;
    }

    /**
     * Returns the data represented by {@code this} as a VCF
     * record with a GT format field.
     * @return the data represented by {@code this} as a VCF
     * record with a GT format field.
     */
    @Override
    public String toString() {
        String missingField = ".";
        StringBuilder sb = new StringBuilder();
        sb.append(marker);
        sb.append(Const.tab);
        sb.append(missingField);
        sb.append(Const.tab);
        sb.append("PASS");
        sb.append(Const.tab);
        sb.append(missingField);
        sb.append(Const.tab);
        sb.append("GT");
        for (int j=0, n=samples.nSamples(); j<n; ++j) {
            sb.append(Const.tab);
            sb.append(allele1(j));
            sb.append(Const.phasedSep);
            sb.append(allele2(j));
        }
        return sb.toString();
    }
}
