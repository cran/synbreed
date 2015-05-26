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

import beagleutil.SampleIds;
import beagleutil.Samples;
import blbutil.Const;
import java.util.BitSet;
import main.NuclearFamilies;

/**
 * <p>Class {@code BitSetGT} represents genotype emission
 * probabilities for a set of samples at a single marker.
 * The genotype emission probabilities are determined by the called
 * genotypes for the samples.
 * </p>
 * <p>Class {@code BitSetGT} is designed to use a relatively small amount
 * of memory.
 * </p>
 * <p>Instances of class {@code BitSetGT} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class BitSetGT implements VcfEmission {

    /**
     * The VCF FORMAT code for genotype data: "GT".
     */
    public static final String GT_FORMAT = "GT";

    private final byte bitsPerAllele;
    private final Samples samples;
    private final Marker marker;
    private final boolean isRefData;

    private final BitSet allele1;
    private final BitSet allele2;
    private final BitSet isMissing1;
    private final BitSet isMissing2;
    private final BitSet isPhased;

    /**
     * Constructs a new {@code LowMemGT} instance representing
     * the specified VCF record's GT format field data.
     *
     * @param rec a VCF record.
     * @param usePhase {@code true} if phase information in the specified
     * VCF file record will be used, and {@code false} if phase
     * information will be ignored.
     *
     * @throws IllegalArgumentException if {@code rec.nSamples()==0}.
     * @throws IllegalArgumentException if the VCF record does not have a
     * "GT" format field.
     * @throws NullPointerException if {@code rec==null}.
     */
    public BitSetGT(VcfRecord rec, boolean usePhase) {
        this(rec);
        setBits(rec, usePhase, bitsPerAllele, allele1, allele2, isMissing1,
                isMissing2, isPhased);
    }

    /**
     * Constructs a new {@code LowMemGT} instance representing
     * the specified VCF record's GT format field data.
     *
     * @param rec a VCF file record.
     * @param fam parent-offspring relationships.
     * @param usePhase {@code true} if phase information in the specified
     * VCF file record will be used, and {@code false} if phase
     * information in the specified VCF file record will be ignored.
     *
     * @throws IllegalArgumentException if
     * {@code rec.nSamples()==0|| rec.samples().equals(fam.samples())==false}.
     * @throws IllegalArgumentException if the VCF record does not have a
     * GT format field.
     * @throws NullPointerException if {@code rec==null || fam==null}.
     */
    public BitSetGT(VcfRecord rec, NuclearFamilies fam, boolean usePhase) {
        this(rec);
        if (rec.samples().equals(fam.samples())==false) {
            throw new IllegalArgumentException("inconsistent samples");
        }
        setBits(rec, usePhase, bitsPerAllele, allele1, allele2, isMissing1,
                isMissing2, isPhased);
        removeMendelianInconsistencies(rec, fam, isPhased, isMissing1,
                isMissing2);
    }

    private BitSetGT(VcfRecord rec) {
        int nSamples = rec.nSamples();
        if (nSamples==0) {
            String s = "missing sample data: " + rec;
            throw new IllegalArgumentException(s);
        }
        if (rec.hasFormat(GT_FORMAT)==false) {
            String s = "missing GT FORMAT: " + rec;
            throw new IllegalArgumentException(s);
        }
        this.bitsPerAllele = bitsPerAllele(rec.marker());
        this.samples = rec.samples();
        this.marker = rec.marker();
        this.isRefData = isRef(rec);

        this.allele1 = new BitSet(nSamples * bitsPerAllele);
        this.allele2 = new BitSet(nSamples * bitsPerAllele);
        this.isMissing1 = new BitSet(nSamples);
        this.isMissing2 = new BitSet(nSamples);
        this.isPhased = new BitSet(nSamples);
    }

    private static boolean isRef(VcfRecord rec) {
        for (int j=0, n=rec.nSamples(); j<n; ++j) {
            if (rec.isPhased(j)==false || rec.gt(j, 0)<0 || rec.gt(j,1)<0) {
                return false;
            }
        }
        return true;
    }

    private static void setBits(VcfRecord rec, boolean usePhase,
            int bitsPerAllele, BitSet allele1, BitSet allele2,
            BitSet isMissing1, BitSet isMissing2, BitSet isPhased) {
        int index1 = 0;
        int index2 = 0;
        for (int j=0, n=rec.nSamples(); j<n; ++j) {
            if (usePhase && rec.isPhased(j)) {
                isPhased.set(j);
            }
            byte a1 = rec.gt(j, 0);
            byte a2 = rec.gt(j, 1);
            if (a1 < 0) {
                isMissing1.set(j);
                index1 += bitsPerAllele;
            }
            else {
                int mask = 1;
                for (int k=0; k<bitsPerAllele; ++k) {
                    if ((a1 & mask)==mask) {
                        allele1.set(index1);
                    }
                    ++index1;
                    mask <<= 1;
                }
            }

            if (a2 < 0) {
                isMissing2.set(j);
                index2 += bitsPerAllele;
            }
            else {
                int mask = 1;
                for (int k=0; k<bitsPerAllele; ++k) {
                    if ((a2 & mask)==mask) {
                        allele2.set(index2);
                    }
                    ++index2;
                    mask <<= 1;
                }
            }
        }
    }

    private static byte bitsPerAllele(Marker marker) {
        int nAllelesM1 = marker.nAlleles() - 1;
        int nStorageBits = Integer.SIZE - Integer.numberOfLeadingZeros(nAllelesM1);
        return (byte) nStorageBits;
    }

    /*
     * Sets phase to unknown for all parent-offspring relationships, and sets
     * all genotypes in a duo or trio genotypes to missing if a Mendelian
     * inconsistency is found.
     */
    private static void removeMendelianInconsistencies(VcfRecord rec,
            NuclearFamilies fam, BitSet isPhased, BitSet isMissing1,
            BitSet isMissing2) {
        for (int j=0, n=fam.nDuos(); j<n; ++j) {
            int p = fam.duoParent(j);
            int o = fam.duoOffspring(j);
            isPhased.clear(p);
            isPhased.clear(o);
            if (duoIsConsistent(rec, p, o) == false) {
                logDuoInconsistency(rec, p, o);
                isMissing1.set(p);
                isMissing2.set(p);
                isMissing1.set(o);
                isMissing2.set(o);
            }
        }
        for (int j=0, n=fam.nTrios(); j<n; ++j) {
            int f = fam.trioFather(j);
            int m = fam.trioMother(j);
            int o = fam.trioOffspring(j);
            isPhased.clear(f);
            isPhased.clear(m);
            isPhased.clear(o);
            if (trioIsConsistent(rec, f, m, o) == false) {
                logTrioInconsistency(rec, f, m, o);
                isMissing1.set(f);
                isMissing2.set(f);
                isMissing1.set(m);
                isMissing2.set(m);
                isMissing1.set(o);
                isMissing2.set(o);
            }
        }
    }

    private static boolean duoIsConsistent(VcfRecord rec, int parent,
            int offspring) {
        byte p1 = rec.gt(parent, 0);
        byte p2 = rec.gt(parent, 1);
        byte o1 = rec.gt(offspring, 0);
        byte o2 = rec.gt(offspring, 1);
        boolean alleleMissing = (p1<0 || p2<0 || o1<0 || o2<0);
        return (alleleMissing || p1==o1 || p1==o2 || p2==o1 || p2==o2);
    }

    private static boolean trioIsConsistent(VcfRecord rec, int father,
            int mother, int offspring) {
        byte f1 = rec.gt(father, 0);
        byte f2 = rec.gt(father, 1);
        byte m1 = rec.gt(mother, 0);
        byte m2 = rec.gt(mother, 1);
        byte o1 = rec.gt(offspring, 0);
        byte o2 = rec.gt(offspring, 1);
        boolean fo1 = (o1<0 || f1<0 || f2<0 || o1==f1 || o1==f2);
        boolean mo2 = (o2<0 || m1<0 || m2<0 || o2==m1 || o2==m2);
        if (fo1 && mo2) {
            return true;
        }
        else {
            boolean fo2 = (o2<0 || f1<0 || f2<0 || o2==f1 || o2==f2);
            boolean mo1 = (o1<0 || m1<0 || m2<0 || o1==m1 || o1==m2);
            return (fo2 && mo1);
        }
    }

    private static void logDuoInconsistency(VcfRecord rec, int parent,
            int offspring) {
        StringBuilder sb = new StringBuilder(80);
        sb.append("WARNING: Inconsistent duo genotype set to missing");
        sb.append(Const.tab);
        sb.append(rec.marker());
        sb.append(Const.colon);
        sb.append(rec.samples().id(parent));
        sb.append(Const.tab);
        sb.append(rec.samples().id(offspring));
        main.Logger.getInstance().println(sb.toString());
    }

    private static void logTrioInconsistency(VcfRecord rec, int father,
            int mother, int offspring) {
        StringBuilder sb = new StringBuilder(80);
        sb.append("WARNING: Inconsistent trio genotype set to missing");
        sb.append(Const.tab);
        sb.append(rec.marker());
        sb.append(Const.tab);
        sb.append(rec.samples().id(father));
        sb.append(Const.tab);
        sb.append(rec.samples().id(mother));
        sb.append(Const.tab);
        sb.append(rec.samples().id(offspring));
        main.Logger.getInstance().println(sb.toString());
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
    public boolean isRefData() {
        return isRefData;
    }

    @Override
    public boolean isMissingData() {
        int n = (isMissing1.cardinality() + isMissing2.cardinality());
        return n==2*samples.nSamples();
    }

    @Override
    public boolean isPhased(int sample) {
        return isPhased.get(sample);
    }

    @Override
    public byte allele1(int sample) {
        return isMissing1.get(sample) ? -1 : allele(allele1, sample);
    }

    @Override
    public byte allele2(int sample) {
        return isMissing2.get(sample) ? -1 : allele(allele2, sample);
    }

    @Override
    public float gl(int sample, byte a1, byte a2) {
        if ( a1 < 0 || a1 >= marker.nAlleles())  {
            String s = "invalid alleles: (" + a1 + "): " + marker;
            throw new IllegalArgumentException(s);
        }
        if ( a2 < 0 || a2 >= marker.nAlleles()) {
            String s = "invalid alleles: (" + a2 + "): " + marker;
            throw new IllegalArgumentException(s);
        }
        if (isMissing1.get(sample) && isMissing2.get(sample)) {
            return 1.0f;
        }
        else if (isMissing1.get(sample) ^ isMissing2.get(sample)) {
            byte obsA1 = allele1(sample);
            byte obsA2 = allele2(sample);
            boolean consistent = (obsA1<0 || obsA1==a1) && (obsA2<0 || obsA2==a2);
            if (isPhased.get(sample)==false && consistent==false) {
                consistent = (obsA1<0 || obsA1==a2) && (obsA2<0 || obsA2==a1);
            }
            return consistent ? 1.0f : 0.0f;
        }
        else {
            byte obsA1 = allele(allele1, sample);
            byte obsA2 = allele(allele2, sample);
            if (isPhased.get(sample)) {
                return (obsA1==a1 && obsA2==a2) ? 1.0f : 0.0f;
            }
            else {
                boolean isConsistent = (obsA1==a1 && obsA2==a2)
                        || (obsA1==a2 && obsA2==a1);
                return isConsistent ? 1.0f : 0.0f;
            }
        }
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
        StringBuilder sb = new StringBuilder();
        sb.append(marker);
        sb.append(Const.tab);
        sb.append(Const.MISSING_DATA_CHAR);
        sb.append(Const.tab);
        sb.append("PASS");
        sb.append(Const.tab);
        sb.append(Const.MISSING_DATA_CHAR);
        sb.append(Const.tab);
        sb.append("GT");
        for (int j=0, n=samples.nSamples(); j<n; ++j) {
            sb.append(Const.tab);
            sb.append(isMissing1.get(j) ? Const.MISSING_DATA_CHAR : allele(allele1, j));
            sb.append(isPhased.get(j) ? Const.phasedSep : Const.unphasedSep);
            sb.append(isMissing2.get(j) ? Const.MISSING_DATA_CHAR : allele(allele2, j));
        }
        return sb.toString();
    }
}
