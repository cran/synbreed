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
import main.NuclearFamilies;

/**
 * <p>Constructs a new {@code MedMemGTGL} instance representing
 * genotype emission probabilities for a set of samples at a single marker.
 * The genotype emission probabilities are determined by either the
 * genotypes likelihoods or the called genotypes for the samples.
 * </p>
 * <p>Instances of class {@code MedMemGT} are immutable.
 * </p>
 */
public final class MedMemGTGL implements VcfEmission {

    private final Samples samples;
    private final Marker marker;
    private final VcfEmission gtData;
    private final VcfEmission glData;
    private final boolean preferGL;

    /**
     * Constructs a new {@code MedMemGTGL} instance representing
     * the specified VCF record's GL, PL, or GT format field data.
     * If both the GL and PL format codes are present, the PL format
     * field data will be ignored.
     *
     * @param rec a VCF record.
     * @param fam parent-offspring relationships.
     * @param usePhase {@code true} if phase information in the specified
     * VCF record will be used when a sample's genotype emission probabilities
     * are determined by a called genotype, and {@code false} if
     * phase information in the specified VCF file record will be ignored.
     * @param maxLR maximum likelihood ratio.  If the likelihood ratio between
     * two possible genotypes is larger than {@code maxLR}, then the
     * smaller likelihood is set to 0.0, unless the change would create a
     * Mendelian inconsistency in a parent-offspring trio or duo.  In such
     * a case the unmodified likelihoods are used for all members of the
     * inconsistent duo or trio.
     * @param preferGL {@code true} if genotype emission probabilities
     * should be determined from the VCF record's GL or PL format field if
     * the GL or PL format code is present, and {@code false} if
     * genotype emission probabilities should be determined from the VCF
     * record's a GT format field if the GT format code is present.
     *
     * @throws IllegalArgumentException if {@code rec.nSamples()==0}
     * @throws IllegalArgumentException if the VCF record does not have a
     * GL, PL, or GT format field
     * @throws IllegalArgumentException if
     * {@code record.samples().equals(fam.samples())==false}
     * @throws IllegalArgumentException if
     * {@code Float.isNaN(maxLR) || maxLR<=1.0f}
     * @throws NullPointerException if {@code vcf==null || fam==null}
     */
    public MedMemGTGL(VcfRecord rec, NuclearFamilies fam, boolean usePhase,
            float maxLR, boolean preferGL) {
        boolean hasGT = rec.hasFormat(BitSetGT.GT_FORMAT);
        boolean hasGL = rec.hasFormat(MedMemGL.GL_FORMAT);
        boolean hasPL = rec.hasFormat(MedMemGL.PL_FORMAT);
        if (hasGT==false && hasGL==false && hasPL==false) {
            String s = "missing GT FORMAT and GL/PL FORMAT: " + rec;
            throw new IllegalArgumentException(s);
        }
        this.samples = rec.samples();
        this.marker = rec.marker();
        this.gtData = hasGT ? new BitSetGT(rec, fam, usePhase) : null;
        this.glData = hasGL || hasPL ? new MedMemGL(rec, fam, maxLR) : null;
        this.preferGL = preferGL;
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
        return gtData!=null && (glData==null || preferGL==false)
                ? gtData.isRefData() : glData.isRefData();
    }

    @Override
    public boolean isMissingData() {
        return gtData!=null && (glData==null || preferGL==false)
                ? gtData.isMissingData() : glData.isMissingData();
    }

    @Override
    public boolean isPhased(int sample) {
        return gtData!=null && (glData==null || preferGL==false)
                ? gtData.isPhased(sample) : false;
    }

    @Override
    public byte allele1(int sample) {
        return gtData!=null && (glData==null || preferGL==false)
                ? gtData.allele1(sample) : -1;
    }

    @Override
    public byte allele2(int sample) {
        return gtData!=null && (glData==null || preferGL==false)
                ? gtData.allele2(sample) : -1;
    }

    @Override
    public float gl(int sample, byte a1, byte a2) {
        if (gtData==null) {
            return glData.gl(sample, a1, a2);
        }
        else if (glData==null) {
            return gtData.gl(sample, a1, a2);
        }
        else if (preferGL) {
            return glData.gl(sample, a1, a2);
        }
        else {
            return gtData.gl(sample, a1, a2);
        }
    }

    /**
     * Returns the data represented by {@code this} as a VCF file
     * record with a GT format field.
     * @return the data represented by {@code this} as a VCF file
     * record with a GT format field.
     */
    @Override
    public String toString() {
        if (gtData==null) {
            return glData.toString();
        }
        else if (glData==null) {
            return gtData.toString();
        }
        else {
            StringBuilder sb = new StringBuilder();
            sb.append(marker);
            sb.append(Const.tab);
            sb.append(Const.MISSING_DATA_CHAR);
            sb.append(Const.tab);
            sb.append("PASS");
            sb.append(Const.tab);
            sb.append(Const.MISSING_DATA_CHAR);
            sb.append(Const.tab);
            sb.append(BitSetGT.GT_FORMAT);
            sb.append(Const.colon);
            sb.append(MedMemGL.GL_FORMAT);
            for (int j=0, n=samples.nSamples(); j<n; ++j) {
                sb.append(Const.tab);
                byte a1 = allele1(j);
                byte a2 = allele2(j);
                if (a1>=0 && a2>=0) {
                    sb.append(a1);
                    sb.append(gtData.isPhased(j) ? Const.phasedSep : Const.unphasedSep);
                    sb.append(a2);
                }
                else {
                    sb.append(Const.MISSING_DATA_CHAR);
                    sb.append(Const.unphasedSep);
                    sb.append(Const.MISSING_DATA_CHAR);
                    int m = marker.nAlleles();
                    for (a2=0; a2<m; ++a2) {
                        for (a1=0; a1<=a2; ++a1) {
                            sb.append(a1==0 && a2==0 ? Const.colon : Const.comma);
                            double d = Math.log10(gl(j, a1, a2));
                            sb.append(MedMemGL.df.format(d));
                        }
                    }
                }
            }
            return sb.toString();
        }
    }
}
