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
import blbutil.StringUtil;
import java.text.DecimalFormat;
import java.util.Arrays;
import main.NuclearFamilies;

/**
 * <p>Class {@code LowMemGT} represents genotype emission
 * probabilities for a set of samples at a single marker.
 * The genotype emission probabilities are determined by the
 * genotype likelihoods for the samples.
 * </p>
 * <p>Instances of class {@code LowMemGT} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class MedMemGL implements VcfEmission {

    /**
     * The VCF FORMAT code for log-scaled genotype likelihood data: "GL".
     */
    public static final String GL_FORMAT = "GL";

    /**
     * The VCF FORMAT code for phred-scaled genotype likelihood data: "GL".
     */
    public static final String PL_FORMAT = "PL";

    static final DecimalFormat df = new DecimalFormat("#.####");

    private static final float MIN_LIKE = 1024*Float.MIN_VALUE;

    private final Marker marker;
    private final Samples samples;
    private final float[] like;

    /**
     * Constructs a new {@code MedMemGL} instance representing
     * the specified VCF record's GL or PL format field data.
     * If both GL and PL format codes are present, the PL format field data
     * will be ignored.
     *
     * @param rec a VCF record.
     * @param fam parent-offspring relationships.
     * @param maxLR maximum likelihood ratio.  If the likelihood ratio between
     * two possible genotypes is larger than {@code maxLR}, then the
     * smaller likelihood is set to 0.0, unless the change would create a
     * Mendelian inconsistency in a parent-offspring trio or duo.  In such
     * a case the unmodified likelihoods are used for all members of the
     * inconsistent duo or trio.
     *
     * @throws IllegalArgumentException if {@code rec.nSamples()==0}.
     * @throws IllegalArgumentException if the VCF record does not have a
     * "GL" or a "PL" format field.
     * @throws IllegalArgumentException if
     * {@code record.samples().equals(fam.samples())==false}
     * @throws IllegalArgumentException if
     * {@code Float.isNaN(maxLR) || maxLR<=1.0f}.
     * @throws NullPointerException if {@code vcf==null || fam==null}
     */
    public MedMemGL(VcfRecord rec, NuclearFamilies fam, float maxLR) {
        if (rec.nSamples()==0) {
            String s = "missing sample data: " + rec;
            throw new IllegalArgumentException(s);
        }
        if (rec.samples().equals(fam.samples())==false) {
            throw new IllegalArgumentException("inconsistent samples");
        }
        if (Float.isNaN(maxLR) || maxLR <= 1.0f) {
            throw new IllegalArgumentException("maxLR: " + maxLR);
        }
        if (rec.hasFormat(GL_FORMAT)==false
                && rec.hasFormat(PL_FORMAT)==false) {
            String s = "missing GL/PL format: " + rec;
            throw new IllegalArgumentException(s);
        }
        float minLR = 1.0f/maxLR;
        this.samples = rec.samples();
        this.marker = rec.marker();
        this.like = likelihoods(rec, fam, minLR);
    }

    private static float[] likelihoods(VcfRecord rec, NuclearFamilies fam,
            float minLR) {
        Marker marker = rec.marker();
        int nGt = marker.nGenotypes();
        boolean[] ba1 = new boolean[rec.marker().nAlleles()];
        boolean[] ba2 = new boolean[rec.marker().nAlleles()];
        float[] origLike = likelihoodsFromGL(rec);
        float[] adjLike = copyAndApplyMinLR(origLike, minLR);
        for (int j=0, n=fam.nDuos(); j<n; ++j) {
            int p = fam.duoParent(j);
            int o = fam.duoOffspring(j);
            if (duoIsConsistent(marker, origLike, minLR, p, o, ba1, ba2) == false) {
                int pBase = p*nGt;
                int oBase = o*nGt;
                for (int gt=0; gt<nGt; ++gt) {
                    adjLike[pBase + gt] = Math.max(MIN_LIKE, origLike[pBase + gt]);
                    adjLike[oBase + gt] = Math.max(MIN_LIKE, origLike[oBase + gt]);
                }
            }
        }
        for (int j=0, n=fam.nTrios(); j<n; ++j) {
            int f = fam.trioFather(j);
            int m = fam.trioMother(j);
            int o = fam.trioOffspring(j);
            if (trioIsConsistent(marker, origLike, minLR, f, m, o, ba1, ba2) == false) {
                int fBase = f*nGt;
                int mBase = m*nGt;
                int oBase = o*nGt;
                for (int gt=0; gt<nGt; ++gt) {
                    adjLike[fBase + gt] = Math.max(MIN_LIKE, origLike[fBase + gt]);
                    adjLike[mBase + gt] = Math.max(MIN_LIKE, origLike[mBase + gt]);
                    adjLike[oBase + gt] = Math.max(MIN_LIKE, origLike[oBase + gt]);
                }
            }
        }
        return adjLike;
    }

    private static float[] likelihoodsFromGL(VcfRecord rec) {
        int nGt = rec.marker().nGenotypes();
        String[] dataGL = rec.hasFormat(GL_FORMAT) ? rec.formatData(GL_FORMAT) : null;
        String[] dataPL = rec.hasFormat(PL_FORMAT) ? rec.formatData(PL_FORMAT) : null;
        double[] doubleLike = new double[nGt];
        float[] floatLike = new float[rec.nSamples()*nGt];
        int floatLikeIndex = 0;
        for (int j=0, n=rec.nSamples(); j<n; ++j) {
            Arrays.fill(doubleLike, 0.0);
            if (dataGL != null) {
                String[] fields = getGL(GL_FORMAT, dataGL, j, nGt, rec);
                for (int k=0; k<nGt; ++k) {
                    doubleLike[k] = GL2Like(fields[k]);
                }
            }
            else if (dataPL != null) {
                String[] fields = getGL(PL_FORMAT, dataPL, j, nGt, rec);
                for (int k=0; k<nGt; ++k) {
                    doubleLike[k] = PL2Like(fields[k]);
                }
            }
            rescaleToMax1(doubleLike);
            for (int gt=0; gt<nGt; ++gt) {
                floatLike[floatLikeIndex++] = (float) doubleLike[gt];
            }
        }
        assert floatLikeIndex==floatLike.length;
        return floatLike;
    }

    private static boolean duoIsConsistent(Marker marker, float[] like,
            float minLR, int parent, int offspring,
            boolean[] parentAlleles, boolean[] offspringAlleles) {
        alleles(marker, like, minLR, parent, parentAlleles);
        alleles(marker, like, minLR, offspring, offspringAlleles);
        for (int j=0, n=marker.nAlleles(); j<n; ++j) {
            if (parentAlleles[j] && offspringAlleles[j]) {
                return true;
            }
        }
        return false;
    }

    private static boolean trioIsConsistent(Marker marker, float[] like,
            float minLR, int father, int mother, int offspring, boolean[] fatherAlleles,
            boolean[] motherAlleles) {
        int nAlleles = marker.nAlleles();
        alleles(marker, like, minLR, father, fatherAlleles);
        alleles(marker, like, minLR, mother, motherAlleles);
        int base = offspring*marker.nGenotypes();
        for (byte a1=0; a1<nAlleles; ++a1) {
            if (fatherAlleles[a1]) {
                for (byte a2=0; a2<nAlleles; ++a2) {
                    if (motherAlleles[a2]) {
                        int gt = BasicGL.genotype(a1, a2);
                        if (like[base + gt] >= minLR) {
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    }

    private static void alleles(Marker marker, float[] like, float minLR,
            int sample, boolean[] alleles) {
        Arrays.fill(alleles, false);
        int base = sample*marker.nGenotypes();
        for (byte a1=0; a1<alleles.length; ++a1) {
            for (byte a2=a1; a2<alleles.length; ++a2) {
                int gt = BasicGL.genotype(a1, a2);
                if ((like[base + gt] >= minLR)) {
                    alleles[a1] |= true;
                    alleles[a2] |= true;
                }
            }
        }
    }

    private static void rescaleToMax1(double[] like) {
        double max = max(like);
        if (max == 0.0f) {
            Arrays.fill(like, 1.0);
        }
        else {
            for (int j=0; j<like.length; ++j) {
                like[j] /= max;
            }
        }
    }

    /* returns max{double[] like, double 0.0} */
    private static double max(double[] like) {
        double max = 0.0;
        for (int k=0; k<like.length; ++k) {
            if (like[k] > max) {
                max = like[k];
            }
        }
        return max;
    }

    /* set likelihoods less than minLR to 0.0 */
    private static float[] copyAndApplyMinLR(float[] like, float minLR) {
        float[] fa = like.clone();
        for (int j=0; j<fa.length; ++j) {
            if (fa[j] < minLR) {
                fa[j] = 0.0f;
            }
        }
        return fa;
    }

    private static String[] getGL(String format, String[] sampleData,
            int sample, int nGt, VcfRecord record) {
        if (sampleData[sample].equals(Const.MISSING_DATA_STRING)) {
            String[] fields = new String[nGt];
            Arrays.fill(fields, "0");
            return fields;
        }
        else {
            String[] subfields = StringUtil.getFields(sampleData[sample],
                    Const.comma);
            if (subfields.length!=nGt) {
                String s = "unexpected number of " + format + " subfields: "
                        + record.sampleFormatData(format, sample) + Const.nl
                        + record;
                throw new IllegalArgumentException(s);
            }
            for (String subfield : subfields) {
                if (subfield.equals(Const.MISSING_DATA_STRING)) {
                    String s = "missing subfield in " + format + " field: "
                        + record.sampleFormatData(format, sample) + Const.nl
                        + record;
                    throw new IllegalArgumentException(s);
                }
            }
            return subfields;
        }
    }

    private static double GL2Like(String gl) {
        return Math.pow(10.0, Double.parseDouble(gl));
    }

    private static double PL2Like(String pl) {
        return VcfRecord.fromPhred(Integer.parseInt(pl));
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
        return false;
    }

    @Override
    public boolean isMissingData() {
        return false;
    }

    @Override
    public boolean isPhased(int sample) {
        return false;
    }

    @Override
    public byte allele1(int sample) {
        return -1;
    }

    @Override
    public byte allele2(int sample) {
        return -1;
    }

    @Override
    public float gl(int sample, byte allele1, byte allele2) {
        int n = marker.nAlleles();
        if (allele1 < 0 || allele2 < 0 || allele1 >= n || allele2 >= n) {
            String s = allele1 + " " + allele2 + " " + n;
            throw new ArrayIndexOutOfBoundsException(s);
        }
        int gt = VcfRecord.gtIndex(allele1, allele2);
        return like[(sample*marker.nGenotypes()) + gt];
    }

    /**
     * Returns the data represented by {@code this} as a VCF file
     * record with a one format field (GL).
     * @return the data represented by {@code this} as a VCF file
     * record with a one format field (GL).
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
        sb.append("GL");
        for (int sample=0, n=samples.nSamples(); sample<n; ++sample) {
            for (int gt=0, m=marker.nGenotypes(); gt<m; ++gt) {
                sb.append(gt==0 ? Const.tab : Const.comma);
                double d = Math.log10(like[(sample*marker.nGenotypes()) + gt]);
                sb.append(df.format(d));
            }
        }
        return sb.toString();
    }
}
