/*
 * Copyright (C) 2016 Brian L. Browning
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

import blbutil.Const;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.math.MathContext;
import java.text.DecimalFormat;
import java.util.Arrays;

/**
 * <p>Class {@code VcfRecBuilder} contains methods for constructing
 * and printing a VCF record in VCF 4.2 format.  The FORMAT field data
 * for each sample is added sequentially to the record via the
 * {@code addSampleData()} method.
 *
 * </p>
 * <p>Instances of class {@code VcfRecBuilder} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class VcfRecBuilder {

    /**
     * The default initial size for the string buffer, which is 50
     * characters.
     */
    public static final int DEFAULT_INIT_SIZE = 50;

    private static final MathContext MC2 = new MathContext(2);
    private static final BigDecimal ONE = new BigDecimal(1.0);
    private static final String[] DS_VALS = dsProbs();
    private static final String[] R2_VALS = r2Probs(DS_VALS);

    private final StringBuilder sb;
    private final R2Estimator r2Est;
    private final double[] gt3Probs;

    private Marker marker;
    private boolean printDS;
    private boolean printGP;
    private int nAlleles;
    private int nGenotypes;
    private int allele1;
    private int allele2;
    private double[] gtProbs;
    private double[] dose;
    private double[] cumAlleleProbs;

    private static String[] dsProbs() {
        DecimalFormat df = new DecimalFormat("#.##");
        String[] probs = new String[201];
        for (int j=0;j<probs.length; ++j) {
            probs[j] = df.format(j/100.0);
        }
        return probs;
    }

    private static String[] r2Probs(String[] dsProbs) {
        DecimalFormat df = new DecimalFormat("0.00");
        String[] probs = Arrays.copyOf(dsProbs, 101);
        for (int j=0;j<probs.length; ++j) {
            if (probs[j].length()!=4) {
                probs[j] = df.format(j/100.0);
            }
        }
        return probs;
    }

    /**
     * Constructs a new {@code VcfRecBuilder} instance with initial buffer
     * size equal to {@code VcfRecBuilder.DEFAULT_INIT_SIZE}.
     */
    public VcfRecBuilder() {
        this(DEFAULT_INIT_SIZE);
    }

    /**
     * Constructs a new {@code VcfRecBuilder} instance with the specified
     * initial buffer size.
     *
     * @param initSize the initial buffer size
     * @throws NegativeArraySizeException if {@code initCapacity < 0}
     */
    public VcfRecBuilder(int initSize) {
        this.sb = new StringBuilder(initSize);
        this.r2Est = new R2Estimator();
        this.gt3Probs = new double[3];
    }

    /**
     * Clears existing data, and sets the current marker to the specified
     * marker.  If the FORMAT field contains a DS or GP subfield,
     * the INFO field will include the AR2 (allele r2), DR2 (dose r2), and
     * AF (ALT allele frequency) subfields.
     * @param marker the marker to which data will be added
     * @param printDS {@code true} if the FORMAT field in the VCF record for
     * this marker will include a DS subfield, and {@code false} otherwise
     * @param printGP {@code true} if the FORMAT field in the VCF record for
     * this marker will include a GP subfield, and {@code false} otherwise
     * @throws NullPointerException if {@code marker == null}
     */
    public void reset(Marker marker, boolean printDS, boolean printGP) {
        this.sb.setLength(0);
        this.r2Est.clear();
        this.marker = marker;
        this.printDS = printDS;
        this.printGP = printGP;
        this.allele1 = -1;
        this.allele2 = -1;
        this.nAlleles = marker.nAlleles();
        this.nGenotypes = marker.nGenotypes();
        this.gtProbs = new double[marker.nGenotypes()];
        this.cumAlleleProbs = new double[marker.nAlleles()];
        this.dose = new double[marker.nAlleles()];
    }

    /**
     * Returns the current marker.  Returns {@code null} if
     * {@code this.reset()} has not been previously invoked.
     * @return the current marker.
     */
    public Marker marker() {
        return marker;
    }

    /**
     * Returns {@code true} if the FORMAT field in the VCF record for
     * this marker includes a DS subfield, and {@code false} otherwise
     * @return {@code true} if the FORMAT field in the VCF record for
     * this marker includes a DS subfield
     */
    public boolean printDS() {
        return printDS;
    }

    /**
     * Returns {@code true} if the FORMAT field in the VCF record for
     * this marker includes a GP subfield, and {@code false} otherwise
     * @return {@code true} if the FORMAT field in the VCF record for
     * this marker includes a GP subfield
     */
    public boolean printGP() {
        return printGP;
    }

    /**
     * Adds the FORMAT field for a sample to the VCF record for the current
     * marker.  If the specified posterior genotype probabilities do not
     * sum to 1.0, the specified array will normalized to sum to 1.0.
     * @param gtypeProbs the posterior genotype probabilities
     * @throws IllegalArgumentException if
     * {@code gtProbs.length != this.marker().nGenotypes()}
     * @throws IllegalArgumentException if any element of the specified
     * array is not a finite non-negative number
     * @throws IllegalStateException if {@code this.marker() == null}
     * @throws NullPointerException if {@code gtProbs == null}
     */
    public void addSampleData(double[] gtypeProbs) {
        if (marker==null) {
            throw new IllegalStateException();
        }
        Arrays.fill(gt3Probs, 0.0);
        Arrays.fill(dose, 0.0);
        int maxGt = maxIndex(gtypeProbs, nGenotypes);
        int gt = 0;
        for (int a2=0; a2<nAlleles; ++a2) {
            for (int a1=0; a1<=a2; ++a1) {
                double gtProb = gtypeProbs[gt];
                if (gtProb < 0) {
                    throw new IllegalArgumentException(String.valueOf(gtProb));
                }
                gtProbs[gt] = gtProb;
                dose[a1] += gtProb;
                dose[a2] += gtProb;
                if (a2==0) {
                    gt3Probs[0] += gtProb;
                }
                else {
                    gt3Probs[(a1==0) ? 1 : 2] += gtProb;
                }
                if (gt==maxGt) {
                    allele1 = a1;
                    allele2 = a2;
                }
                ++gt;
            }
        }
        addToCumAlleleProbs(dose);
        r2Est.addSampleData(gt3Probs);
        appendFormatData(false);
    }

    private int maxIndex(double[] da, int expLength) {
        if (da.length != expLength) {
            throw new IllegalArgumentException(String.valueOf(da.length));
        }
        int maxIndex = 0;
        double sum = 0;
        for (int j=0; j<da.length; ++j) {
            if (da[j] < 0 || Double.isFinite(da[j])==false) {
                throw new IllegalArgumentException(String.valueOf(da[j]));
            }
            sum += da[j];
            if (da[j] > da[maxIndex]) {
                maxIndex = j;
            }
        }
        if (sum != 1.0) {
            for (int j=0; j<da.length; ++j) {
                da[j] /= sum;
            }
        }
        return maxIndex;
    }

    /**
     * Adds the FORMAT field for a sample to the VCF record for the current
     * marker.   If either of the specified posterior allele probabilities
     * does not sum to 1.0, it will be normalized to sum to 1.0.
     * @param alProbs1 the posterior allele probabilities for the individual's
     * first allele
     * @param alProbs2 the posterior allele probabilities for the individual's
     * second allele
     * @throws IllegalArgumentException if
     * {@code alProbs1.length != this.marker().nAlleles()}
     * @throws IllegalArgumentException if
     * {@code alProbs2.length != this.marker().nAlleles()}
     * @throws IllegalArgumentException if any element of the specified
     * array is not a finite non-negative number
     * @throws IllegalStateException if {@code this.marker() == null}
     * @throws NullPointerException if
     * {@code alProbs1 == null || alProbs2 == null}
     */
    public void addSampleData(double[] alProbs1, double[] alProbs2) {
        if (marker==null) {
            throw new IllegalStateException();
        }
        Arrays.fill(gt3Probs, 0.0);
        allele1 = maxIndex(alProbs1, nAlleles);
        allele2 = maxIndex(alProbs2, nAlleles);
        dose[0] = alProbs1[0] + alProbs2[0];
        gtProbs[0] = alProbs1[0] * alProbs2[0];
        gt3Probs[0] = gtProbs[0];
        int gt = 1;
        for (int a2=1; a2<alProbs1.length; ++a2) {
            dose[a2] = alProbs1[a2] + alProbs2[a2];
            for (int a1=0; a1<=a2; ++a1) {
                double gtProb = alProbs1[a1]*alProbs2[a2];
                if (a1!=a2) {
                    gtProb += alProbs1[a2]*alProbs2[a1];
                }
                gtProbs[gt++] = gtProb;
                gt3Probs[(a1==0) ? 1 : 2] += gtProb;
            }
        }
        addToCumAlleleProbs(dose);
        r2Est.addSampleData(gt3Probs);
        boolean isPhased = true;
        appendFormatData(isPhased);
    }

    private void addToCumAlleleProbs(double[] dose) {
        for (int j=0; j<dose.length; ++j) {
            cumAlleleProbs[j] += dose[j];
        }
    }

    private void appendFormatData(boolean isPhased) {
        sb.append(Const.tab);
        sb.append(allele1);
        sb.append(isPhased ? Const.phasedSep : Const.unphasedSep);
        sb.append(allele2);
        if (printDS) {
            for (int j=1; j<nAlleles; ++j) {
                sb.append( (j==1) ? Const.colon : Const.comma );
                sb.append(DS_VALS[(int) Math.rint(100*dose[j])]);
            }
        }
        if (printGP) {
            for (int j=0; j<gtProbs.length; ++j) {
                sb.append(j==0? Const.colon : Const.comma);
                sb.append(DS_VALS[(int) Math.rint(100*gtProbs[j])]);
            }
        }
    }

    /**
     * Prints the current VCF record for the current marker to the specified
     * {@code PrintWriter}.  If the FORMAT field contains a DS or GP subfield,
     * the INFO field will include the AR2 (allele r2), DR2 (dose r2), and
     * AF (ALT allele frequency) subfields.  Invocation of this method has
     * no effect if {@code this.reset()} has not previously been invoked.
     * @param out the {@code PrintWriter} to which the VCF record will be
     * printed
     * @param isImputed {@code true} if the printed VCF record will
     * have an IMP flag in the INFO field and {@code false} otherwise
     * @throws NullPointerException if {@code out == null}
     */
    public void writeRec(PrintWriter out, boolean isImputed) {
        if (marker!=null) {
            printMarker(marker, out);
            out.print(Const.tab);
            out.print(Const.MISSING_DATA_CHAR);     // QUAL
            out.print(Const.tab);
            out.print("PASS");                      // FILTER
            out.print(Const.tab);
            printInfo(out, isImputed);              // INFO
            out.print(Const.tab);
            out.print(format(printDS, printGP));    // FORMAT
            out.println(sb);
        }
    }

    private void printInfo(PrintWriter out, boolean isImputed) {
        if (printDS || printGP) {
            out.print("AR2=");
            out.print(R2_VALS[(int) Math.rint(100*r2Est.allelicR2())]);
            out.print(";DR2=");
            out.print(R2_VALS[(int) Math.rint(100*r2Est.doseR2())]);
            for (int j=1; j<nAlleles; ++j) {
                out.print( (j==1) ? ";AF=" : Const.comma);
                out.print(formatProb(cumAlleleProbs[j]/(2*r2Est.nGenotypes())));
            }
            if (isImputed) {
                out.print(";IMP");
            }
        }
        else {
            out.print(Const.MISSING_DATA_CHAR);
        }
    }

    private static String format(boolean printDS, boolean printGP) {
        if (printDS) {
            return printGP ? "GT:DS:GP" : "GT:DS";
        }
        else {
            return "GT";
        }
    }

    private static void printMarker(Marker marker, PrintWriter out) {
        out.print(marker.chrom());
        out.print(Const.tab);
        out.print(marker.pos());
        int nIds = marker.nIds();
        if (nIds==0) {
            out.print(Const.tab);
            out.print(Const.MISSING_DATA_CHAR);
        }
        else {
            for (int j=0; j<nIds; ++j) {
                out.print(j==0 ? Const.tab : Const.semicolon);
                out.print(marker.id(j));
            }
        }
        int nAlleles = marker.nAlleles();
        if (nAlleles==1) {
            out.print(Const.tab);
            out.print(marker.allele(0));
            out.print(Const.tab);
            out.print(Const.MISSING_DATA_CHAR);
        }
        else {
            for (int j=0; j<nAlleles; ++j) {
                out.print(j<2 ? Const.tab : Const.comma);
                out.print(marker.allele(j));
            }
        }
    }

    private static String formatProb(double p) {
        if (p>=0 && p <= 0.5) {
            return new BigDecimal(p).round(MC2).toString();
        }
        else if (p <= 1.0) {
            return new BigDecimal(p-1.0).round(MC2).add(ONE).toString();
        }
        else {
            throw new IllegalArgumentException(String.valueOf(p));
        }
    }
}
