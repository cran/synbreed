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
import blbutil.Const;
import blbutil.FileUtil;
import haplotype.HapPair;
import haplotype.BasicSampleHapPairs;
import haplotype.SampleHapPairs;
import java.io.File;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.List;
import main.GenotypeValues;
import main.GprobsStatistics;

/**
 * Class {@code VcfWriter} contains static methods for writing data in
 * VCF 4.1 format.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class VcfWriter {

    private static final String PASS = "PASS";
    private static final DecimalFormat df3 = new DecimalFormat("#.###");

    private static final String fileformat = "##fileformat=VCFv4.1";
    private static final String afInfo = "##INFO=<ID=AF,Number=A,Type=Float,"
            + "Description=\"Estimated Allele Frequencies\">";
    private static final String ar2Info = "##INFO=<ID=AR2,Number=1,Type=Float,"
            + "Description=\"Allelic R-Squared: estimated correlation between "
            + "most probable ALT dose and true ALT dose\">";
    private static final String dr2Info = "##INFO=<ID=DR2,Number=1,Type=Float,"
            + "Description=\"Dosage R-Squared: estimated correlation between "
            + "estimated ALT dose [P(RA) + 2*P(AA)] and true ALT dose\">";

    private static final String gtFormat = "##FORMAT=<ID=GT,Number=1,Type=String,"
            + "Description=\"Genotype\">";
    private static final String dsFormat = "##FORMAT=<ID=DS,Number=1,Type=Float,"
            +"Description=\"estimated ALT dose [P(RA) + P(AA)]\">";
    private static final String glFormat = "##FORMAT=<ID=GL,Number=G,Type=Float,"
            + "Description=\"Log10-scaled Genotype Likelihood\">";
    private static final String gpFormat = "##FORMAT=<ID=GP,Number=G,Type=Float,"
            + "Description=\"Estimated Genotype Probability\">";

    private static final String shortChromPrefix= "#CHROM" + Const.tab + "POS"
            + Const.tab + "ID" + Const.tab + "REF" + Const.tab + "ALT"
            + Const.tab + "QUAL" + Const.tab + "FILTER" + Const.tab + "INFO";

    private static final String longChromPrefix =
            shortChromPrefix + Const.tab + "FORMAT";


    private VcfWriter() {
        // private constructor prevents instantiation
    }

    /**
     * Writes the specified data to the specified VCF file.
     * @param source a description of the data source, or {@code null} if
     * no description is available.
     * @param haps the sample haplotype pairs.
     * @param vcfFile the file to which VCF output will be written.
     * @throws NullPointerException if {@code haps==null || vcfFile==null}.
     */
    public static void write(String source, SampleHapPairs haps,
            File vcfFile) {
        try (PrintWriter out=FileUtil.bgzipPrintWriter(vcfFile)) {
            writeMetaLinesGT(haps.samples().ids(), source, out);
            appendRecords(haps, 0, haps.nMarkers(), out);
        }
    }

    /**
     * Writes the specified data to the specified VCF file.  This method
     * permits an individual to have multiple haplotype pairs.
     * @param source a description of the data source, or {@code null} if
     * no description is available.
     * @param hapPairList a list of haplotype pairs.
     * @param vcfFile the file to which VCF output will be written.
     * @throws NullPointerException if {@code haps==null || vcfFile==null}.
     */
    public static void write(String source, List<HapPair> hapPairList,
            File vcfFile) {
        try (PrintWriter out=FileUtil.bgzipPrintWriter(vcfFile)) {
            String[] sampleIds = sampleIds(hapPairList);
            Markers markers = BasicSampleHapPairs.checkAndExtractMarkers(hapPairList);
            writeMetaLinesGT(sampleIds, source, out);
            for (int j=0; j<markers.nMarkers(); ++j) {
                printFixedFieldsGT(markers.marker(j), out);
                for (int hp=0; hp<sampleIds.length; ++hp) {
                    out.print(Const.tab);
                    out.print(hapPairList.get(hp).allele1(j));
                    out.print(Const.phasedSep);
                    out.print(hapPairList.get(hp).allele2(j));
                }
                out.println();
            }
        }
    }

    private static String[] sampleIds(List<HapPair> hapPairList) {
        String[] sa = new String[hapPairList.size()];
        for (int j=0; j<sa.length; ++j) {
            sa[j] = SampleIds.instance().id(hapPairList.get(j).idIndex());
        }
        return sa;
    }

    /**
     * Writes VCF meta-information lines to the specified {@code PrintWriter}.
     * The meta-information lines assume that the only FORMAT field is the GT
     * field.
     * @param sampleIds the sample identifiers.
     * @param source a description of the data source, or {@code null} if
     * no description is to be printed.
     * @param out the {@code PrintWriter} to which VCF meta-information
     * lines will be written.
     * @throws NullPointerException if {@code sampleIds==null}, if
     * {@code out==null}, or if {@code sampleIds[j]==null} for any
     * {@code 0<=j<sampleIds.length}.
     */
    public static void writeMetaLinesGT(String[] sampleIds, String source,
            PrintWriter out) {
        boolean printGT = true;
        boolean printGP = false;
        boolean printGL = false;
        writeMetaLines(sampleIds, source, printGT, printGP, printGL, out);
    }

    /**
     * Writes VCF meta-information lines to the specified {@code PrintWriter}.
     * @param sampleIds the sample identifiers.
     * @param source a description of the data source, or {@code null} if
     * no description is to be printed.
     * @param printGT {@code true} if there is a GT FORMAT field and
     * {@code false} otherwise.
     * @param printGP {@code true} if there is a GP FORMAT field and
     * {@code false} otherwise.
     * @param printGL {@code true} if there is a GL FORMAT field and
     * {@code false} otherwise.
     * @param out the {@code PrintWriter} to which VCF meta-information lines
     * will be written.
     * @throws NullPointerException if {@code sampleIds==null}, if
     * {@code out==null}, or if {@code sampleIds[j]==null} for any
     * {@code 0<=j<sampleIds.length}.
     */
    public static void writeMetaLines(String[] sampleIds, String source,
            boolean printGT, boolean printGP, boolean printGL, PrintWriter out) {
        out.print(fileformat);
        out.print(Const.nl);
        out.print("##filedate=");
        out.print(now());
        out.print(Const.nl);
        if (source != null) {
            out.print("##source=\"");
            out.print(source);
            out.println("\"");
        }
        if (printGP) {
            out.println(afInfo);
            out.println(ar2Info);
            out.println(dr2Info);
        }
        if (printGT) {
            out.println(gtFormat);
        }
        if (printGL) {
            out.println(glFormat);
        }
        if (printGP) {
            out.println(dsFormat);
            out.println(gpFormat);
        }
        out.print(longChromPrefix);
        for (String id : sampleIds) {
            if (id==null) {
                throw new NullPointerException("id==null");
            }
            out.print(Const.tab);
            out.print(id);
        }
        out.println();
    }

    private static String now() {
        String dateFormat = "yyyyMMdd";
        Calendar cal = Calendar.getInstance();
        SimpleDateFormat sdf = new SimpleDateFormat(dateFormat);
        return sdf.format(cal.getTime());
    }

    /**
     * Writes the specified haplotypes and posterior genotype
     * probabilities as VCF records to the specified {@code PrintWriter}.
     * @param haps the sample haplotype pairs.
     * @param gv the scaled sample posterior genotype probabilities.
     * @param start the starting marker index (inclusive).
     * @param end the ending marker index (exclusive).
     * @param out the {@code PrintWriter} to which VCF records will
     * be written.
     *
     * @throws IllegalArgumentException if
     * {@code haps.markers().equals(gv.markers()==false}.
     * @throws IndexOutOfBoundsException if
     * {@code start<0 || start>end || end>haps.nMarkers()}.
     * @throws NullPointerException if
     * {@code haps==null || gv==null || out==null}.
     */
    public static void appendRecords(SampleHapPairs haps,
            GenotypeValues gv, int start, int end, PrintWriter out) {
        if (start > end) {
            throw new IllegalArgumentException("start=" + start + " end=" + end);
        }
        if (haps.markers().equals(gv.markers())==false) {
            throw new IllegalArgumentException("inconsistent markers");
        }
        float[] sumAndAltDose = new float[2];
        for (int marker=start; marker<end; ++marker) {
            printFixedFields(gv, marker, out);
            for (int hp=0, n=haps.nSamples(); hp<n; ++hp) {
                out.print(Const.tab);
                out.print(haps.allele1(marker, hp));
                out.print(Const.phasedSep);
                out.print(haps.allele2(marker, hp));
                int sampleIdIndex = haps.idIndex(hp);
                int sampleIndex = gv.samples().index(sampleIdIndex);
                if (sampleIndex < 0) {
                    throw new IllegalArgumentException("inconsistent samples");
                }
                int nGenotypes = gv.marker(marker).nGenotypes();
                sumAndAltDose(gv, marker, sampleIndex, sumAndAltDose);
                float sum = sumAndAltDose[0];
                float altDoseSum = sumAndAltDose[1];
                if (sum==0.0f) {
                    out.print(Const.colon);
                    out.print(Const.MISSING_DATA_CHAR);
                    out.print(Const.colon);
                    out.print(Const.MISSING_DATA_CHAR);
                }
                else {
                    out.print(Const.colon);
                    out.print(df3.format(altDoseSum/sum));
                    for (int gt=0; gt<nGenotypes; ++gt) {
                        out.print(gt==0 ? Const.colon : Const.comma);
                        double v = gv.value(marker, sampleIndex, gt)/sum;
                        out.print(df3.format(v));
                    }
                }
            }
            out.println();
        }
    }

    /**
     * Writes the specified haplotypes as VCF records to the specified
     * {@code PrintWriter}.
     * @param haps the sample haplotype pairs.
     * @param start the starting marker index (inclusive).
     * @param end the ending marker index (exclusive).
     * @param out the {@code PrintWriter} to which VCF records will
     * be written.
     *
     * @throws IndexOutOfBoundsException if
     * {@code start<0 || start>end || end>haps.nMarkers()}.
     * @throws NullPointerException if
     * {@code haps==null || out==null}.
     */
    public static void appendRecords(SampleHapPairs haps,
            int start, int end, PrintWriter out) {
        if (start > end) {
            throw new IllegalArgumentException("start=" + start + " end=" + end);
        }
        for (int marker=start; marker<end; ++marker) {
            printFixedFieldsGT(haps.marker(marker), out);
            for (int hp=0, n=haps.nSamples(); hp<n; ++hp) {
                out.print(Const.tab);
                out.print(haps.allele1(marker, hp));
                out.print(Const.phasedSep);
                out.print(haps.allele2(marker, hp));
            }
            out.println();
        }
    }


    private static void sumAndAltDose(GenotypeValues gv,
            int marker, int sample, float[] sumAndAltAlleleSum) {
        float sum = 0.0f;
        float altSum = 0.0f;
        int nAlleles = gv.marker(marker).nAlleles();
        int gt = 0;
        for (int a2=0; a2<nAlleles; ++a2) {
            for (int a1=0; a1<a2; ++a1) {
                float f = gv.value(marker, sample, gt++);
                sum += f;
                altSum += (a1==0) ? f : (2.0f)*f;
            }
            float f = gv.value(marker, sample, gt++);
            sum += f;
            altSum += (a2==0) ? 0.0 : (2.0f)*f;
        }
        assert gt==gv.marker(marker).nGenotypes();
        sumAndAltAlleleSum[0] = sum;
        sumAndAltAlleleSum[1] = altSum;
    }

    /**
     * Prints the fixed fields for the specified marker to the specified
     * {@code PrintWriter}.  The fixed fields include only the GT FORMAT
     * subfield.
     *
     * @param marker the marker whose fixed fields will be written.
     * @param out the {@code PrintWriter} to which VCF records will
     * be written.
     *
     * @throws NullPointerException if {@code marker==null || out==null}.
     */
    public static void printFixedFieldsGT(Marker marker, PrintWriter out) {
        out.print(marker);
        out.print(Const.tab);
        out.print(Const.MISSING_DATA_CHAR); // QUAL
        out.print(Const.tab);
        out.print(PASS);                    // FILTER
        out.print(Const.tab);
        out.print(Const.MISSING_DATA_CHAR); // INFO
        out.print(Const.tab);
        out.print("GT");
    }

    private static void printFixedFields(GenotypeValues gv, int marker,
            PrintWriter out) {
        GprobsStatistics gpm = new GprobsStatistics(gv, marker);
        float[] alleleFreq = gpm.alleleFreq();
        out.print(gv.marker(marker));
        out.print(Const.tab);
        out.print(Const.MISSING_DATA_CHAR); // QUAL
        out.print(Const.tab);
        out.print(PASS);                    // FILTER
        out.print(Const.tab);
        out.print("AR2=");                  // INFO
        out.print(df3.format(gpm.allelicR2()));
        out.print(";DR2=");
        out.print(df3.format(gpm.doseR2()));
        for (int j=1; j<alleleFreq.length; ++j) {
            out.print( (j==1) ? ";AF=" : Const.comma);
            out.print(df3.format(alleleFreq[j]));
        }
        out.print(Const.tab);
        out.print("GT:DS:GP");
    }
}
