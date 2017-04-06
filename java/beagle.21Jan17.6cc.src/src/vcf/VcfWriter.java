/*__LICENSE__*/
package vcf;

import blbutil.Const;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import main.AlleleProbs;
import main.GenotypeValues;

/**
 * <p>Class {@code VcfWriter} contains static methods for writing data in
 * VCF 4.2 format.
 * </p>
 * <p>Instances of class {@code VcfWriter} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class VcfWriter {

    private static final String fileformat = "##fileformat=VCFv4.2";

    private static final String afInfo = "##INFO=<ID=AF,Number=A,Type=Float,"
            + "Description=\"Estimated ALT Allele Frequencies\">";
    private static final String ar2Info = "##INFO=<ID=AR2,Number=1,Type=Float,"
            + "Description=\"Allelic R-Squared: estimated squared correlation between "
            + "most probable REF dose and true REF dose\">";
    private static final String dr2Info = "##INFO=<ID=DR2,Number=1,Type=Float,"
            + "Description=\"Dosage R-Squared: estimated squared correlation between "
            + "estimated REF dose [P(RA) + 2*P(RR)] and true REF dose\">";
    private static final String impInfo = "##INFO=<ID=IMP,Number=0,Type=Flag,"
            + "Description=\"Imputed marker\">";

    private static final String gtFormat = "##FORMAT=<ID=GT,Number=1,Type=String,"
            + "Description=\"Genotype\">";
    private static final String dsFormat = "##FORMAT=<ID=DS,Number=A,Type=Float,"
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
     * Writes VCF meta-information lines and header line to the specified
     * {@code PrintWriter}. Only one FORMAT subfield, the GT subfield,
     * is described in the meta-information lines.
     * @param sampleIds the sample identifiers
     * @param source a description of the data source, or {@code null} if
     * no description is to be printed
     * @param out the {@code PrintWriter} to which VCF meta-information
     * lines will be written
     * @throws NullPointerException if {@code out == null}
     * @throws NullPointerException if
     * {@code sampleIds == null}, or if {@code sampleIds[j] == null} for any
     * {@code j} satisfying {@code (0 <= j && j < <sampleIds.length)}
     */
    public static void writeMetaLinesGT(String[] sampleIds, String source,
            PrintWriter out) {
        boolean printGT = true;
        boolean printGP = false;
        boolean printGL = false;
        writeMetaLines(sampleIds, source, printGT, printGP, printGL, out);
    }

    /**
     * Writes VCF meta-information lines and header line to the specified
     * {@code PrintWriter}.
     * @param sampleIds the sample identifiers
     * @param source a description of the data source, or {@code null} if
     * no description is to be printed
     * @param printGT {@code true} if the meta-information lines
     * will describe the GT FORMAT subfield and {@code false} otherwise
     * @param printGP {@code true} if the meta-information lines
     * will describe the GP FORMAT subfield and {@code false} otherwise
     * @param printGL {@code true} if the meta-information lines
     * will describe the GL FORMAT subfield and {@code false} otherwise
     * @param out the {@code PrintWriter} to which VCF meta-information lines
     * will be written.
     * @throws NullPointerException if {@code out == null}
     * @throws NullPointerException if
     * {@code sampleIds == null}, or if {@code sampleIds[j] == null} for any
     * {@code j} satisfying {@code (0 <= j && j < sampleIds.length)}
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
            out.println(impInfo);
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
     * Writes the specified genotype data  as VCF records to the specified
     * {@code PrintWriter}.
     * @param gv the scaled sample posterior genotype probabilities
     * @param start the starting marker index (inclusive)
     * @param end the ending marker index (exclusive)
     * @param out the {@code PrintWriter} to which VCF records will
     * be written.
     *
     * @throws IllegalArgumentException if
     * {@code haps.markers().equals(gv.markers()) == false}
     * @throws IndexOutOfBoundsException if
     * {@code (start < 0 || start > end || end > haps.nMarkers())}
     * @throws NullPointerException if
     * {@code (gv == null || out == null)}
     */
    public static void appendRecords(GenotypeValues gv, int start, int end,
            PrintWriter out) {
        if (start > end) {
            throw new IllegalArgumentException("start=" + start + " end=" + end);
        }
        boolean printDS = true;
        boolean printGP = true;
        boolean isImputed = false;
        VcfRecBuilder vrb = new VcfRecBuilder(12*gv.nSamples());
        for (int m=start; m<end; ++m) {
            Marker marker = gv.marker(m);
            vrb.reset(marker, printDS, printGP);
            double[] gprobs = new double[marker.nGenotypes()];
            for (int s=0, n=gv.nSamples(); s<n; ++s) {
                double sum = 0.0;
                for (int gt=0; gt<gprobs.length; ++gt) {
                    gprobs[gt] = gv.value(m, s, gt);
                    sum += gprobs[gt];
                }
                for (int gt=0; gt<gprobs.length; ++gt) {
                    gprobs[gt] /= sum;
                }
                vrb.addSampleData(gprobs);
            }
            vrb.writeRec(out, isImputed);
        }
    }

    /**
     * Writes the data in alProbs for markers with index between
     * {@code start} (inclusive) and {@code end} (exclusive) to the specified
     * {@code PrintWriter}.
     * @param alProbs the estimated haplotype allele probabilities
     * @param isImputed an array of length {@code alProbs.nMarkers()}
     * whose {@code j}-th element is {@code true} if the corresponding
     * marker is imputed, and {@code false} otherwise
     * @param start the starting marker index (inclusive)
     * @param end the ending marker index (exclusive)
     * @param printDS {@code true} if the DS field should be printed, and
     * {@code false} otherwise
     * @param printGP {@code true} if the GP field should be printed, and
     * {@code false} otherwise
     * @param out the {@code PrintWriter} to which VCF records will be written
     * @throws IllegalArgumentException if
     * {@code isImputed.length != alProbs.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code (start < 0 || start > end || end > alProbs.nMarkers())}
     * @throws NullPointerException if
     * {@code alProbs == null || isImputed == null || out == null}
     */
    public static void appendRecords(AlleleProbs alProbs, boolean[] isImputed,
            int start, int end, boolean printDS, boolean printGP,
            PrintWriter out) {
        if (isImputed.length != alProbs.nMarkers()) {
            throw new IllegalArgumentException("inconsistent data");
        }
        if (start > end) {
            throw new IllegalArgumentException("start=" + start + " end=" + end);
        }
        VcfRecBuilder vrb = new VcfRecBuilder(4*alProbs.nSamples());
        for (int m=start; m<end; ++m) {
            Marker marker = alProbs.marker(m);
            vrb.reset(marker, printDS, printGP);
            double[] a1 = new double[marker.nAlleles()];
            double[] a2 = new double[marker.nAlleles()];
            for (int sample=0, n=alProbs.nSamples(); sample<n; ++sample) {
                for (int j=0; j<a1.length; ++j) {
                    a1[j] = alProbs.alProb1(m, sample, j);
                    a2[j] = alProbs.alProb2(m, sample, j);
                }
                vrb.addSampleData(a1, a2);
            }
            vrb.writeRec(out, isImputed[m]);
        }
    }

    /**
     * Writes the data in alProbs for markers with index between
     * {@code start} (inclusive) and {@code end} (exclusive) to the specified
     * {@code PrintWriter}.
     * @param alProbs the estimated haplotype allele probabilities
     * @param isImputed an array of length {@code alProbs.nMarkers()}
     * whose {@code j}-th element is {@code true} if the corresponding
     * marker is imputed, and {@code false} otherwise
     * @param start the starting marker index (inclusive)
     * @param end the ending marker index (exclusive)
     * @param printDS {@code true} if the DS field should be printed, and
     * {@code false} otherwise
     * @param printGP {@code true} if the GP field should be printed, and
     * {@code false} otherwise
     * @param out the {@code PrintWriter} to which VCF records will be written
     * @throws IllegalArgumentException if
     * {@code isImputed.length != alProbs.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code (start < 0 || start > end || end > alProbs.nMarkers())}
     * @throws NullPointerException if
     * {@code alProbs == null || isImputed == null || out == null}
     */
    public static void blockedAppendRecords(AlleleProbs alProbs,
            boolean[] isImputed, int start, int end,
            boolean printDS, boolean printGP, PrintWriter out) {
        if (isImputed.length != alProbs.nMarkers()) {
            throw new IllegalArgumentException("inconsistent data");
        }
        if (start > end) {
            throw new IllegalArgumentException("start=" + start + " end=" + end);
        }
        int nSamples = alProbs.nSamples();
        int size = 8;
        int initCapacity = 4*alProbs.nSamples();
        VcfRecBuilder[] recBuilders = new VcfRecBuilder[size];
        double[][] a1 = new double[size][];
        double[][] a2 = new double[size][];
        for (int j=0; j<recBuilders.length; ++j) {
            recBuilders[j] = new VcfRecBuilder(initCapacity);
        }
        for (int ii=start; ii<end; ii+=size) {
            int ni = Math.min(end, ii+size);
            for (int i=ii; i<ni; ++i) {
                Marker marker = alProbs.marker(i);
                recBuilders[i-ii].reset(marker, printDS, printGP);
                a1[i-ii] = new double[marker.nAlleles()];
                a2[i-ii] = new double[marker.nAlleles()];
                for (int jj=0; jj<nSamples; jj+=size) {
                    int nj = Math.min(jj+size, nSamples);
                    for (int j=jj; j<nj; ++j) { // j = sample
                        for (int k=0; k<a1[i-ii].length; ++k) {
                            a1[i-ii][k] = alProbs.alProb1(i, j, k);
                            a2[i-ii][k] = alProbs.alProb2(i, j, k);
                        }
                        recBuilders[i-ii].addSampleData(a1[i-ii], a2[i-ii]);
                    }
                }
            }
            for (int i=ii; i<ni; ++i) {
                recBuilders[i-ii].writeRec(out, isImputed[i]);
            }
        }
    }

    /**
     * Prints the first 9 VCF record fields for the specified marker to
     * the specified {@code PrintWriter}.  Only one VCF FORMAT subfield,
     * the GT subfield, is printed.
     *
     * @param marker a marker
     * @param out the {@code PrintWriter} to which the first 9 VCF record
     * fields will be written
     *
     * @throws NullPointerException if {@code marker == null || out == null}
     */
    public static void printFixedFieldsGT(Marker marker, PrintWriter out) {
        out.print(marker);
        out.print(Const.tab);
        out.print(Const.MISSING_DATA_CHAR); // QUAL
        out.print(Const.tab);
        out.print("PASS");                    // FILTER
        out.print(Const.tab);
        out.print(Const.MISSING_DATA_CHAR); // INFO
        out.print(Const.tab);
        out.print("GT");
    }
}
