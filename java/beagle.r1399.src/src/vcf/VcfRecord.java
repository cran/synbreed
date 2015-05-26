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
import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 * <p>Class {@code VcfRecord} represents a VCF record.
 * </p>
 * <p>Instances of class {@code VcfRecord} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class VcfRecord {

    private final String vcfRecord;
    private final VcfHeader vcfHeader;
    private final Marker marker;
    private final String qual;
    private final String filter;
    private final String info;
    private final String format;
    private final String[] sampleData;

    private final double qualityScore;

    private final String[] failedFilters;
    private final boolean filtersApplied;
    private final String[] infoFields;

    private final String[] formatFields;
    private final String[][] sampleFormatFields;
    private final boolean hasGTFormat;
    private final byte[][] alleles;
    private final boolean[] isPhased;

    private final Map<String, Integer> formatMap;

    /**
     * Returns the VCF genotype index for the specified pair of alleles.
     * @param a1 the first allele.
     * @param a2 the second allele.
     * @return the VCF genotype index for the specified pair of alleles.
     * @throws IllegalArgumentException if {@code a1<0 || a2<0}.
     */
    public static int gtIndex(int a1, int a2) {
        if (a1 < 0) {
            throw new IllegalArgumentException("a1<0: " + a1);
        }
        if (a2 < 0) {
            throw new IllegalArgumentException("a2<0: " + a2);
        } else if (a1 < a2) {
            return (a2 * (a2 + 1)) / 2 + a1;
        } else {
            return (a1 * (a1 + 1)) / 2 + a2;
        }
    }

    /**
     * Creates a new {@code VcfRecord} instance.
     *
     * @param vcfRecord a VCF version 4.1 record.
     * @param vcfHeader meta-information lines and header line for the
     * specified VCF record.
     *
     * @throws IllegalArgumentException if format error is
     * detected in any fixed field or in any non-excluded sample field.
     * @throws IllegalArgumentException if there are not
     * {@code vcfHeader.nHeaderFields()} tab-delimited fields in the
     * specified VCF record.
     * @throws NullPointerException if
     * {@code vcfRecord==null || vcfHeader==null}.
     */
    public VcfRecord(String vcfRecord, VcfHeader vcfHeader) {
        this.vcfRecord = vcfRecord;
        if (vcfHeader==null) {
            throw new NullPointerException("vcfHeader==null");
        }
        String[] fields = getAndCheckFields(vcfHeader, vcfRecord);
        this.vcfHeader = vcfHeader;
        this.marker = new Marker(vcfRecord);
        this.qualityScore = fromPhred(quality(fields[5]));
        this.qual = fields[5];
        this.filter = fields[6];
        this.info = fields[7];
        this.format = (fields.length > 8) ? fields[8] : "";
        this.filtersApplied = fields[6].equals(Const.MISSING_DATA_STRING)==false;
        this.failedFilters = failedFilters(fields[6]);
        this.infoFields = StringUtil.getFields(fields[7], Const.semicolon);
        this.formatFields = fields.length > 8 ? formats(fields[8]) : new String[0];
        this.formatMap = formatToIndexMap(vcfHeader, vcfRecord, formatFields);

        this.sampleData = Arrays.copyOfRange(fields, Math.min(fields.length, 9),
                fields.length);
        this.hasGTFormat = formatMap.containsKey("GT");
        int n = vcfHeader.samples().nSamples();
        this.alleles = new byte[2][n];
        this.isPhased = new boolean[n];
        this.sampleFormatFields = new String[n][formatFields.length];
        storePerSampleData(sampleData);
    }

    private static String[] getAndCheckFields(VcfHeader vcfHeader,
            String vcfRecord) {
        String[] fields = StringUtil.getFields(vcfRecord, Const.tab);
        if (vcfHeader.nHeaderFields() != fields.length) {
            File f = vcfHeader.file();
            String src = "File source: " + (f!=null ? f : "stdin or unknown");
            String s = "Header line has " + vcfHeader.nHeaderFields()
                    + " fields, but data line has " + fields.length + " fields"
                    + Const.nl + src + Const.nl + Arrays.toString(fields);
            throw new IllegalArgumentException(s);
        }
        return fields;
    }

    /**
     * Return {@code true} if all characters in the specified
     * string are letters or digits and returns {@code false} otherwise.
     * @param s a string.
     * @return {@code true} if all characters in the specified
     * string are letters or digits and returns {@code false} otherwise.
     */
    public static boolean isAlphanumeric(String s) {
        for (int j=0, n=s.length(); j<n; ++j) {
            if (Character.isLetterOrDigit(s.charAt(j))==false) {
                return false;
            }
        }
        return true;
    }

    /**
     * Returns {@code -10.0*Math.log10(d)}.
     * @param d a double to convert to the Phred scale.
     * @return {@code -10.0*Math.log10(d)}
     */
    public static double toPhred(double d) {
        return -10.0*Math.log10(d);
    }

   /**
     * Returns {@code Math.pow(10.0, -phred/10.0)}.
     * @param phred a double to convert from the Phred scale.
     * @return {@code Math.pow(10.0, -phred/10.0)}.
     */
    public static double fromPhred(double phred) {
        return Math.pow(10.0, -phred/10.0);
    }

    private double quality(String quality) {
        if (quality.equals(Const.MISSING_DATA_STRING)) {
            return Double.NaN;
        }
        else {
            return Double.parseDouble(quality);
        }
    }

    private String[] failedFilters(String filters) {
        if (filters.isEmpty()) {
            String s = "missing FILTER field: " + vcfRecord;
            throw new IllegalArgumentException(s);
        }
        if (filters.equals(Const.MISSING_DATA_STRING) || filters.equals("PASS")) {
            return new String[0];
        }
        String[] fields = StringUtil.getFields(filters, Const.semicolon);
        for (String f : fields) {
            if (f.isEmpty()) {
                String s = "missing filter in filter list: " + filters;
                throw new IllegalArgumentException(s);
            }
        }
        return fields;
    }

    private String[] formats(String formats) {
        if (formats.equals(Const.MISSING_DATA_STRING)) {
            return new String[0];
        }
        String[] fields =  StringUtil.getFields(formats, Const.colon);
        for (String f : fields) {
            if (f.isEmpty()) {
                String s = "missing format in format list: " + vcfRecord;
                throw new IllegalArgumentException(s);
            }
//            //  Commented-out alpha-numeric check to avoid throwing an
//            //    exception when FORMAT field code is not alphanumeric.
//            if (isAlphanumeric(f)==false) {
//                 String s = "format must be alphanumeric (" + f + "): " + vcfRecord;
//                 throw new IllegalArgumentException(s);
//            }
        }
        return fields;
    }

    private static Map<String, Integer> formatToIndexMap(VcfHeader vcfHeader,
            String vcfRecord, String[] formatFields) {
        if (vcfHeader.samples().nSamples()==0) {
            return Collections.emptyMap();
        }
        Map<String, Integer> map =
                new HashMap<>(formatFields.length);
        for (int j=0; j<formatFields.length; ++j) {
            map.put(formatFields[j], j);
        }
        if (map.containsKey("GT") && map.get("GT")!=0) {
            String s = "GT format is not first format: " + vcfRecord;
            throw new IllegalArgumentException(s);
        }
        return map;
    }

    private void storePerSampleData(String[] sampleFields) {
    /* Only GT field (optional first data field) is checked for validity */
        int filteredIndex = 0;
        for (int j=0; j<sampleFields.length; ++j) {
            if (vcfHeader.filter(j)==false) {
                String[] fields = parseSampleField(sampleFields, j);
                if (hasGTFormat) {
                    String gt = fields[0];
                    int sepIndex = separatorIndex(gt);
                    alleles[0][filteredIndex] = allele(gt.substring(0, sepIndex));
                    alleles[1][filteredIndex] = allele(gt.substring(sepIndex+1));
                    isPhased[filteredIndex]
                            = (gt.charAt(sepIndex)==Const.phasedSep);
                }
                else {
                    alleles[0][filteredIndex] = -1;
                    alleles[1][filteredIndex] = -1;
                    isPhased[filteredIndex] = false;
                }
                sampleFormatFields[filteredIndex++] = fields;
            }
        }
        assert filteredIndex==sampleFormatFields.length;
    }

    private String[] parseSampleField(String[] sampleFields, int sampleIndex) {
        String sampleField = sampleFields[sampleIndex];
        if (sampleField.isEmpty()) {
            String s = "Missing data for sample " +
                    vcfHeader.samples().id(sampleIndex) + ": " + vcfRecord;
            throw new IllegalArgumentException(s);
        }
        if (sampleField.equals(Const.MISSING_DATA_STRING)) {
            String[] fields = new String[formatFields.length];
            Arrays.fill(fields, Const.MISSING_DATA_STRING);
            if (hasGTFormat) {
                fields[0] = "./.";
            }
            return fields;
        }
        else {
            String[] fields = StringUtil.getFields(sampleField, Const.colon);
            for (String f : fields) {
                if (f.isEmpty()) {
                    String s = "empty sub-field for sample "
                            + vcfHeader.samples().id(sampleIndex)
                            + ": " + sampleField;
                    throw new IllegalArgumentException(s);
                }
            }
            if (fields.length < formatFields.length) {
                String[] newFields = Arrays.copyOf(fields, formatFields.length);
                for (int k=fields.length; k<newFields.length; ++k) {
                    newFields[k] = Const.MISSING_DATA_STRING;
                }
                fields = newFields;
            }
            if (fields.length > formatFields.length) {
                String s = "Expected at most " + formatFields.length
                        + " sub-fields for sample "
                        + vcfHeader.samples().id(sampleIndex) + ": " + vcfRecord;
                throw new IllegalArgumentException(s);
            }
            return fields;
        }
    }

    /**
     * Returns the index of the genotype separator;
     */
    private int separatorIndex(String gt) {
        int index = gt.indexOf(Const.unphasedSep);
        if (index == -1) {
            index = gt.indexOf(Const.phasedSep);
            if (index== -1) {
                String s = "missing genotype separator ("
                        + gt + "): " + vcfRecord;
                throw new IllegalArgumentException(s);
            }
        }
        return index;
    }

    private byte allele(String allele) {
        if (allele.isEmpty()) {
            String s = "Missing a sample allele: " + vcfRecord;
            throw new IllegalArgumentException(s);
        }
        if (allele.equals(Const.MISSING_DATA_STRING)) {
            return -1;
        }
        int a = Integer.parseInt(allele);
        if (a < 0) {
            String s = "allele cannot be negative (" + a + "): " + vcfRecord;
            throw new IllegalArgumentException(s);
        }
        if (a >= marker.nAlleles()) {
            String s = "allele " + a + " is not defined: " + vcfRecord;
            throw new IllegalArgumentException(s);
        }
        if (a > Byte.MAX_VALUE) {
            String s = "Marker cannot have more than " + Byte.MAX_VALUE
                    + " alternate alleles: " + vcfRecord;
            throw new IllegalArgumentException(s);
        }
        return (byte) a;
    }

    /**
     * Returns the QUAL field.
     * @return the QUAL field.
     */
    public String qual() {
        return qual;
    }

    /**
     * Returns the probability that the asserted presence or absence of
     * alternate bases is incorrect.  If ALT is ”.” (no variant) then this is
     * P(variant), and if ALT is not ”.” then this P(no variant).
     *
     * @return the probability that the asserted presence or absence of
     * alternate bases is incorrect.
     */
    public double qualityScore() {
        return qualityScore;
    }

    /**
     * Returns {@code true} if filters were applied, and returns
     * {@code false} if the FILTER field has missing data.
     *
     * @return {@code true} if filters were applied, and returns
     * {@code false} if the FILTER field has missing data.
     */
    public boolean filtersApplied() {
        return filtersApplied;
    }

    /**
     * Returns the FILTER field.
     * @return the FILTER field.
     */
    public String filter() {
        return filter;
    }

    /**
     * Returns the number of failed filters.
     * @return the number of failed filters.
     */
    public int nFailedFilters() {
        return failedFilters.length;
    }

    /**
     * Returns the specified failed filter code.
     *
     * @param index a failed filter index.
     * @return the specified failed filter code.
     *
     * @throws IndexOutOfBoundsException if
     * {@code index<0 || index>=this.nFailedFilters()}.
     */
    public String failedFilters(int index) {
        return failedFilters[index];
    }

    /**
     * Returns the INFO field.
     * @return the INFO field.
     */
    public String info() {
        return info;
    }

    /**
     * Returns the number of INFO sub-fields.
     * @return the number of INFO sub-fields.
     */
    public int nInfoFields() {
        return infoFields.length;
    }

    /**
     * Returns the specified INFO sub-field.
     * @param subfield an INFO sub-field index.
     * @return the specified INFO sub-field.
     *
     * @throws IndexOutOfBoundsException if
     * {@code index<0 || index>=this.nInfoFields()}.
     */
    public String infoField(int subfield) {
        return infoFields[subfield];
    }

    /**
     * Returns the FORMAT field.  Returns the empty string ("") if the FORMAT
     * field is missing.
     * @return the FORMAT field.
     */
    public String format() {
        return format;
    }

    /**
     * Returns the number of FORMAT subfields.
     * @return the number of FORMAT subfields.
     */
    public int nFormatFields() {
        return formatFields.length;
    }

    /**
     * Returns the specified FORMAT subfield.
     * @param index a FORMAT subfield index.
     * @return the specified FORMAT subfield.
     *
     * @throws IndexOutOfBoundsException if
     * {@code index<0 || index>=this.nFormatFields()}.
     */
    public String formatField(int index) {
        if (formatFields==null) {
            throw new IllegalArgumentException("No format exists");
        }
        return formatFields[index];
    }

    /**
     * Returns {@code true} if the specified FORMAT subfield is
     * present, and returns {@code false} otherwise.
     * @param formatCode a FORMAT sub-field code.
     * @return {@code true} if the specified FORMAT subfield is
     * present, and returns {@code false} otherwise.
     */
    public boolean hasFormat(String formatCode) {
        return formatMap.get(formatCode)!=null;
    }

    /**
     * Returns the index of the specified FORMAT subfield code if the
     * specified subfield is defined for the VCF record and returns -1
     * otherwise.
     * @param formatCode the string format code.
     * @return the index of the specified FORMAT subfield code if the
     * specified subfield is defined for the VCF record and returns -1
     * otherwise.
     */
    public int formatIndex(String formatCode) {
        Integer index = formatMap.get(formatCode);
        return (index==null) ? -1 : index.intValue();
    }

    /**
     * Returns the specified sample allele.
     * @param sample a sample index.
     * @param allele an allele index.
     * @return the specified sample allele.
     *
     * @throws IndexOutOfBoundsException if
     * {@code sample<0 || sample>=this.nSamples()}.
     * @throws IndexOutOfBoundsException if
     * {@code allele<0 || allele>=2}.
     */
    public byte gt(int sample, int allele) {
        return alleles[allele][sample];
    }

    /**
     * Returns {@code true} if the genotype for the specified sample is
     * phased and returns {@code false} otherwise.
     * @param sample a sample index.
     * @return  {@code true} if the genotype for the specified sample is
     * phased and returns {@code false} if the genotype is unphased.
     *
     * @throws IndexOutOfBoundsException if
     * {@code sample<0 || sample>=this.nSamples()}.
     */
    public boolean isPhased(int sample) {
        return isPhased[sample];
    }

    /**
     * Returns the data for the specified sample.
     * @param sample a sample index
     * @return the data for the specified sample.
     *
     * @throws IndexOutOfBoundsException if
     * {@code sample<0 || sample>=this.nSamples()}.
     */
    public String sampleData(int sample) {
        return sampleData[sample];
    }

    /**
     * Returns the specified FORMAT subfield data for the specified sample.
     * @param formatCode a FORMAT subfield code.
     * @param sample a sample index
     * @return the specified FORMAT subfield data for the specified sample.
     *
     * @throws IllegalArgumentException if
     * {@code this.hasFormat(formatCode)==false}.
     * @throws IndexOutOfBoundsException if
     * {@code sample<0 || sample>=this.nSamples()}.
     */
    public String sampleFormatData(String formatCode, int sample) {
        Integer formatIndex = formatMap.get(formatCode);
        if (formatIndex==null) {
            String s = "missing format data: " + formatCode;
            throw new IllegalArgumentException(s);
        }
        return sampleFormatFields[sample][formatIndex];
    }

    /**
     * Returns the specified FORMAT subfield data for the specified sample.
     * @param subfield a FORMAT subfield index.
     * @param sample a sample index.
     * @return the specified FORMAT subfield data for the specified sample.
     *
     * @throws IndexOutOfBoundsException if
     * {@code subfield<0 || subfield>=this.nFormatFields()}.
     * @throws IndexOutOfBoundsException if
     * {@code sampleIndex<0 || sampleIndex>=this.nSamples()}.
     */
    public String sampleFormatData(int subfield, int sample) {
        return sampleFormatFields[sample][subfield];
    }

    /**
     * Returns an array of length {@code this.nSamples()}
     * containing the specified FORMAT subfield data for each sample.  The
     * {@code k}-th element of the array is the specified FORMAT subfield data
     * for the {@code k}-th sample.
     * @param formatCode a format-field code.
     * @return an array of length {@code this.nSamples()}
     * containing the specified FORMAT subfield data for each sample.
     *
     * @throws IllegalArgumentException if
     * {@code this.hasFormat(formatCode)==false}.
     */
    public String[] formatData(String formatCode) {
        Integer formatIndex = formatMap.get(formatCode);
        if (formatIndex==null) {
            String s = "missing format data: " + formatCode;
            throw new IllegalArgumentException(s);
        }
        return formatData(formatIndex);
    }

   /**
     * Returns an array of length {@code this.nSamples()}
     * containing the specified FORMAT subfield data for each sample.  The
     * {@code k}-th element of the returned array is the specified FORMAT
     * subfield data for the {@code k}-th sample.
     * @param subfield a FORMAT subfield index.
     * @return an array of length {@code this.nSamples()}
     * containing the specified FORMAT subfield data for each sample.
     *
     * @throws IndexOutOfBoundsException if
     * {@code subfield<0 || subfield>=this.nFormatFields()}.
     */
    public String[] formatData(int subfield) {
        String[] sa = new String[sampleFormatFields.length];
        for (int j=0; j<sa.length; ++j) {
            sa[j] = sampleFormatFields[j][subfield];
        }
        return sa;
    }

    /**
     * Returns the samples. The returned samples are the filtered samples
     * after all sample exclusions.
     *
     * @return the samples.
     */
    public Samples samples() {
        return vcfHeader.samples();
    }

    /**
     * Returns the number of samples.  The number of samples is the
     * number of filtered samples after all sample exclusions.
     *
     * @return the number of samples.
     */
    public int nSamples() {
        return vcfHeader.samples().nSamples();
    }

    /**
     * Returns the VCF header line and meta-information lines.
     * @return the VCF header line and meta-information lines.
     */
    public VcfHeader vcfHeader() {
        return vcfHeader;
    }

    /**
     * Returns the marker.
     * @return the marker.
     */
    public Marker marker() {
        return marker;
    }

    /**
     * Returns the VCF record.
     * @return the VCF record.
     */
    @Override
    public String toString() {
        return vcfRecord;
    }
}
