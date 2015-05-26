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
import blbutil.FileIterator;
import blbutil.Filter;
import blbutil.StringUtil;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * <p>Class {@code VcfHeader} represents the VCF meta-information lines and
 * the VCF header line that precede the first VCF record.
 * </p>
 * <p>Instances of class {@code VcfHeader} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class VcfHeader  {

    private static final String SHORT_HEADER_PREFIX= "#CHROM" + Const.tab + "POS"
            + Const.tab + "ID" + Const.tab + "REF" + Const.tab + "ALT"
            + Const.tab + "QUAL" + Const.tab + "FILTER" + Const.tab + "INFO";
    /**
     * A string equal to the first nine tab-delimited fields of a VCF header
     * line.
     */
    public static final String HEADER_PREFIX =
            SHORT_HEADER_PREFIX + Const.tab + "FORMAT";

    private static final int nShortHeaderFields
            = StringUtil.countFields(SHORT_HEADER_PREFIX, Const.tab);
    private static final int nLongHeaderFields
            = StringUtil.countFields(HEADER_PREFIX, Const.tab);

    private final File file;   // null if source is standard input or unknown
    private final VcfMetaInfo[] metaInfoLines;
    private final String headerLine;
    private final int nHeaderFields;
    private final Samples samples;
    private final boolean[] filter;

    /**
     * Constructs a new {@code VcfHeader} object from the VCF
     * meta-information lines and the VCF header line returned by the
     * specified {@code FileIterator}.
     * @param it an iterator that returns lines of a VCF file.
     * @param sampleFilter a sample filter.
     *
     * @throws IllegalArgumentException if any of the meta-information and
     * header lines returned by the specified {@code FileIterator}
     * do not conform to the VCF specification.
     *
     * @throws NullPointerException if {@code it==null || sampleFilter==null}.
     */
    public VcfHeader(FileIterator<String> it, Filter<String> sampleFilter) {
        List<VcfMetaInfo> metaInfo = new ArrayList<>(20);
        String candidateHeader = null;
        while (it.hasNext() && candidateHeader==null) {
            String line = it.next().trim();
            if (line.startsWith(VcfMetaInfo.PREFIX)) {
                metaInfo.add(new VcfMetaInfo(line));
            }
            else if (line.startsWith(SHORT_HEADER_PREFIX)) {
                candidateHeader = line;
            }
            else {
                headerError(line, it.file());
            }
        }
        assert candidateHeader != null;
        String[] headerFields = StringUtil.getFields(candidateHeader, Const.tab);
        int nSamples = nSamples(headerFields, candidateHeader, it.file());

        this.file = it.file();
        this.metaInfoLines = metaInfo.toArray(new VcfMetaInfo[0]);
        this.headerLine = candidateHeader;
        this.nHeaderFields = headerFields.length;
        this.filter = new boolean[nSamples];
        this.samples = getSamplesAndSetFilter(headerFields, sampleFilter, filter);
    }

    private static boolean headerError(String line, File file) {
        String src = "File source: " + (file!=null ? file : "stdin or unknown");
        if (line.startsWith("#")==false) {
            String s = "Missing line (#CHROM ...) after meta-information lines"
                    + Const.nl + src + Const.nl + line;
            throw new IllegalArgumentException(s);
        }
        String[] fields = StringUtil.getFields(line, nShortHeaderFields);
        fields = Arrays.copyOf(fields, nShortHeaderFields);
        String[] headerFields = StringUtil.getFields(SHORT_HEADER_PREFIX,
                Const.tab);
        if (Arrays.equals(fields, headerFields)) {
            String s = "Header line is white-space delimited, but not tab "
                    + "delimited" + Const.nl + src + Const.nl + line;
            throw new IllegalArgumentException(s);
        }
        else {
            String s = "Error in meta-information line or header line"
                    + Const.nl + src + Const.nl + line;
            throw new IllegalArgumentException(s);
        }
    }

    private static int nSamples(String[] headerFields, String header, File file) {
        if (headerFields.length==nShortHeaderFields) {
            return 0;
        }
        else if (headerFields.length==nLongHeaderFields) {
            if (headerFields[nLongHeaderFields-1].equals("FORMAT")==false) {
                String src = "File source: "
                        + (file!=null ? file : "stdin or unknown");
                String s = "Ninth field of header line is not \"FORMAT\""
                        + Const.nl + src + Const.nl + header;
                throw new IllegalArgumentException(s);
            }
            return 0;
        }
        else {
            return (headerFields.length - nLongHeaderFields);
        }
    }

    private static Samples getSamplesAndSetFilter(String[] headerFields,
            Filter<String> sampleFilter, boolean[] filter) {
        if (filter.length==0) {
            return new Samples(new int[0]);
        }
        else {
            List<String> filteredIds = new ArrayList<>(filter.length);
            for (int j=0; j<filter.length; ++j) {
                String id = headerFields[nLongHeaderFields + j];
                if (sampleFilter.accept(id)) {
                    filteredIds.add(id);
                }
                else {
                    filter[j] = true;
                }
            }
            return Samples.fromIds(filteredIds.toArray(new String[0]));
        }
    }

    /**
     * Returns the file from which data are read, or returns
     * {@code null} if the source is standard input or if
     * the file source is unknown.
     * @return the file from which data are read, or returns
     * {@code null} if the source is standard input or if
     * the file source is unknown.
     */
    public File file() {
        return file;
    }

    /**
     * Returns the number of VCF meta-information lines. VCF meta-information
     * lines are lines that precede the VCF header line.  The first
     * two characters of a VCF meta-information line must be "##".
     *
     * @return the number of VCF meta-information lines.
     */
     public int nMetaInfoLines() {
         return metaInfoLines.length;
     }

    /**
      * Returns the specified VCF meta-information line.

      * @param index a VCF meta-information line index.
      * @return the specified VCF meta-information line.
      *
      * @throws IndexOutOfBoundsException if
      * {@code index < 0 || index >= this.nMetaInfoLines()}.
      */
     public VcfMetaInfo metaInfoLine(int index) {
         return metaInfoLines[index];
     }

     /**
      * Returns the VCF header line.  The VCF header line begins with "#CHROM".
      * @return the VCF header line.
      */
     public String headerLine() {
         return headerLine;
     }

     /**
      * Returns {@code true} if the specified sample in the VCF
      * header line is excluded, and returns {@code false} otherwise.
      *
      * @param index a sample index for the list of unfiltered samples.
      *
      * @return {@code true} if the specified sample in the VCF
      * header line is excluded, and returns {@code false} otherwise.
      *
      * @throws IndexOutOfBoundsException if
      * {@code index < 0 || index >= this.nUnfilteredSamples()}.
      */
     public boolean filter(int index) {
         return filter[index];
     }

     /**
      * Returns the number of fields in the VCF header line before sample
      * exclusions.
      * @return the number of fields in the VCF header line before sample
      * exclusions.
      */
     public int nHeaderFields() {
         return nHeaderFields;
     }

     /**
      * Returns the number of samples before sample exclusions.
      * @return the number of samples before sample exclusions.
      */
     public int nUnfilteredSamples() {
         return filter.length;
     }

    /**
     * Return the filtered list of samples after all sample exclusions.
     * @return the filtered list of samples after all sample exclusions.
     */
    public Samples samples() {
        return samples;
    }

    /**
     * Returns {@code this.sample().ids()}.
     * @return {@code this.sample().ids()}.
     */
    public String[] sampleIds() {
        return samples.ids();
    }

    /**
     * Returns the VCF meta-information lines and the VCF header line used to
     * construct {@code this}.
     * @return the VCF meta-information lines and the VCF header line used to
     * construct {@code this}.
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(400);
        for (int j=0; j<metaInfoLines.length; ++j) {
            sb.append(metaInfoLines[j]);
            sb.append(Const.nl);
        }
        sb.append(headerLine);
        sb.append(Const.nl);
        return sb.toString();
    }
}
