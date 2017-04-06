/*
 * Copyright (C) 2015 Brian L. Browning
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
import blbutil.Utilities;
import java.io.File;
import java.util.Arrays;

/**
 * <p>Class {@code VcfRecGTParser} parses VCF records and extracts the GT format
 * field.
 * </p>
 * <p>Instances of class {@code VcfRecGTParser} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class VcfRecGTParser {

    private final VcfHeader vcfHeader;
    private final String vcfRec;
    private final Marker marker;
    private final int nSamples;

    private int character;   // points to delimiter at start of a field
    private int currentSample;
    private int unfilteredSample;

    private int allele1;
    private int allele2;
    private boolean isPhased;

    /**
     * Constructs a new {@code VcfRecGTParser} object from the specified VCF
     * record.
     * @param vcfHeader the VCF meta-information lines and header line
     * @param vcfRec the VCF record
     * @throws IllegalArgumentException if
     * {@code vcfHeader.nSamples() == 0}
     * @throws IllegalArgumentException if a format error is detected in the
     * {@code vcfRecord}
     * @throws NullPointerException if
     * {@code vcfHeader == null || vcfRec == null}
     */
    public VcfRecGTParser(VcfHeader vcfHeader, String vcfRec) {
        if (vcfHeader.nSamples()==0) {
            throw new IllegalArgumentException("nSamples==0");
        }
        this.vcfHeader = vcfHeader;
        this.vcfRec = vcfRec;
        this.marker = new BasicMarker(vcfRec);
        this.nSamples = vcfHeader.nSamples();

        this.unfilteredSample = -1;
        this.character = -1;
        this.currentSample = -1;
        skipFixedFields();
        nextSample();
    }

    private void skipFixedFields() {
        for (int j=0; j<8; ++j) {
            character = vcfRec.indexOf(Const.tab, character + 1);
        }
        if (vcfRec.startsWith("GT:", character + 1) == false
                && vcfRec.startsWith("GT\t", character + 1) == false) {
            throw new IllegalArgumentException("invalid VCF rec: " + vcfRec);
        }
        character = vcfRec.indexOf(Const.tab, character + 1);
    }

    /**
     * Returns the VCF meta-information lines and header line.
     * @return the VCF meta-information lines and header line
     */
    public VcfHeader vcfHeader() {
        return vcfHeader;
    }

    /**
     * Returns the VCF record that is being parsed.
     * @return the VCF record that is being parsed
     */
    public String vcfRecord() {
        return vcfRec;
    }

    /**
     * Returns the marker.
     * @return the marker
     */
    public Marker marker() {
        return marker;
    }

    /**
     * Returns the index of the current sample, or
     * -1 if {@code this.nextSample()} has not yet been invoked.
     * @return the index of the current sample, or
     * -1 if {@code this.nextSample()} has not yet been invoked
     */
    public int currentSample() {
        return currentSample;
    }

    /**
     * Returns the first allele of the genotype for the current sample.
     * @return the first allele of the genotype for the current sample
     */
    public int allele1() {
        return allele1;
    }

    /**
     * Returns the second allele of the genotype for the current sample.
     * @return the second allele of the genotype for the current sample
     */
    public int allele2() {
        return allele2;
    }

    /**
     * Returns {@code true} if the genotype for the current sample is phased,
     * and returns {@code false} otherwise.
     * @return {@code true} if the genotype for the current sample is phased
     */
    public boolean isPhased() {
        return isPhased;
    }

    /**
     * Returns the list of samples.
     * @return the list of samples
     */
    public Samples samples() {
        return vcfHeader.samples();
    }

    /**
     * Returns the number of samples.
     * @return the number of samples
     */
    public int nSamples() {
        return nSamples;
    }

    /**
     * Increases the current sample index by one.
     * @throws IndexOutOfBoundsException if
     * {@code (this.currentSample() + 1 == this.nSamples())} immediately prior
     * to the invocation of this method
     * @throws IllegalArgumentException if a format error is detected in
     * {@code this.vcfRecord()}
     */
    public void nextSample() {
        ++currentSample;
        if (currentSample == nSamples) {
            throw new IndexOutOfBoundsException(String.valueOf(currentSample));
        }
        int nextUnfilteredSample = vcfHeader.unfilteredSampleIndex(currentSample);
        while (++unfilteredSample < nextUnfilteredSample) {
            if (character == -1) {
                throwFieldCountError();
            }
            character = vcfRec.indexOf(Const.tab, character + 1);
        }
        if (character == -1) {
            throwFieldCountError();
        }
        int end1 = end1(vcfRec, character + 1);
        int end2 = end2(vcfRec, end1 + 1);
        this.allele1 = allele(vcfRec, marker.nAlleles(), character + 1, end1);
        this.allele2 = allele(vcfRec, marker.nAlleles(), end1 + 1, end2);
        this.isPhased = vcfRec.charAt(end1) != Const.unphasedSep;
        character = vcfRec.indexOf(Const.tab, end2);
    }

    /* returns exclusive end */
    private static int end1(String rec, int start) {
        if (start==rec.length()) {
            throwGTFormatError(rec, rec.length());
        }
        int index = start;
        while (index < rec.length()) {
            char c = rec.charAt(index);
            if (c == Const.unphasedSep || c == Const.phasedSep) {
                return index;
            }
            else if (c == Const.colon || c == Const.tab) {
                throwGTFormatError(rec, index+1);
            }
            ++index;
        }
        if (index==rec.length()) {
            throwGTFormatError(rec, rec.length());
        }
        return index;
    }

    /* returns exclusive end */
    private static int end2(String rec, int start) {
        int index = start;
        while (index < rec.length()) {
            char c = rec.charAt(index);
            if (c == Const.colon || c == Const.tab) {
                return index;
            }
            ++index;
        }
        return index;
    }

    private static int allele(String vcfRecord, int nAlleles, int start,
            int end) {
        if (start==end) {
            String s = "Missing sample allele: " + vcfRecord;
            throw new IllegalArgumentException(s);
        }
        int a;
        if (start + 1 == end) {
            char c = vcfRecord.charAt(start);
            switch (c) {
                case '.' : a = -1; break;
                case '0' : a = 0; break;
                case '1' : a = 1; break;
                case '2' : a = 2; break;
                case '3' : a = 3; break;
                case '4' : a = 4; break;
                case '5' : a = 5; break;
                case '6' : a = 6; break;
                case '7' : a = 7; break;
                case '8' : a = 8; break;
                case '9' : a = 9; break;
                default: String s = "invalid allele (" + c + "): "
                        + vcfRecord;
                        throw new IllegalArgumentException(s);
            }
        }
        else {
            a = Integer.parseInt(vcfRecord.substring(start, end));
            if (a < 0) {
                String s = "invalid allele (" + a + "): " + Const.nl + vcfRecord;
                throw new IllegalArgumentException(s);
            }
        }
        if (a >= nAlleles) {
            String s = "invalid allele (" + a + "): " + Const.nl + vcfRecord;
            throw new IllegalArgumentException(s);
        }
        return a;
    }

    private static void throwGTFormatError(String rec, int index) {
        StringBuilder sb = new StringBuilder(1000);
        sb.append("ERROR: Missing one or both alleles for a genotype:");
        sb.append(Const.nl);
        sb.append(rec.substring(0, index));
        sb.append(Const.nl);
        sb.append("Exiting Program");
        sb.append(Const.nl);
        Utilities.exit(sb.toString());
    }

    private void throwFieldCountError() {
        File f = vcfHeader.file();
        String[] fields = StringUtil.getFields(vcfRec, Const.tab);
        StringBuilder sb = new StringBuilder(1000);
        sb.append("VCF header line has ");
        sb.append(vcfHeader.nHeaderFields());
        sb.append(" fields, but data line has ");
        sb.append(fields.length);
        sb.append(" fields");
        sb.append(Const.nl);
        sb.append("File source: ");
        sb.append((f!=null ? f.toString() : "stdin"));
        sb.append(Const.nl);
        sb.append(Arrays.toString(fields));
        sb.append(Const.nl);
        Utilities.exit(sb.toString());
    }
}
