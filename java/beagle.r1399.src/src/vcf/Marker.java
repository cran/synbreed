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

import beagleutil.ChromIds;
import blbutil.Const;
import blbutil.StringUtil;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * <p>Class {@code Marker} represents a marker.
 * </p>
 * <p>Instances of class {@code Marker} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Marker implements Comparable<Marker> {

    private static final int MIN_NUMBER_FIELDS = 8;

    private static final String[] EMPTY_ID_ARRAY = new String[0];
    private static final Map<String, String[]> allelesMap
            = new HashMap<>(24);

    private final int chromIndex;
    private final int pos;
    private final String[] ids;
    private final String id;
    private final String[] alleles;
    private final int nGenotypes;
    private final int end;

    /**
     * Constructs a new {@code Marker} instance from the specified
     * VCF record.
     * @param vcfRecord a VCF record.
     * @throws IllegalArgumentException if the specified VCF file
     * record has fewer than 8 tab-delimited fields, or if the
     * first 5 tab-delimited fields have an incorrect format.
     */
    @SuppressWarnings("RedundantStringConstructorCall")
    public Marker(String vcfRecord) {
        String[] fields = StringUtil.getFields(vcfRecord, Const.tab,
                MIN_NUMBER_FIELDS+1);
        if (fields.length < MIN_NUMBER_FIELDS) {
            String s = "VCF record does not contain at least "
                    + MIN_NUMBER_FIELDS + " tab-delimited fields: "
                    + vcfRecord;
            throw new IllegalArgumentException(s);
        }
        // Store minimal required data, not entire VCF record
        for (int j=0; j<5; ++j) {
            fields[j] = new String(fields[j]);
        }
        checkCHROM(fields[0], vcfRecord);
        checkPOS(fields[1], vcfRecord);
        String[] markerIds = checkID(fields[2], vcfRecord);
        checkREF(fields[3], vcfRecord);
        String[] altAlleles = checkALT(fields[4], vcfRecord);

        this.chromIndex = ChromIds.instance().indexOf(fields[0]);
        this.pos = Integer.parseInt(fields[1]);
        this.ids = markerIds;
        this.id = markerIds.length>0 ? markerIds[0] :
                (fields[0] + Const.colon + pos);
        this.alleles = alleles(fields[3], altAlleles);
        this.nGenotypes = alleles.length*(alleles.length+1)/2;
        this.end = extractEnd(fields[7]);
    }

    /**
     * Returns the marker obtained from the specified marker by changing
     * the marker's alleles to the alleles on the opposite chromosome strand.
     * @param marker a marker.
     * @return the marker obtained from the specified marker by changing
     * the marker's alleles to the alleles on the opposite chromosome strand.
     * @throws NullPointerException if {@code marker==null}.
     */
    public static Marker flipStrand(Marker marker) {
        return new Marker(marker);
    }

    /* Private constructor used by flipStrand(Marker) method */
    private Marker(Marker markerOnReverseStrand) {
        Marker m = markerOnReverseStrand;
        this.chromIndex = m.chromIndex;
        this.pos = m.pos;
        this.ids = m.ids;
        this.id = m.id;
        this.alleles = m.alleles.clone();
        for (int j=0; j<m.alleles.length; ++j) {
            if (alleles[j].charAt(0)!='<') {    // not a symbolic allele
                this.alleles[j] = flipAllele(m.alleles[j]);
            }
        }
        this.nGenotypes = m.nGenotypes;
        this.end = m.end;
    }

    private static String flipAllele(String allele) {
        char[] ca = new char[allele.length()];
        for (int j=0; j<ca.length; ++j) {
            ca[j] = flipBase(allele.charAt(j));
        }
        return new String(ca);
    }

    private static char flipBase(char c) {
        switch (c) {
            case 'A' : return 'T';
            case 'C' : return 'G';
            case 'G' : return 'C';
            case 'T' : return 'A';
            case 'N' : return 'N';
            case '*' : return '*';
            default: assert false; return 0;
        }
    }

    private static void checkCHROM(String chrom, String vcfRecord) {
        if (chrom.isEmpty() || chrom.equals(Const.MISSING_DATA_STRING)) {
            String s = "missing CHROM field: " + vcfRecord;
            throw new IllegalArgumentException(s);
        }
        for (int j=0, n=chrom.length(); j<n; ++j) {
            char c = chrom.charAt(j);
            if (c==Const.colon || Character.isWhitespace(c)) {
                String s = "invalid character in CHROM field ['" + c
                        + "']: " + vcfRecord;
                throw new IllegalArgumentException(s);
            }
        }
    }

    private static void checkPOS(String pos, String vcfRecord) {
        for (int j=0, n=pos.length(); j<n; ++j) {
            if (Character.isDigit(pos.charAt(j))==false) {
                String s = "invalid POS field [" + pos + "]: " + vcfRecord;
                throw new IllegalArgumentException(s);
            }
        }
    }

    private static String[] checkID(String id, String vcfRecord) {
        if (id.isEmpty()) {
            String s = "missing ID field: " + vcfRecord;
            throw new IllegalArgumentException(s);
        }
        if (id.equals(Const.MISSING_DATA_STRING)) {
            return EMPTY_ID_ARRAY;
        }
        String[] sa = StringUtil.getFields(id, Const.semicolon);
        for (String s : sa) {
            for (int j=0, n=s.length(); j<n; ++j) {
                char c = s.charAt(j);
                if (Character.isWhitespace(c)) {
                    String msg = "marker identifier (" + s
                            + ") contains white-space: " + vcfRecord;
                    throw new IllegalArgumentException(msg);
                }
            }
        }
        return sa;
    }

    private static void checkREF(String ref, String vcfRecord) {
        if (ref.isEmpty()) {
            String s = "missing REF field: " + vcfRecord;
            throw new IllegalArgumentException(s);
        }
        for  (int j=0, n=ref.length(); j<n; ++j) {
            char c = Character.toUpperCase(ref.charAt(j));
            if ((c=='A' || c=='C' || c=='G' || c=='T' || c=='N')==false) {
                String s = "REF allele [" + ref + "] is not a sequence"
                        + " of A, C, T, G, or N characters" + Const.nl
                        + vcfRecord;
                throw new IllegalArgumentException(s);
            }
        }
    }

    private static String[] checkALT(String alt, String vcfRecord) {
        if (alt.isEmpty()) {
            String s = "missing ALT field: " + vcfRecord;
            throw new IllegalArgumentException(s);
        }
        String[] altAlleles = EMPTY_ID_ARRAY;
        if (alt.equals(Const.MISSING_DATA_STRING)==false) {
            altAlleles = StringUtil.getFields(alt, Const.comma);
        }
        if (altAlleles.length >= Byte.MAX_VALUE - 1) {
            String s = "More than " + (Byte.MAX_VALUE - 1) + " ALT alleles: "
                    + vcfRecord;
            throw new IllegalArgumentException(s);
        }
        for (String s : altAlleles) {
            checkAltAllele(s, vcfRecord);
        }
        return altAlleles;
    }

    private static void checkAltAllele(String alt, String vcfRecord) {
        int n = alt.length();
        if (n >= 2 && alt.charAt(0)=='<' && alt.charAt(n-1)=='>') {
            for (int j=1; j<n-1; ++j) {
                char c = alt.charAt(j);
                if (Character.isWhitespace(c) || c==Const.comma || c=='<'
                        || c=='>') {
                    String s = "invalid allele (" + alt + "): " + vcfRecord;
                    throw new IllegalArgumentException(s);
                }
            }
        }
        else {
            for (int j=0; j<n; ++j) {
                char c = Character.toUpperCase(alt.charAt(j));
                if ((c=='A' || c=='C' || c=='G' || c=='T'
                        || c=='N' || c=='*')==false) {
                    String s = "ALT allele [" + alt + "] is not a sequence of"
                            + " A, C, T, G, N, or '*' characters" + Const.nl
                            + vcfRecord;
                    throw new IllegalArgumentException(s);
                }
            }
        }
    }

    private static String[] alleles(String ref, String[] altAlleles) {
        if (isSNV(ref, altAlleles)) {
            String key = ref;
            for (String a : altAlleles) {
                key += a;
            }
            String[] alleles = allelesMap.get(key);
            if (alleles==null) {
                alleles = createAllelesArray(ref, altAlleles);
                allelesMap.put(key, alleles);
            }
            return alleles;
        }
        else {
            return createAllelesArray(ref, altAlleles);
        }
    }

    private static boolean isSNV(String ref, String[] altAlleles) {
        if (ref.length()!=1) {
            return false;
        }
        for (String a : altAlleles) {
            if (a.length()!=1) {
                return false;
            }
        }
        return true;
    }


    private static String[] createAllelesArray(String ref, String[] altAlleles) {
        String[] alleles = new String[altAlleles.length + 1];
        alleles[0] = ref;
        System.arraycopy(altAlleles, 0, alleles, 1, altAlleles.length);
        return alleles;
    }

    /*
     * Returns value of first END key in the specified INFO field, or
     * returns -1 if there is no END key in INFO field.
     */
    private static int extractEnd(String info) {
        String[] fields = StringUtil.getFields(info, Const.semicolon);
        String key = "END=";
        for (String field : fields) {
            if (field.startsWith(key)) {
                String value = field.substring(4);
                for (int j=0, n=value.length(); j<n; ++j) {
                    char c = value.charAt(j);
                    if (Character.isDigit(c)==false) {
                        String s = "INFO END field has non-numeric value: "
                                + info;
                        throw new IllegalArgumentException(s);
                    }
                }
                return Integer.parseInt(value);
            }
        }
        return -1;
    }

    /**
     * Returns the chromosome.
     * @return the chromosome.
     */
    public String chrom() {
        return ChromIds.instance().id(chromIndex);
    }

    /**
     * Returns the chromosome index.
     * @return the chromosome index.
     */
    public int chromIndex() {
        return chromIndex;
    }

    /**
     * Returns the chromosome coordinate.
     * @return the chromosome coordinate.
     */
    public int pos() {
        return pos;
    }

    /**
     * Returns the number of marker identifiers.
     * @return the number of marker identifiers.
     */
    public int nIds() {
        return ids.length;
    }

    /**
     * Returns the specified marker identifier.
     * @param index a marker identifier index.
     * @return the specified marker identifier.
     *
     * @throws IndexOutOfBoundsException if
     * {@code index<0 || index>=this.nIds()}.
     */
    public String id(int index) {
        return ids[index];
    }

    /**
     * Returns the first marker identifier if there is at least
     * one identifier in the VCF record ID field, and returns
     * {@code this.chr() + ":" + this.pos()} otherwise.
     *
     * @return a marker identifier.
     */
    public String id() {
        return id;
    }

    /**
     * Returns the number of alleles for the marker, including the REF
     * allele.
     * @return the number of alleles for the marker, including the REF
     * allele.
     */
    public int nAlleles() {
        return alleles.length;
    }

    /**
     * Returns the number of distinct genotypes:
     * {@code this.nAlleles()*(1 + this.nAlleles())/2}.
     *
     * @return the number of distinct genotypes:
     * {@code this.nAlleles()*(1 + this.nAlleles())/2}.
     */
    public int nGenotypes() {
        return nGenotypes;
    }

    /**
     * Returns the specified allele.  The reference allele has index 0.
     * @param index an allele index.
     * @return the specified allele.
     *
     * @throws IndexOutOfBoundsException if
     * {@code index<0 || index>=this.nAlleles()}.
     */
    public String allele(int index) {
        return alleles[index];
    }

    /**
     * Returns the INFO END field, or -1 if there is no INFO END field.
     *
     * @return the INFO END field, or -1 if there is no INFO END field.
     */
    public int end() {
        return end;
    }

    /**
     * Returns a string equal to the first five tab-delimited fields
     * of a VCF record corresponding to this marker.
     *
     * @return a string equal to the first five tab-delimited fields
     * of a VCF record corresponding to this marker.
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(50);
        sb.append(chrom());
        sb.append(Const.tab);
        sb.append(pos);
        if (ids.length==0) {
            sb.append(Const.tab);
            sb.append(Const.MISSING_DATA_CHAR);
        }
        else {
            for (int j=0; j<ids.length; ++j) {
                sb.append(j==0 ? Const.tab : Const.semicolon);
                sb.append(ids[j]);
            }
        }
        if (alleles.length==1) {
            sb.append(Const.tab);
            sb.append(alleles[0]);
            sb.append(Const.tab);
            sb.append(Const.MISSING_DATA_CHAR);
        }
        else {
            for (int j=0; j<alleles.length; ++j) {
                sb.append(j<2 ? Const.tab : Const.comma);
                sb.append(alleles[j]);
            }
        }
        return sb.toString();
    }

    /**
     * <p>Returns the hash code value for this object. The hash code does not
     * depend on value of the VCF record ID field.
     * The hash code is defined by the following calculation:
     * </p>
     * <pre>
     *   int hash = 5;
     *   hash = 29 * hash + this.chromIndex();
     *   hash = 29 * hash + this.pos();
     *   for (int j=0, n=this.nAlleles(); j&le;n; ++j) {
     *       hash = 29 * hash + alleles[j].hashCode();
     *   }
     *   hash = 29 * hash + end;
     * </pre>
     *
     * @return the hash code value for this marker.
     */
    @Override
    public int hashCode() {
        int hash = 5;
        hash = 29 * hash + chromIndex;
        hash = 29 * hash + this.pos;
        for (int j=0; j<alleles.length; ++j) {
            hash = 29 * hash + alleles[j].hashCode();
        }
        hash = 29 * hash + end;
        return hash;
    }

    /**
     * Returns {@code true} if the specified object is a
     * {@code Marker} with the same chromosome,
     * position, allele lists, and INFO END field, and
     * returns {@code false} otherwise.  Equality does not
     * depend on value of the VCF record ID field.
     *
     * @param obj object to be compared with {@code this} for equality.
     *
     * @return {@code true} if the specified object is a
     * {@code Marker} with the same chromosome,
     * position, and allele lists, and INFO END field, and
     * {@code false} otherwise.
     */
    @Override
    public boolean equals(Object obj) {
        if (this==obj) {
            return true;
        }
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Marker other = (Marker) obj;
        if (this.chromIndex != other.chromIndex) {
            return false;
        }
        if (this.pos != other.pos) {
            return false;
        }
        if (!Arrays.equals(this.alleles, other.alleles)) {
            return false;
        }
        return this.end == other.end;
    }

    /**
     * Compares this marker with the specified marker
     * for order, and returns a negative integer, 0, or a positive integer
     * depending on whether this marker is less than, equal to,
     * or greater than the specified marker.  Comparison is
     * on chromosome index ({@code chromIndex()}), position,
     * allele identifier lists, and end value in that order.  Allele
     * identifier lists are compared for lexicographical order, and
     * alleles are compared using the {@code String compareTo()} method.
     *
     * @param other the {@code Marker} to be compared.
     * @return a negative integer, 0, or a positive integer
     * depending on whether this marker is less than, equal,
     * or greater than the specified marker.
     */
    @Override
    public int compareTo(Marker other) {
        if (this.chromIndex != other.chromIndex) {
            return (this.chromIndex < other.chromIndex) ? -1 : 1;
        }
        if (this.pos != other.pos) {
            return (this.pos < other.pos) ? -1 : 1;
        }
        int n = Math.min(this.alleles.length, other.alleles.length);
        for (int j=0; j<n; ++j) {
            int cmp = this.alleles[j].compareTo(other.alleles[j]);
            if (cmp != 0) {
                return cmp;
            }
        }
        if (this.alleles.length != other.alleles.length) {
            return (this.alleles.length < other.alleles.length) ? -1 : 1;
        }
        if (this.end != other.end) {
            return (this.end < other.end) ? -1 : 1;
        }
        return 0;
    }
}
