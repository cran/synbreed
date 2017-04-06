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
package ibd;

import haplotype.SampleHapPairs;
import vcf.Marker;

/**
 * <p>Class {@code Haplotype} represents a haplotype segment.
 * </p>
 * Instances of class {@code Haplotype} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Haplotype {

    private final int start;    // inclusive
    private final int end;      // exclusive
    private final SampleHapPairs haps;
    private final int hapIndex;

    /**
     * Constructs a new {@code Haplotype} instance. The haplotype will
     * include all markers in the specified {@code SampleHapPairs} parameter.
     *
     * @param haps sample haplotype pairs
     * @param hap a haplotype index
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap >= haps.nHaps()}
     * @throws NullPointerException if {@code haps == null}
     */
    public Haplotype(SampleHapPairs haps, int hap) {
        this(haps, hap, 0, haps.nMarkers());
    }


    /**
     * Constructs a new {@code Haplotype} instance.
     *
     * @param haps sample haplotype pairs
     * @param hap a haplotype index
     * @param start the starting marker index for the haplotype segment
     * (inclusive)
     * @param end the ending marker index for the haplotype segment (exclusive)
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap >= haps.nHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code start < 0 || start > end || end > haps.nMarkers()}
     * @throws NullPointerException if {@code haps == null}
     */
    public Haplotype(SampleHapPairs haps, int hap, int start, int end) {
        if (start < 0 || start > end || end > haps.nMarkers()) {
            String s = "start=" + start + " end=" + end + " haps.nMarkers()="
                    + haps.nMarkers();
            throw new IndexOutOfBoundsException(s);
        }
        if (hap < 0 || hap >= haps.nHaps()) {
            throw new IllegalArgumentException("hapIndex=" + hap);
        }
        this.start = start;
        this.end = end;
        this.haps = haps;
        this.hapIndex = hap;
    }

    /**
     * Returns a new {@code Haplotype} instance that is
     * obtained by restricting this haplotype to the specified marker interval.
     * @param start the starting marker index for the haplotype segment
     * (inclusive)
     * @param end the ending marker index for the haplotype segment (exclusive)
     * @return the restricted haplotype segment
     * @throws IndexOutOfBoundsException if
     * {@code start < 0 || start > end || end > this.length()}
     */
    public Haplotype restrict(int start, int end) {
        if (start < 0 || start > end || end > this.length()) {
            String s = "start=" + start + " end=" + end + " this.length()="
                    + this.length();
            throw new IndexOutOfBoundsException(s);
        }
        int newStart = this.start + start;
        int newEnd = this.start + end;
        return new Haplotype(haps, hapIndex, newStart, newEnd);
    }

    /**
     * Returns the number of alleles in this haplotype segment.
     * @return the number of alleles in this haplotype segment
     */
    public int length() {
        return end - start;
    }

    /**
     * Returns the index of the haplotype.
     * @return the index of the haplotype
     */
    public int hapIndex() {
        return hapIndex;
    }

    /**
     * Returns the sample haplotype pairs.
     * @return the sample haplotype pairs
     */
    public SampleHapPairs sampleHapPairs() {
        return haps;
    }

    /**
     * Returns the specified marker. The first marker on the haplotype segment
     * has index 0.
     * @param index a marker index
     * @return the specified marker
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.length()}
     */
    public Marker marker(int index) {
        int i = start + index;
        if (i < start || i >= end) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        return haps.marker(i);
    }

    /**
     * Returns the specified allele on the haplotype.  The first allele on
     * the haplotype segment has index 0.
     *
     * @param index a marker index
     * @return the specified allele on the haplotype
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.length()}
     */
    public int allele(int index) {
        int i = start + index;
        if (i < start || i >= end) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        return haps.allele(i, hapIndex);
    }

    /**
     * Compares the specified object with this {@code  Haplotype} for
     * equality.  Returns {@code  true} if the specified object
     * is a {@code  Haplotype} that represents the same haplotype segment
     * as {@code  this}, and returns {@code  false} otherwise.
     * @param obj the object to be compared for equality with this
     * {@code  Haplotype}
     * @return {@code  true} if the specified object is an {@code  Haplotype}
     * that represents the same haplotype segment as {@code  this}
     */
    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Haplotype other = (Haplotype) obj;
        int length = this.length();
        if (length != other.length()) {
            return false;
        }
        for (int j=0; j<length; ++j) {
            if (this.allele(j) != other.allele(j)) {
                return false;
            }
        }
        for (int j=0; j<length; ++j) {
            if (false==this.marker(j).equals(other.marker(j))) {
                return false;
            }
        }
        return true;
    }

    /**
     * <p>Returns a hash code value for the object.
     * </p>
     * <p>The hash code is defined by the following calculation:
     * </p>
     * <pre>
        int hash = 17;
        for (int j=0; j&lt;this.length(); ++j) {
            hash += 29 * hash + haps.allele(j, this.hapIndex());
            hash += 29 * hash + haps.marker(j).hashCode();
        }
     * </pre>
     * @return a hash code value for the object
     */
    @Override
    public int hashCode() {
        int hash = 17;
        for (int j = start; j<end; ++j) {
            hash += 29 * hash + haps.allele(j, hapIndex);
            hash += 29 * hash + haps.marker(j).hashCode();
        }
        return hash;
    }

    /**
     * Returns a string representation of {@code this}.  The
     * exact details of the representation are unspecified and
     * subject to change.
     *
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(3 * (end - start));
        sb.append('[');
        if (end > start) {
            sb.append(allele(0));
        }
        for (int j = 1, n = end - start; j < n; ++j) {
            sb.append(", ");
            sb.append(allele(j));
        }
        sb.append(']');
        return sb.toString();
    }

}
