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

import beagleutil.IntInterval;
import blbutil.Const;

/**
 * <p>Class {@code HapSegment} represents a marker interval
 * for a haplotype.
 * </p>
 * Instances of class {@code HapSegment} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class HapSegment implements Comparable<HapSegment>, IntInterval {

    private final int hap;
    private final int start;
    private final int end;

    /**
     * Constructs a new {@code HapSegment} instance.
     * @param hap the haplotype index
     * @param start the start marker index (inclusive)
     * @param end the end marker index (inclusive)
     * @throws IllegalArgumentException if {@code start > end}
     */
    public HapSegment(int hap, int start, int end) {
        if (start > end) {
            throw new IllegalArgumentException(String.valueOf(start));
        }
        this.hap = hap;
        this.start = start;
        this.end = end;
    }

    /**
     * Returns the first haplotype index.
     * @return the first haplotype index
     */
    public int hap() {
        return hap;
    }

    /**
     * Returns the start marker index (inclusive).
     * @return the start marker index (inclusive)
     */
    @Override
    public int start() {
        return start;
    }

    /**
     * Returns the end marker index (inclusive).
     * @return the end marker index (inclusive)
     */
    @Override
    public int end() {
        return end;
    }

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(70);
        sb.append(hap);
        sb.append(Const.tab);
        sb.append(start);
        sb.append(Const.tab);
        sb.append(end);
        return sb.toString();
    }

    /**
     * <p>Returns the hash code value for this object. The hash code is defined
     * by the following calculation:
     * </p>
     * <pre>
     *  int hash = 5;
     *  hash = 89 * hash + this.hap();
     *  hash = 89 * hash + this.start();
     *  hash = 89 * hash + this.end();
     </pre>
     * @return the hash code value for this object
     */
    @Override
    public int hashCode() {
        int hash = 5;
        hash = 89*hash + this.hap;
        hash = 89*hash + this.start;
        hash = 89*hash + this.end;
        return hash;
    }

    /**
     * Compares the specified object with this {@code HapSegment} for
     * equality. Returns {@code true} if the specified object is a
     * {@code HapSegment} instance and if this {@code HapSegment} is
     * equal to the specified {@code HapSegment}, and returns
     * {@code false}  otherwise.  Two {@code HapSegment}  instances
     * are equal if they have equal haplotype indices,
     * equal starting marker indices, and equal ending marker indices.
     * @param o the reference object with which to compare.
     * @return {@code true} if the specified object is an
     * {@code HapSegment} instance and if this {@code HapSegment} is
     * equal to the specified {@code HapSegment}
     */
    @Override
    public boolean equals(Object o) {
        if (o==null) {
            return false;
        }
        if (getClass()!=o.getClass()) {
            return false;
        }
        final HapSegment other=(HapSegment) o;
        return (this.hap==other.hap && this.start==other.start
                && this.end==other.end);
    }

    /**
     * Compares this object with the specified object for order.  Returns a
     * negative integer, zero, or a positive integer as this object is less
     * than, equal to, or greater than the specified object.
     * {@code HapSegment} instances are ordered first by
     * {@code this.start()}, then by {@code this.end()},
     * and finally by {@code this.hap()}.
     * @param hs the {@code HapSegment} to be compared
     * @return a negative integer, zero, or a positive integer as this
     * {@code HapSegment} is less than, equal to, or greater than the
     * specified {@code HapSegment}
     * @throws NullPointerException if {@code o == null}
     */
    @Override
    public int compareTo(HapSegment hs) {
        if (this.start != hs.start) {
            return (this.start < hs.start) ? -1 : 1;
        }
        else if (this.end != hs.end) {
            return (this.end < hs.end) ? -1 : 1;
        }
        if (this.hap != hs.hap) {
            return (this.hap < hs.hap) ? -1 : 1;
        }
        return 0;
    }
}
