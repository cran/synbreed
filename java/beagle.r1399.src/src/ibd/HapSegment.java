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
 * for an ordered pair of haplotype indices.
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
     * Constructs a {@code HapSegment} instance.
     * @param hap the haplotype index.
     * @param start the start marker index (inclusive).
     * @param end the end marker index (inclusive).
     */
    public HapSegment(int hap, int start, int end) {
        this.hap = hap;
        this.start = start;
        this.end = end;
    }

    /**
     * Returns the first haplotype index.
     * @return the first haplotype index.
     */
    public int hap() {
        return hap;
    }

    /**
     * Returns the start marker index (inclusive).
     * @return the start marker index (inclusive).
     */
    @Override
    public int start() {
        return start;
    }

    /**
     * Returns the end marker index (inclusive).
     * @return the end marker index (inclusive).
     */
    @Override
    public int end() {
        return end;
    }

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}.
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
     * Compares this object with the specified object for order.  Returns a
     * negative integer, zero, or a positive integer as this object is less
     * than, equal to, or greater than the specified object.
     * {@code HapSegment} instances are ordered first by
     * {@code this.start()}, then by {@code this.end()},
     * and finally by {@code this.hap()}.
     * @param o the object to be compared
     * @return a negative integer, zero, or a positive integer as this object
     * is less than, equal to, or greater than the specified object.
     * @throws NullPointerException if {@code o==null}.
     */
    @Override
    public int compareTo(HapSegment o) {
        if (this.start != o.start) {
            return (this.start < o.start) ? -1 : 1;
        }
        else if (this.end != o.end) {
            return (this.end < o.end) ? -1 : 1;
        }
        if (this.hap != o.hap) {
            return (this.hap < o.hap) ? -1 : 1;
        }
        return 0;
    }
}
