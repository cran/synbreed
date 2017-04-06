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

import blbutil.Const;
import blbutil.IntPair;
import java.text.DecimalFormat;
import vcf.Marker;

/**
 * <p>Class {@code IbdSegment} represents a pair of IBD haplotype segments.
 * </p>
 * <p>Instances of class {@code IbdSegment} are immutable.
 * </p>
 *
 * @author Brian L Browning {@code <browning@uw.edu>}
 */
public final class IbdSegment {

    private static final DecimalFormat df2 = new DecimalFormat("0.00");

    private final IntPair hapPair;
    private final Marker start;     // inclusive
    private final Marker end;       // inclusive
    private final float score;
    private final int startIndex;   // inclusive; -1 if missing
    private final int endIndex;     // inclusive; -1 if missing

    /**
     * Constructs an new {@code IbdSegment} instance from the specified data.
     * @param hapPair an ordered pair of haplotype indices
     * @param start the starting marker for the IBD segment (inclusive)
     * @param end the ending marker for the IBD segment (inclusive)
     * @param score the score for the IBD segment
     * @param startIndex the starting marker index (inclusive) or -1 if
     * the starting marker index is unknown
     * @param endIndex the ending marker index (inclusive) or -1 if
     * the ending marker index is unknown
     *
     * @throws IllegalArgumentException if
     * {@code hapPair.first() < 0 || hapPair.second() <= hapPair.first()}
     * @throws IllegalArgumentException if
     * {@code start.chromIndex() != end.chromIndex() || end.pos() < start.pos()}
     * @throws IllegalArgumentException if
     * {@code startIndex < -1 || endIndex < -1}
     * @throws IllegalArgumentException if {@code Float.isNaN(score) == true}
     * @throws NullPointerException if
     * {@code hapPair == null || start == null || end == null}
     */
    public IbdSegment(IntPair hapPair, Marker start, Marker end, float score,
            int startIndex, int endIndex) {
        checkArguments(hapPair, start, end, score, startIndex, endIndex);
        this.hapPair = hapPair;
        this.start = start;
        this.end = end;
        this.score = score;
        this.startIndex = startIndex;
        this.endIndex = endIndex;
    }

    private void checkArguments(IntPair hapPair, Marker start, Marker end,
            float score, int startIndex, int endIndex) {
        if (hapPair.first()<0 || hapPair.second()<=hapPair.first()) {
            throw new IllegalArgumentException(hapPair.toString());
        }
        if ( (start.chromIndex()!=end.chromIndex())
                || (end.pos() < start.pos()) ){
            String s = Const.nl + start + Const.nl + end;
            throw new IllegalArgumentException(s);
        }
        if (Float.isNaN(score)) {
           throw new IllegalArgumentException(String.valueOf(score));
        }
        if (startIndex < -1) {
            throw new IllegalArgumentException(String.valueOf(startIndex));
        }
        if (endIndex < -1) {
            throw new IllegalArgumentException(String.valueOf(endIndex));
        }
    }

    /**
     * Compares the specified object with this {@code IbdSegment} for
     * equality. Returns {@code true} if the specified object is an
     * {@code IbdSegment} instance and if this {@code IbdSegment} is
     * equal to the specified {@code IbdSegment}, and returns
     * {@code false}  otherwise.  Two {@code IbdSegment}  instances
     * are equal if they have equal ordered pairs of haplotype indices,
     * equal starting and ending markers, and equal scores.
     *
     * @param o the reference object with which to compare
     *
     * @return {@code true} if this {@code IbdSegment} is
     * equal to the specified object.
     */
    @Override
    public boolean equals(Object o) {
        if (this==o) {
            return true;
        }
        if ((o instanceof IbdSegment)==false) {
            return false;
        }
        IbdSegment other = (IbdSegment) o;
        if (this.hapPair.equals(other.hapPair)==false) {
            return false;
        }
        if (this.start.equals(other.start)==false) {
            return false;
        }
        if (this.end.equals(other.end)==false) {
            return false;
        }
        return Float.floatToIntBits(this.score)
                == Float.floatToIntBits(other.score);
    }

    /**
     * <p>Returns the hash code value for this object. The hash code does not
     * depend on the values of {@code this.startIndex()} or
     * {@code this.endIndex()}. The hash code is defined by the following
     * calculation:
     * </p>
     * <pre>
     *  int hash = 5;
     *  hash = 67 * hash + this.hapPair().hashCode();
     *  hash = 67 * hash + this.start().hashCode();
     *  hash = 67 * hash + this.end().hashCode();
     *  hash = 67 * hash + Float.floatToIntBits(this.score());
     </pre>
     * @return the hash code value for this object
     */
    @Override
    public int hashCode() {
        int hash = 5;
        hash = 67 * hash + this.hapPair.hashCode();
        hash = 67 * hash + this.start.hashCode();
        hash = 67 * hash + this.end.hashCode();
        hash = 67 * hash + Float.floatToIntBits(this.score);
        return hash;
    }

    /**
     * Returns the first haplotype index.
     * @return the first haplotype index
     */
    public int hap1() {
        return hapPair.first();
    }

    /**
     * Returns the second haplotype index.
     * @return the second haplotype index
     */
    public int hap2() {
        return hapPair.second();
    }

    /**
     * Returns the ordered pair of haplotype indices.
     * @return the ordered pair of haplotype indices
     */
    public IntPair hapPair() {
        return hapPair;
    }


    /**
     * Returns the starting marker (inclusive).
     * @return the starting marker (inclusive)
     */
    public Marker start() {
        return start;
    }

    /**
     * Returns the ending marker (inclusive).
     * @return the ending marker (inclusive)
     */
    public Marker end() {
        return end;
    }

    /**
     * Returns the IBD segment score.
     * @return the IBD segment score
     */
    public float score() {
        return score;
    }

    /**
     * Returns the starting marker index (inclusive) or -1 if the starting
     * marker index is unknown.
     * @return the starting marker index (inclusive)
     */
    public int startIndex() {
        return startIndex;
    }

    /**
     * Returns the ending marker index (inclusive) or  -1 if the ending
     * marker index is unknown.
     * @return the ending marker index (inclusive)
     */
    public int endIndex() {
        return endIndex;
    }

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append(hapPair.first());
        sb.append(Const.tab);
        sb.append(hapPair.second());
        sb.append(Const.tab);
        sb.append(start.chrom());
        sb.append(Const.tab);
        sb.append(start.pos());
        sb.append(Const.tab);
        sb.append(end.pos());
        sb.append(Const.tab);
        sb.append(df2.format(score));
        return sb.toString();
    }
}
