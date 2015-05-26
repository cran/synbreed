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
package blbutil;

/**
 * <p>Class {@code  IntPair} represents an ordered pair of integers.
 * </p>
 * Instances of class {@code  IntPair} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class IntPair implements Comparable<IntPair> {

    private final int first;
    private final int second;

    /**
     * Constructs an {@code  IntPair} instance that represents the
     * specified ordered pair of integers.
     * @param first the first element of the ordered pair of integers.
     * @param second the second element of the ordered pair of integers.
     */
    public IntPair(int first, int second) {
        this.first = first;
        this.second = second;
    }

    /**
     * Returns the first integer in the ordered pair of integers.
     * @return the first integer in the ordered pair of integers.
     */
    public int first() {
        return first;
    }

    /**
     * Returns the second integer in the ordered pair of integers.
     * @return the second integer in the ordered pair of integers.
     */
    public int second() {
        return second;
    }

    /**
     * Compares the specified object with this {@code  IntPair} for
     * equality.  Returns {@code  true} if the specified object
     * is an {@code  IntPair} that represents the same ordered
     * pair of integers as {@code  this}, and returns {@code  false}
     * otherwise.
     * @param obj the object to be compared for equality with this
     * {@code  IntPair}.
     * @return {@code  true} if the specified object is an {@code  IntPair}
     * that represents the same ordered pair of integers as {@code  this},
     * and returns {@code  false} otherwise.
     */
    @Override
    public boolean equals(Object obj) {
        if (obj == this) {
            return true;
        }
        if (!(obj instanceof IntPair)) {
            return false;
        }
        IntPair other = (IntPair) obj;
        return (this.first==other.first) && (this.second==other.second);
    }

     /**
     * Returns a hash code value for the object.
     *
     * <p>The hash code is defined by the following calculation:
     * </p>
     * <pre>
        int hash = 5;
        hash = 29 * hash + this.first;
        hash = 29 * hash + this.second;
     * </pre>
     * @return a hash code value for the object.
     */
    @Override
    public int hashCode() {
        int hash = 5;
        hash = 29 * hash + this.first;
        hash = 29 * hash + this.second;
        return hash;
    }

    /**
     * Returns a string representation of {@code this}.  The string
     * representation is {@code  "[i1, i2]"} where
     * {@code i1} and {@code i2} are the first and second integers
     * in the ordered pair of integers represented by {@code this}.
     * @return a string representation of {@code this}.
     */
    @Override
    public String toString() {
        return "[" + first + ", " + second + "]";
    }

    /**
     * Returns -1, 0, or 1 depending on whether {@code this} is
     * less than, equal, or greater than the specified {@code IntPair}
     * object. {@code IntPair} instances are ordered using
     * lexicographical order.
     * @param other an {@code IntPair} instance to be compared to
     * {@code this}.
     * @return -1, 0, or 1 depending on whether {@code this} is
     * less than, equal, or greater than the specified {@code IntPair}
     * object.
     */
    @Override
    public int compareTo(IntPair other) {
        if (this.first < other.first) {
            return -1;
        }
        if (this.first > other.first) {
            return 1;
        }
        if (this.second < other.second) {
            return -1;
        }
        if (this.second > other.second) {
            return 1;
        }
        return 0;
    }
}
