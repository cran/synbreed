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

import java.util.Arrays;

/**
 * <p>Interface {@code IntArray} represents an immutable {@code int[]} array.
 * </p>
 * Instances of class {@code IntArray} are required to be immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface IntArray {

    /**
     * Returns the number of elements.
     * @return the number of elements.
     */
    int size();

    /**
     * Returns the specified array element.
     * @param index an array index
     * @return the specified array element
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 ||  index >= this.size()}
     */
    int get(int index);

    /**
     * Returns a string representation of this {@code IntArray} by applying
     * {@code java.utils.Arrays.toString()} to an equivalent {@code int[]}
     * object.
     *
     * @return a string representation of this {@code IntArray}
     */
    @Override
    String toString();

    /**
     * Returns a string representation of this {@code IntArray} by applying
     * {@code java.utils.Arrays.toString()} to an equivalent {@code int[]}
     * object.
     *
     * @return a string representation of this {@code IntArray}.
     */
    default String asString() {
        int[] ia = new int[size()];
        for (int j=0; j<ia.length; ++j) {
            ia[j] = get(j);
        }
        return Arrays.toString(ia);
    }

    /**
     * Returns a new {@code IntArray} instance that has the same
     * sequence of nonnegative integers as the specified array.
     * @param ia the array of non-negative integers to be copied
     * @param min the minimum element in the specified array
     * @param max the maximum element in the specified array
     * @return a new {@code IntArray} instance that has
     * the same sequence of integers as the specified array
     * @throws IllegalArgumentException if {@code minElement > maxElement}
     * @throws IllegalArgumentException if an out-of-range
     * element is detected
     * @throws NullPointerException if {@code ia == null}
     */
    static IntArray create(int[] ia, int min, int max) {
        if (min > max) {
            throw new IllegalArgumentException("min > max");
        }
        if (min >= 0) {
            if (max < 128) {
                return new ByteIndexArray(ia);
            }
            if (max < 256) {
                return new ShiftedByteIndexArray(ia);
            }
            else if (max < 65536) {
                return new CharIndexArray(ia);
            }
        }
        return new WrappedIntArray(ia);
    }
}
