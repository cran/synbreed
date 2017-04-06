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
 * <p>Class {@code ShiftedByteIndexArray} represents an immutable
 * {@code int[]} array that is stored as a {@code byte[]} array.
 * </p>
 * Instances of {@code ShiftedByteIndexArray} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class ShiftedByteIndexArray implements IntArray {

    private static final int shift = 128;
    private final byte[] ba;

    /**
     * Constructs a new {@code ShiftedByteIndexArray} instance.
     * @param ia an array of integers
     * @throws IllegalArgumentException if
     * {@code ia[j] < 0 || ia[j] > 255} for any index {@code j}
     * satisfying  {@code j >= 0 && j < ia.length}
     * @throws NullPointerException if {@code ia == null}
     */
    public ShiftedByteIndexArray(int[] ia) {
        this(ia, 0, ia.length);
    }

    /**
     * Constructs a new {@code ShiftedByteIndexArray} instance from the
     * specified subarray.
     * @param ia an array of integers
     * @param start the first element to be included (inclusive)
     * @param end the last element to be included (exclusive)
     * @throws IllegalArgumentException if
     * {@code ia[j] < 0 || ia[j] > 255} for any index {@code j}
     * satisfying {@code j >= start && j < end}
     * @throws IndexOutOfBoundsException if {@code start < 0 or end > ia.length}
     * @throws IllegalArgumentException if {@code end > start}
     * @throws NullPointerException if {@code ia == null}
     */
    public ShiftedByteIndexArray(int[] ia, int start, int end) {
        if (start > end) {
            throw new IllegalArgumentException("start > end");
        }
        this.ba = new byte[end - start];
        for (int j=start; j<end; ++j) {
            if (ia[j] < 0 || ia[j] > 255) {
                throw new IllegalArgumentException(String.valueOf(ia[j]));
            }
            ba[j - start] = (byte) (ia[j] - shift);
        }
    }

    @Override
    public int size() {
        return ba.length;
    }

    @Override
    public int get(int index) {
        return ba[index] + shift;
    }

    @Override
    public String toString() {
        return this.asString();
    }
}
