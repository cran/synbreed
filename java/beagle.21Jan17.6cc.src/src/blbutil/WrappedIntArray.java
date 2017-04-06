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
 * <p>Class {@code WrappedIntArray} represents an immutable
 * {@code int[]} array.
 * </p>
 * Instances of {@code WrappedIntArray} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class WrappedIntArray implements IntArray {

    private final int[] ia;

    /**
     * Constructs a new {@code CharCompressedIntArray} instance.
     * @param ia an array of integers
     * @throws NullPointerException if {@code ia == null}
     */
    public WrappedIntArray(int[] ia) {
        this.ia = ia.clone();
    }

    @Override
    public int size() {
        return ia.length;
    }

    @Override
    public int get(int index) {
        return ia[index];
    }

    @Override
    public String toString() {
        return Arrays.toString(ia);
    }
}
