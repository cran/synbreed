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
 * <p>Class {@code IndexSet} is a set that stores non-negative indices that are
 * less than or equal to a specified maximum value.
 * </p>
 * <p>Class {@code IndexSet} supports a {@code clear()} method, but it does not
 * support a {@code remove()} method.
 * </p>
 * <p>Class {@code IndexSet} is not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class IndexSet {

    private final boolean[] inSet;
    private final int[] indices;
    private int size = 0;

    /**
     * Creates a new instance of {@code IndexSet} that can contain
     * non-negative integer indices that are less than or equal to the specified
     * maximum value.
     *
     * @param max the maximum element that is permitted in the set.
     * @throws IllegalArgumentException if {@code max < 0}
     */
    public IndexSet(int max) {
        if (max < 0) {
            throw new IllegalArgumentException(String.valueOf(max));
        }
        this.inSet = new boolean[max+1];
        this.indices = new int[max+1];
    }

    /**
     * Adds the specified element to the set.
     *
     * @param element an element to add to this set.
     * @return {@code true} if the set was changed by the operation, and
     * {@code false} otherwise.
     *
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index > this.maxPermittedIndex()}
     */
    public boolean add(int element) {
        if (inSet[element]==false) {
            indices[size++] = element;
            inSet[element]=true;
            return true;
        }
        else {
            return false;
        }
    }

    /**
     * Returns {@code true} if the set contains the specified element,
     * and returns {@code false} otherwise.
     * @param element an element
     * @return {@code true} if the set contains the specified element
     *
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index > this.maxPermittedIndex()}
     */
    public boolean contains(int element) {
        return inSet[element];
    }

    /**
     * Returns the number of elements in this set.
     *
     * @return the number of elements in this set
     */
    public int size() {
        return size;
    }

    /**
     * Returns the maximum permitted element in the set.
     *
     * @return the maximum permitted element in the set
     */
    public int maxPermittedElement() {
        return indices.length-1;
    }

    /**
     * Removes all elements from the set.
     */
    public void clear() {
        for (int j=0, n=size; j<n; ++j) {
            inSet[indices[j]] = false;
        }
        size = 0;
    }

    /**
     * Returns the specified element in an enumeration of the elements in the
     * set.
     * @param enumIndex an index of an element in the enumeration
     * @return the specified element in an enumeration of the elements in the
     * set
     * @throws IndexOutOfBoundsException if
     * {@code enumIndex < 0 || enumIndex >= this.size()}
     */
    public int enumeratedValue(int enumIndex) {
        if (enumIndex>=size) {
            throw new IndexOutOfBoundsException(String.valueOf(enumIndex));
        }
        return indices[enumIndex];
    }

    /**
     * Returns an array containing the elements in this set.
     * @return an array containing the elements in this set
     */
    public int[] toArray() {
        return Arrays.copyOf(indices, size);
    }

    /**
     * Returns {@code java.util.Arrays.toString(this.toArray())}.
     *
     * @return {@code java.util.Arrays.toString(this.toArray())}
     */
    @Override
    public String toString() {
        return Arrays.toString(toArray());
    }
}
