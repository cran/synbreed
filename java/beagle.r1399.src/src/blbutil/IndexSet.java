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
 * <p>Class {@code IndexSet} is a set that stores non-negative indices that are
 * less than or equal to a specified maximum value.
 * </p>
 * Class {@code IndexSet} supports a {@code clear()} method, but does not
 * support a {@code remove()} method.
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
     * @param maxIndex the maximum index value..
     */
    public IndexSet(int maxIndex) {
        if (maxIndex < 0) {
            throw new IllegalArgumentException("maxIndex<0: " + maxIndex);
        }
        this.inSet = new boolean[maxIndex+1];
        this.indices = new int[maxIndex+1];
    }

    /**
     * Adds the specified index and value the set.  Returns {@code true}
     * if the set was changed by the operation, and returns {@code false}
     * otherwise.
     *
     * @param index the index to add to this set.
     * @return {@code true} if the set was changed by the operation, and
     * {@code false} otherwise.
     *
     * @throws IndexOutOfBoundsException if
     * {@code index<0 || index>this.maxIndex()}.
     */
    public boolean add(int index) {
        if (inSet[index]==false) {
            indices[size++] = index;
            inSet[index]=true;
            return true;
        }
        else {
            return false;
        }
    }

    /**
     * Returns {@code true} if the set contains the specified index,
     * and returns {@code false} otherwise.
     * @param index an index
     * @return {@code true} if the set contains the specified index,
     * and {@code false} otherwise.
     *
     * @throws IndexOutOfBoundsException if
     * {@code index<0 || index>this.maxIndex()}.
     */
    public boolean contains(int index) {
        return inSet[index];
    }

    /**
     * Returns the size of this set.
     *
     * @return the size of this set.
     */
    public int size() {
        return size;
    }

    /**
     * Returns the maximum permitted index in the set.
     *
     * @return the maximum permitted index in the set.
     */
    public int maxIndex() {
        return indices.length-1;
    }

    /**
     * Removes all indices from the set.
     */
    public void clear() {
        for (int j=0, n=size; j<n; ++j) {
            inSet[indices[j]] = false;
        }
        size = 0;
    }

    /**
     * Returns the specified index in an enumeration of the indices in the set.
     * @param index an index.
     * @return the specified index in an enumeration of the indices in the set.
     * @throws IndexOutOfBoundsException if
     * {@code index<0 || index>=this.size()}
     */
    public int enumIndex(int index) {
        if (index>=size) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        return indices[index];
    }

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     *
     * @return a string representation of {@code this}.
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(80);
        sb.append("[ size=");
        sb.append(size);
        sb.append(" {");
        for (int j=0; j<size; ++j) {
            sb.append(enumIndex(j));
            if (j+1 < size) {
                sb.append(", ");
            }
        }
        sb.append("} ]");
        return sb.toString();
    }
}
