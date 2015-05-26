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
 * Class {@code IntList} represents a list of integers.
 * Class {@code IntList} supports a {@code clear()} method, but does not
 * support a {@code remove()} method.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class IntList {

    /**
     * The default initial capacity of an {@code IntList}, which is 10.
     */
    public static final int DEFAULT_INIT_CAPACITY = 10;

    private int size;
    private int[] values;

    /**
     * Constructs an {@code IntList} object with the default
     * initial capacity.
     *
     * @see #DEFAULT_INIT_CAPACITY
     */
    public IntList() {
        this(DEFAULT_INIT_CAPACITY);
    }

    /**
     * Constructs an {@code IntList} object with the specified
     * initial capacity.
     *
     * @param initCapacity the initial capacity of this list
     * @throws IllegalArgumentException if {@code initCapacity<0}.
     */
    public IntList(int initCapacity) {
        if (initCapacity < 0) {
            String s = "initCapacity < 0: " + initCapacity;
            throw new IllegalArgumentException(s);
        }
        this.size = 0;
        this.values = new int[initCapacity];
    }

    /**
     * Returns the integer at the specified position in this list.
     * @param index the index of the returned integer.
     * @return the integer at the specified position in this list.
     * @throws IndexOutOfBoundsException if
     * {@code index<0 || index>=size}.
     */
    public int get(int index) {
        if (index < 0 && index >= size) {
            throw new IndexOutOfBoundsException("index=" + index);
        }
        return values[index];
    }

    /**
     * Returns the number of elements in this list.
     * @return the number of elements in this list.
     */
    public int size() {
        return size;
    }

    /**
     * Returns {@code true} if this list has no elements, and returns
     * {@code false} otherwise.
     * @return {@code true} if this list has no elements, and returns
     * {@code false} otherwise.
     */
    public boolean isEmpty() {
        return size==0;
    }

    /**
     * Returns an integer array containing the sequence of elements in this
     * list.
     * @return an integer array containing the sequence of elements in this
     * list.
     */
    public int[] toArray() {
        return Arrays.copyOfRange(values, 0, size);
    }

    /**
     * Adds the specified integer to the end of this list.
     *
     * @param element the integer to be added to the end of this list.
     */
    public void add(int element) {
        if (size==values.length) {
            int newCapacity = (values.length * 3)/2 + 1;
            this.values = Arrays.copyOf(this.values, newCapacity);
        }
        this.values[size++] = element;
    }

    /**
     * Removes all elements from this list.
     */
    public void clear() {
        this.size = 0;
    }

    /**
     * Returns a string representation of this list.  The
     * exact details of the representation are unspecified and
     * subject to change.
     *
     * @return a string representation of this list.
     */
    @Override
    public String toString() {
        return Arrays.toString(Arrays.copyOf(values, size));
    }
}
