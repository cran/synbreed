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
 * <p>Class {@code IndexMap} is a map whose keys are a bounded set of
 * non-negative integers and whose values are integers.
 * </p>
 * <p>Class {@code IndexMap} supports a {@code clear()} method, but it does not
 * support a {@code remove()} method.
 * </p>
 * <p>Class {@code IndexMap} is not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class IndexMap {

    private final int nil;
    private final int[] values;
    private final int[] keys;
    private int size = 0;

    /**
     * Creates a new instance of {@code IndexMap} whose {@code nil()} method
     * will return the specified {@code nil} value.
     * @param maxKey the maximum key
     * @param nil the value that will be returned by the instance's
     * {@code get()} method if a key has no assigned value
     * @throws IllegalArgumentException if {@code maxKey < 0}
     */
    public IndexMap(int maxKey, int nil) {
        if (maxKey < 0) {
            throw new IllegalArgumentException(String.valueOf(maxKey));
        }
        this.nil = nil;
        this.values = new int[maxKey+1];
        this.keys = new int[maxKey+1];
        Arrays.fill(values, nil);
    }

    /**
     * Returns the value that is returned by {@code this.get()} if
     * a key has no assigned value.
     * @return the value that is returned by {@code this.get()} if
     * a key has no assigned value
     */
    public int nil() {
        return nil;
    }

    /**
     * Adds the specified key and value to the map. If the map
     * contains a value for the specified key when the method is invoked,
     * the old value is replaced by the specified value.
     *
     * @param key the key
     * @param value the value
     * @return the previous value associated with {@code key}, or
     * {@code this.nil()} if no such previous value exists
     *
     * @throws IllegalArgumentException if {@code value == this.nil()}
     * @throws IndexOutOfBoundsException if
     * {@code key < 0 || key > this.maxKey()}
     */
    public int put(int key, int value) {
        if (value==nil) {
            throw new IllegalArgumentException("value==nil()");
        }
        int prevValue = values[key];
        if (prevValue == nil) {
            keys[size++] = key;
        }
        this.values[key] = value;
        return prevValue;
    }

    /**
     * Returns the value associated with the specified key, or
     * {@code this.nil()} if the specified key is not contained in the map.
     * @param key the key
     * @return the value associated with the specified key, or
     * {@code this.nil()} if the specified key is not contained in the map.
     *
     * @throws IndexOutOfBoundsException if
     * {@code key < 0 || key > this.maxKey()}
     */
    public int get(int key) {
        return values[key];
    }

    /**
     * Returns the number of key-value pairs in the map.
     *
     * @return the number of key-value pairs in the map
     */
    public int size() {
        return size;
    }

    /**
     * Returns the maximum key.
     *
     * @return the maximum key
     */
    public int maxKey() {
        return keys.length-1;
    }

    /**
     * Removes all key-value pairs from the map.
     */
    public void clear() {
        for (int j=0, n=size; j<n; ++j) {
            values[keys[j]] = nil;
        }
        size = 0;
    }

    /**
     * Returns the specified key in an enumeration of the keys in the map.
     * @param index an index of an element in the enumeration
     * @return the specified key in an enumeration of the keys-value
     * pairs in the map
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public int enumeratedKey(int index) {
        if (index>=size) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        return keys[index];
    }

    /**
     * Returns the value associated with the specified key
     * in an enumeration of the keys in the map.
     * If {@code (index >= 0 && index < this.size())}, then the returned value
     * will satisfy:
     * {@code this.get(this.enumeratedKey(index)==this.enumeratedValue(index)}.
     * @param index an index of an element in the enumeration
     * @return the value associated with the specified key
     * in an enumeration of the keys in the map
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public int enumeratedValue(int index) {
        if (index>=size) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        return values[keys[index]];
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
        sb.append("size=");
        sb.append(size);
        sb.append(" {");
        for (int j=0; j<size; ++j) {
            sb.append(enumeratedKey(j));
            sb.append(" : ");
            sb.append(enumeratedValue(j));
            if (j+1 < size) {
                sb.append(Const.comma);
            }
        }
        sb.append("}");
        return sb.toString();
    }
}
