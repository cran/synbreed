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
package beagleutil;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * <p>Class {@code ThreadSafeIndexer} indexes objects.
 * </p>
 * Instances of class {@code ThreadSafeIndexer} are thread-safe.
 *
 * @param <T> the type parameter.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class ThreadSafeIndexer<T> {

    /**
     * The default initial capacity, which is 500.
     */
    public static final int DEFAULT_INIT_CAPACITY = 500;

    private final List<T> list ;
    private final Map<T, Integer> map;

    /**
     * Creates a new {@code ThreadSafeIndexer} instance with the default
     * initial capacity.
     *
     * @see #DEFAULT_INIT_CAPACITY
     */
    public ThreadSafeIndexer() {
        this(DEFAULT_INIT_CAPACITY);
    }

    /**
     * Creates a new {@code ThreadSafeIndexer}instance with the specified
     * initial capacity.
     * @param initCapacity the initial capacity
     * @throws IllegalArgumentException if {@code initCapacity < 1}
     */
    public ThreadSafeIndexer(int initCapacity) {
        if (initCapacity < 1) {
            throw new IllegalArgumentException(String.valueOf(initCapacity));
        }
        this.list = new ArrayList<>(initCapacity);
        this.map = new HashMap<>(initCapacity);
    }

    /**
     * Returns the index of the specified object.  If the object
     * is not yet indexed, the object will be indexed. Indices
     * are assigned in consecutive order beginning with 0.
     * @param object the object whose index will be retrieved
     * @return the index of the specified object
     * @throws NullPointerException if {@code object==null}
     */
    public synchronized int getIndex(T object) {
        if (object==null) {
            throw new NullPointerException();
        }
        if (map.keySet().contains(object)) {
            return map.get(object);
        }
        else {
            int idIndex = list.size();
            list.add(object);
            map.put(object, idIndex);
            return idIndex;
        }
    }

    /**
     * Returns the index of the specified object, or returns
     * {@code -1} if the specified object is not indexed.
     *
     * @param object an object
     * @return the index of the specified object, or
     * {@code -1} if the specified object is not indexed
     *
     * @throws NullPointerException if {@code object == null}.
     */
    public synchronized int getIndexIfIndexed(T object) {
        if (object==null) {
            throw new NullPointerException();
        }
        if (map.keySet().contains(object)) {
            return map.get(object);
        }
        else {
            return -1;
        }
    }

    /**
     * Returns the number of indexed objects.
     * @return the number of indexed objects
     */
    public synchronized int size() {
        return list.size();
    }

    /**
     * Returns the object with the specified index.
     * @param index an object index
     * @return the object with the specified index
     * @throws IndexOutOfBoundsException if
     * {@code  index<0 || index>=this.size()}
     */
    public synchronized T item(int index) {
        return list.get(index);
    }

    /**
     * Returns an listed of indexed objects. The returned list will
     * have size {@code this.size()}, and it will satisfy
     * {@code this.items().get(k).equals(this.item(k))==true}
     * for {@code  0 <= k && k < this.size()}
     *
     * @return an array of objects
     */
    public synchronized List<T> items() {
        return new ArrayList<>(list);
    }

    /**
     * Returns {@code this.items().toString()}.
     * @return a string representation of {@code this}
     */
    @Override
    public synchronized String toString() {
        return this.items().toString();
    }
}
