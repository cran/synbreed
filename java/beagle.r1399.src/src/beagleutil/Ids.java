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
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * <p>Class {@code Ids} stores a list of identifiers.
 * </p>
 * Class {@code Ids} is thread-safe.
 *
 * @author Brian L. Browning
 */
public final class Ids {

    /**
     * The default initial capacity, which is 5000.
     */
    public static final int DEFAULT_INIT_CAPACITY = 5000;
    private final List<String> idList ;
    private final Map<String, Integer> idMap;

    /**
     * Creates an {@code Ids} instance with the default initial
     * capacity.
     *
     * @see #DEFAULT_INIT_CAPACITY
     */
    public Ids() {
        this(DEFAULT_INIT_CAPACITY);
    }

    /**
     * Creates an {@code Ids} instance with the specified initial
     * capacity.
     * @param initCapacity the initial capacity.
     * @throws IllegalArgumentException if {@code initCapacity < 1}.
     */
    public Ids(int initCapacity) {
        if (initCapacity < 1) {
            throw new IllegalArgumentException("initCapacity: " + initCapacity);
        }
        this.idList = new ArrayList<>(initCapacity);
        this.idMap = new HashMap<>(initCapacity);
    }

    /**
     * Returns the index of the specified identifier.  If
     * {@code this.indexIfIndexed(id)==-1}, an index will be
     * assigned to the specified identifier.  Identifier indices
     * are assigned in consecutive order beginning with 0.
     * @param id an identifier.
     * @return the index of the specified identifier.
     * @throws IllegalArgumentException if {@code id.isEmpty()}.
     * @throws NullPointerException if {@code id==null}.
     */
    public synchronized int indexOf(String id) {
        if (id.isEmpty()) {
            throw new IllegalArgumentException("id.isEmpty()==true");
        }
        if (idMap.keySet().contains(id)) {
            return idMap.get(id);
        }
        else {
            int idIndex = idList.size();
            idList.add(id);
            idMap.put(id, idIndex);
            return idIndex;
        }
    }

    /**
     * Returns the index of the specified identifier, or returns
     * {@code -1} if the specified identifier is not indexed.
     *
     * @param id an identifiers.
     * @return the index of the specified identifier, or
     * {@code -1} if the specified identifier is not indexed.
     *
     * @throws IllegalArgumentException if {@code id.isEmpty()}.
     * @throws NullPointerException if {@code id==null}.
     */
    public synchronized int indexIfIndexed(String id) {
        if (id.isEmpty()) {
            throw new IllegalArgumentException("id.isEmpty()==true");
        }
        if (idMap.keySet().contains(id)) {
            return idMap.get(id);
        }
        else {
            return -1;
        }
    }

    /**
     * Returns the number of identifiers.
     * @return the number of identifiers.
     */
    public synchronized int size() {
        return idList.size();
    }

    /**
     * Returns the specified identifier.
     * @param index a identifier index.
     * @return the specified identifier.
     * @throws IndexOutOfBoundsException if
     * {@code  index<0 || index>=this.size()}.
     */
    public synchronized String id(int index) {
        return idList.get(index);
    }

    /**
     * Returns an array of identifiers. The returned array will
     * have length {@code this.size()}, and it will satisfy
     * {@code this.ids()[k].equals(this.id(k))==true}
     * for {@code  0<=k<this.size()}.
     *
     * @return an array of identifiers.
     */
    public synchronized String[] ids() {
        return idList.toArray(new String[0]);
    }

    /**
     * Returns a string representation of {@code this}.
     * The returned string is equal to
     * {@code java.util.Arrays.toString(this.ids())}.
     * @return a string representation of {@code this}.
     */
    @Override
    public String toString() {
        return Arrays.toString(this.ids());
    }
}
