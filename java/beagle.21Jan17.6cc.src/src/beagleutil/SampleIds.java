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

import java.util.Arrays;

/**
 * <p>Class {@code SampleIds} is a singleton class that represents a
 * list of sample identifiers.
 * </p>
 * The singleton instance of {@code SampleIds} is thread-safe.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class SampleIds {

    private static final SampleIds sampleIds = new SampleIds();

    private final ThreadSafeIndexer<String> instance;

    private SampleIds() {
        // private constructor to restrict instantiation.
        int initCapacity = 5000;
        this.instance = new ThreadSafeIndexer<>(initCapacity);
    }

    /**
     * Returns the singleton {@code SampleIds} instance.
     * @return the singleton {@code SampleIds} instance
     */
    public static SampleIds instance() {
        return sampleIds;
    }

    /**
     * Returns the index of the specified sample identifier.  If
     * the sample identifier is not yet indexed, the sample identifier
     * will be indexed.  Sample identifier indices are assigned in
     * consecutive order beginning with 0.
     * @param id a sample identifier
     * @return the index of the specified sample identifier
     * @throws IllegalArgumentException if {@code id.isEmpty()}
     * @throws NullPointerException if {@code id == null}
     */
    public int getIndex(String id) {
        if (id.isEmpty()) {
            throw new IllegalArgumentException("id.isEmpty()");
        }
        return instance.getIndex(id);
    }

    /**
     * Returns the index of the specified sampled identifier, or returns
     * {@code -1} if the specified sample identifier is not indexed.
     *
     * @param id a sample identifiers
     * @return the index of the specified sampled identifier, or
     * {@code -1} if the specified sample identifier is not indexed
     *
     * @throws IllegalArgumentException if {@code id.isEmpty()}
     * @throws NullPointerException if {@code id == null}
     */
    public int getIndexIfIndexed(String id) {
        if (id.isEmpty()) {
            throw new IllegalArgumentException("id.isEmpty()");
        }
        return instance.getIndexIfIndexed(id);
    }

    /**
     * Returns the number of indexed sample identifiers.
     * @return the number of indexed samples identifiers
     */
    public int size() {
        return instance.size();
    }

    /**
     * Returns the sample identifier with the specified index.
     * @param index a sample identifier index
     * @return the specified sample identifier
     * @throws IndexOutOfBoundsException if
     * {@code  index < 0 || index >= this.size()}
     */
    public String id(int index) {
        return instance.item(index);
    }

    /**
     * Returns the list of indexed sample identifiers as an array.
     * The returned array will have length {@code this.size()}, and
     * it will satisfy
     * {@code this.ids()[k].equals(this.id(k)) == true}
     * for {@code  0 <= k && k < this.size()}.
     *
     * @return an array of sample identifiers
     */
    public String[] ids() {
        return instance.items().toArray(new String[0]);
    }

    /**
     * Returns {@code java.util.Arrays.toString(this.ids())}.
     *
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        return Arrays.toString(this.ids());
    }
}
