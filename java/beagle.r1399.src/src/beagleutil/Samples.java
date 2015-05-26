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
 * <p>Class {@code Samples} stores a list of samples.
 * </p>
 * Class {@code Samples} is immutable.
 *
 * @author Brian L. Browning
 */
public final class Samples {

    private static final SampleIds sampleIds = SampleIds.instance();
    private final int[] idIndexToIndex;
    private final int[] indexToIdIndex;

    /**
     * Constructs a new instance of {@code Samples} corresponding to
     * the specified list of sample identifier indices (see the
     * {@code beagleutil.Ids} class).
     * @param idIndices an array of sample identifier indices.
     *
     * @throws IllegalArgumentException if the specified array
     * has two or more elements that are equal.
     * @throws IndexOutOfBoundsException if any element of the specified
     * array is negative or greater than or equal to
     * {@code beagleutil.Ids.size()}.
     * @throws NullPointerException if {@code indices==null}.
     */
    public Samples(int[] idIndices) {
        int[] copy = idIndices.clone();
        this.idIndexToIndex = idIndexToIndex(copy);
        this.indexToIdIndex = copy;
    }

    private static int[] idIndexToIndex(int[] idIndices) {
        int[] idIndexToIndex = new int[sampleIds.size()];
        Arrays.fill(idIndexToIndex, -1);
        for (int j=0; j<idIndices.length; ++j) {
            int idIndex = idIndices[j];
            if (idIndexToIndex[idIndex] != -1) {
                String s = "duplicate sample: " + sampleIds.id(idIndex) +
                        " (ID index: " + idIndex + ")";
                throw new IllegalArgumentException(s);
            }
            else {
                idIndexToIndex[idIndex] = j;
            }
        }
        return idIndexToIndex;
    }

    /**
     * Constructs and returns a }Samples} instance
     * corresponding to the specified list of sample identifiers.
     * @param ids an array of sample identifiers.
     * @return a }Samples} instance
     * corresponding to the specified list of sample identifiers.
     *
     * @throws IllegalArgumentException if the specified array
     * has two or more elements that are equal as strings.
     * @throws NullPointerException if {@code ids==null}.
     */
    public static Samples fromIds(String[] ids) {
        int[] indices = new int[ids.length];
        for (int j=0; j<ids.length; ++j) {
            indices[j] = sampleIds.indexOf(ids[j]);
        }
        return new Samples(indices);
    }

    /**
     * Returns a hash code value for the object.
     * @return a hash code value for the object.
     */
    @Override
    public int hashCode() {
        return Arrays.hashCode(this.indexToIdIndex);
    }

    /**
     * Returns {@code true} if the specified object is a
     * {@code Samples} object which represents the same ordered
     * list of samples as {@code this}, and returns {@code false}
     * otherwise.
     * @param obj the object to be tested for equality with {@code this}.
     * @return {@code true} if the specified object is a
     * {@code Samples} object which represents the same ordered
     * list of samples as {@code this}, and returns {@code false}
     * otherwise.
     */
    @Override
    public boolean equals(Object obj) {
        if (this==obj) {
            return true;
        }
        if (obj == null || this.getClass() != obj.getClass()) {
            return false;
        }
        final Samples other = (Samples) obj;
        return Arrays.equals(this.indexToIdIndex, other.indexToIdIndex);
    }

    /**
     * Returns the sample identifier index corresponding to the specified
     * index in this list of samples.
     * @param index a sample index.
     * @return the sample identifier index corresponding to the specified
     * index in this list of samples.
     * @throws IndexOutOfBoundsException if
     * {@code index<0 || index>=this.nSamples()}.
     */
    public int idIndex(int index) {
        return indexToIdIndex[index];
    }

    /**
     * Returns the sample index that corresponds to the
     * specified sample identifier index, or returns {@code -1}
     * if there is no corresponding sample in this list of samples.
     * @param idIndex a sample identifier index.
     * @return the sample index that corresponds to the
     * specified sample identifier index, or returns {@code -1}
     * if there is no corresponding sample in this list of samples.
     * @throws IndexOutOfBoundsException if {@code index<0}.
     */
    public int index(int idIndex) {
        if (idIndex >= idIndexToIndex.length) {
            return -1;
        }
        return idIndexToIndex[idIndex];
    }

    /**
     * Returns the sample index that corresponds to the
     * specified sample identifier, or returns {@code -1}
     * if there is no corresponding sample in this list of samples.
     * @param id a sample identifier.
     * @return the sample index that corresponds to the
     * specified sample identifier index, or returns {@code -1}
     * if there is no corresponding sample in this list of samples.
     * @throws NullPointerException if {@code id==null}.
     */
    public int index(String id) {
        int idIndex = SampleIds.instance().indexIfIndexed(id);
        if (idIndex != -1) {
            return index(idIndex);
        }
        else {
            return -1;
        }
    }

    /**
     * Returns the number of samples in this list.
     * @return the number of samples in this list.
     */
    public int nSamples() {
        return indexToIdIndex.length;
    }

    /**
     * Returns the sample identifier for the specified sample in this list
     * of samples.
     * @param index a sample index
     * @return the sample identifier for the specified sample in this list
     * of samples.
     * @throws IndexOutOfBoundsException if
     * {@code index<0 || index>=this.nSamples()}.
     */
    public String id(int index) {
        return sampleIds.id(indexToIdIndex[index]);
    }

    /**
     * Returns this list of samples as an array of sample identifiers.
     * The returned array has length {@code this.nSamples()}, and it
     * satisfies {@code this.ids()[index].equals(this.id(index))} for
     * {@code 0<=index<this.nSamples()}.
     * @return this list of samples as an array of sample identifiers.
     */
    public String[] ids() {
        String[] ids = new String[indexToIdIndex.length];
        for (int j=0; j<ids.length; ++j) {
            ids[j] = sampleIds.id(indexToIdIndex[j]);
        }
        return ids;
    }

    /**
     * Returns a string representation of {@code this}.
     * The exact details of the representation are unspecified and
     * subject to change.
     * @return a string representation of {@code this}.
     */
    @Override
    public String toString() {
        return Arrays.toString(indexToIdIndex);
    }
}
