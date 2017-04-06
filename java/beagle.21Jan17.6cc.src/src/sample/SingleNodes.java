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
package sample;

/**
 * <p>Class {@code SingleNodes} stores ordered node pairs and associated values.
 * </p>
 * <p>Instances of class {@code SingleNodes} are not thread safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class SingleNodes {

    private static final float loadFactor = 0.75f;

    private int size;
    private int capacity; // required to be a power of 2
    private int rehashThreshold;

    private int[] index;
    private int[] node1;
    private int[] node2;
    private float[] value;

    /**
     * Creates a new instance of {@code SingleNodes} that has an
     * initial value of 0 for each ordered node pair. The first node
     * has index 0.
     */
    public SingleNodes() {
        this.size = 0;
        this.capacity = (1<<10);
        this.rehashThreshold = (int) (loadFactor * capacity);
        this.index = new int[capacity];
        this.node1 = new int[capacity];
        this.node2 = new int[capacity];
        this.value = new float[capacity];
    }

    private static long hash1(int node1, int node2) {
        long hash = 5;
        hash = 71 * hash + node1;
        hash = 71 * hash + node2;
        return hash;
    }

    private static long hash2(int node1, int node2) {
        long hash = 7;
        hash = 97 * hash + node1;
        hash = 97 * hash + node2;
        return hash;
    }

    /*
     * Return the storage index for specified node pair.  If the key is not
     * currently stored in the hash table, the index at which the value
     * should be stored is returned.
     */
    private int index(int node1, int node2) {
        long h1 = hash1(node1, node2);
        long h2 = hash2(node1, node2);
        if ((h2 & 1)==0) {
            // h2 must be relatively prime to maxSize, which is a power of 2
            ++h2;
        }
        long l = h1;
        for (int k=0; k<capacity; ++k) {
            int i = (int) (l % capacity);
            if (value[i]==0.0
                    || (this.node1[i]==node1 && this.node2[i]==node2)) {
                return i;
            }
            l += h2;
        }
        assert false;
        return -1;
    }

    /*
     * Increases the capacity of the internal hash table.
     */
    private void rehash() {
        assert this.size>=this.rehashThreshold;
        int newMaxSize = 2*capacity;
        if (newMaxSize<0) {
            throw new IllegalStateException("hash table overflow");
        }
        int[] oldIndex = index;
        int[] oldNode1 = node1;
        int[] oldNode2 = node2;
        float[] oldValue = value;

        capacity = newMaxSize;
        index = new int[newMaxSize];
        node1 = new int[newMaxSize];
        node2 = new int[newMaxSize];
        value = new float[newMaxSize];

        for (int j=0; j<size; ++j) {
            int oldInd = oldIndex[j];
            int newIndex = index(oldNode1[oldInd], oldNode2[oldInd]);
            index[j] = newIndex;
            node1[newIndex] = oldNode1[oldInd];
            node2[newIndex] = oldNode2[oldInd];
            value[newIndex] = oldValue[oldInd];
        }
        rehashThreshold = (int) (loadFactor * capacity);
    }

    /**
     * Adds the specified positive value to the stored value of the specified
     * node pair.
     *
     * @param node1 the first node
     * @param node2 the second node
     * @param value the value
     *
     * @throws IllegalArgumentException if {@code node1 < 0 || node2 < 0}
     * @throws IllegalArgumentException if
     * {@code value <= 0 || (Double.isFinite(value) == false)}
     */
    public void sumUpdate(int node1, int node2, float value) {
        if (node1 < 0) {
            throw new IllegalArgumentException(String.valueOf(node1));
        }
        if (node2 < 0) {
            throw new IllegalArgumentException(String.valueOf(node2));
        }
        if (value <= 0 || (Double.isFinite(value)==false) ) {
            throw new IllegalArgumentException(String.valueOf(value));
        }
        int i = index(node1, node2);
        boolean addNode = (this.value[i]==0f);
        this.value[i] += value;
        if (addNode) {
            this.index[size++] = i;
            this.node1[i] = node1;
            this.node2[i] = node2;
            if (this.size>=this.rehashThreshold) {
                rehash();
            }
        }
    }

    /**
     * Returns the number of node pairs with non-zero value.
     * @return the number of node pairs with non-zero value
     */
    public int size() {
        return size;
    }

    private void checkSize(int index) {
        if (index>=size()) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
    }

    /**
     * Returns the first node of the specified node pair in the list of
     * node pairs with non-zero value.  Repeated invocations of this
     * method with the same parameter will return the same value if
     * node values are not modified between invocations. If
     * {@code (index >= 0 && index < this.size())}, then the following
     * expression will always evaluate to {@code true}:<br>
     * {@code (this.value(this.enumNode1(index),
     * this.enumNode2(index)) == this.enumValue(index))}.
     *
     * @param index an index in a list of node pairs with non-zero
     * value
     * @return the first node of the specified node pair in a list of
     * node pairs with non-zero value
     *
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public int enumNode1(int index) {
        checkSize(index);
        return node1[this.index[index]];
    }

    /**
     * Returns the second node of the specified node pair in a list of
     * node pairs with non-zero value.  Repeated invocations of this
     * method with the same parameter will return the same value if
     * node values are not modified between invocations. If
     * {@code (index >= 0 && index < this.size())}, then the following
     * expression will always evaluate to {@code true}:<br>
     * {@code (this.value(this.enumNode1(index),
     * this.enumNode2(index)) == this.enumValue(index))}.
     *
     * @param index an index in a list of node pairs with non-zero value
     * @return the second node of the specified node pair in a list of
     * node pairs with non-zero value
     *
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public int enumNode2(int index) {
        checkSize(index);
        return node2[this.index[index]];
    }

    /**
     * Returns the value of the specified node pair in a list of
     * node pairs with non-zero value.  Repeated invocations of this
     * method with the same parameter will return the same value if
     * node values are not modified between invocations. If
     * {@code (index >= 0 && index < this.size())}, then the following
     * expression will always evaluate to {@code true}:<br>
     * {@code (this.value(this.enumNode1(index),
     * this.enumNode2(index)) == this.enumValue(index))}.
     *
     * @param index an index in a list of node pairs with non-zero value
     * @return the value of the specified ordered node pair in a list of
     * node pairs with non-zero value
     *
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public float enumValue(int index) {
        checkSize(index);
        return value[this.index[index]];
    }

    /**
     * Returns the value of the specified node pair.
     *
     * @param node1 the first node
     * @param node2 the second node
     * @return the value of the specified node pair
     * @throws IllegalArgumentException if {@code node1 < 0 || node2 < 0}
     */
    public float value(int node1, int node2) {
        if (node1 < 0) {
            throw new IllegalArgumentException(String.valueOf(node1));
        }
        if (node2 < 0) {
            throw new IllegalArgumentException(String.valueOf(node2));
        }
        return value[index(node1, node2)];
    }

    /**
     * Sets the value of each ordered node pair to 0.
     */
    public void clear() {
        for (int j=0; j<this.size; ++j) {
            value[index[j]] = 0f;
        }
        size = 0;
    }

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     *
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(80);
        sb.append("size=");
        sb.append(size);
        for (int j=0; j<size; ++j) {
            sb.append(" (");
            sb.append(j);
            sb.append(": node1=");
            sb.append(enumNode1(j));
            sb.append(" node2=");
            sb.append(enumNode2(j));
            sb.append(" value=");
            sb.append(enumValue(j));
            sb.append(") ");
        }
        return sb.toString();
    }
}
