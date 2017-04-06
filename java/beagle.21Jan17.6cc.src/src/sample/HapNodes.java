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
 * <p>Class {@code HapNodes} stores nodes and associated values.
 * </p>
 * <p>Instances of class {@code HapNodes} are not thread safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class HapNodes {

    private static final float loadFactor = 0.75f;

    private int size;
    private int capacity; // required to be a power of 2
    private int rehashThreshold;

    private int[] index;
    private int[] node;
    private float[] value;

    /**
     * Creates a new instance of {@code HapNodes} that has an
     * initial value of 0 for each node.
     */
    public HapNodes() {
        this.size = 0;
        this.capacity = (1<<10);
        this.rehashThreshold = (int) (loadFactor * capacity);
        this.index = new int[capacity];
        this.node = new int[capacity];
        this.value = new float[capacity];
    }

    /*
     * Return the storage index for specified node.  If the key is not
     * currently stored in the hash table, the index at which the value
     * should be stored is returned.
     */
    private int index(int node) {
        long l = (71 * 5) + node;
        long h2 = (97 * 7) + node;
        if ((h2 & 1)==0) {
            // h2 must be relatively prime to maxSize, which is a power of 2
            ++h2;
        }
        for (int k=0; k<capacity; ++k) {
            int i = (int) (l % capacity);
            if (value[i]==0.0 || (this.node[i]==node)) {
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
        int[] oldNode = node;
        float[] oldValue = value;

        capacity = newMaxSize;
        index = new int[newMaxSize];
        node = new int[newMaxSize];
        value = new float[newMaxSize];

        for (int j=0; j<size; ++j) {
            int oldInd = oldIndex[j];
            int newIndex = index(oldNode[oldInd]);
            index[j] = newIndex;
            node[newIndex] = oldNode[oldInd];
            value[newIndex] = oldValue[oldInd];
        }
        rehashThreshold = (int) (loadFactor * capacity);
    }

    /**
     * Adds the specified value to the stored value of the specified
     * node.
     *
     * @param node the node
     * @param value the value
     *
     * @throws IllegalArgumentException if {@code node < 0}
     * @throws IllegalArgumentException if
     * {@code value <= 0 || (Double.isFinite(value) == false)}
     */
    public void sumUpdate(int node, float value) {
        if (node < 0) {
            throw new IllegalArgumentException(String.valueOf(node));
        }
        if (value <= 0 || (Double.isFinite(value)==false) ) {
            throw new IllegalArgumentException(String.valueOf(value));
        }
        int i = index(node);
        boolean addNode = (this.value[i]==0f);
        this.value[i] += value;
        if (addNode) {
            this.index[size++] = i;
            this.node[i] = node;
            if (this.size>=this.rehashThreshold) {
                rehash();
            }
        }
    }

    /**
     * Returns the number of nodes with non-zero value.
     * @return the number of nodes with non-zero value
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
     * Returns the specified node in a list of nodes with non-zero value.
     * Repeated invocations of this method with the same parameter will
     * return the same value if node values are not modified between
     * invocations. If {@code (index >= 0 && index < this.size())}, then the
     * following expression will always evaluate to {@code true}:<br>
     * {@code (this.value(this.enumNode(index)) == this.enumValue(index))}.
     *
     * @param index an index in a list of nodes with non-zero value
     * @return the specified node in the list of nodes with non-zero value
     *
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public int enumNode(int index) {
        checkSize(index);
        return node[this.index[index]];
    }

    /**
     * Returns the value of the specified node in a list of nodes with
     * non-zero value. Repeated invocations of this method with the same
     * parameter will return the same value if node values are not modified
     * between invocations. If {@code (index >= 0 && index < this.size())}, then
     * the following expression will always evaluate to {@code true}:<br>
     * {@code (this.value(this.enumNode(index)) == this.enumValue(index))}.
     *
     * @param index an index in a list of nodes with non-zero value
     * @return the value of the specified node in a list of nodes with
     * non-zero value
     *
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public float enumValue(int index) {
        checkSize(index);
        return value[this.index[index]];
    }

    /**
     * Returns the specified node value.
     *
     * @param node the first node
     * @return the specified node value
     * @throws IllegalArgumentException if {@code node < 0}
     */
    public float value(int node) {
        if (node < 0) {
            throw new IllegalArgumentException(String.valueOf(node));
        }
        return value[index(node)];
    }

    /**
     * Sets the value of each node to 0.
     */
    public void clear() {
        for (int j=0; j<this.size; ++j) {
            value[index[j]] = 0f;
        }
        size = 0;
    }

    /**
     * Returns a string representation of {@code this}. The exact
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
            sb.append(": node=");
            sb.append(enumNode(j));
            sb.append(" value=");
            sb.append(enumValue(j));
            sb.append(") ");
        }
        return sb.toString();
    }
}
