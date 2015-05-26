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
 * Class {@code HapNodes} stores nodes and associated values.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class HapNodes {

    private static final double loadFactor = 0.75;

    private int[] index;
    private int[] node;
    private double[] value;
    private int size;
    private int maxSize; // required to be a power of 2.
    private int rehashThreshold;


    /**
     * Creates a new instance of {@code HapNodes} that has an
     * initial value of 0 for each node.
     */
    public HapNodes() {
        this.size = 0;
        this.maxSize = (1<<10);
        this.rehashThreshold = (int) (loadFactor * maxSize);
        this.index = new int[maxSize];
        this.node = new int[maxSize];
        this.value = new double[maxSize];
    }

    /*
     * Return the storage index for specified node.  If the key is not
     * currently stored in the hash table, the index at which the value
     * should be stored is returned.
     */
    private int index(int node) {
        int h1 = node;
        int h2 = 97 + node;
        if ((h2 & 1)==0) {
            // h2 must be relatively prime to maxSize, which is a power of 2
            ++h2;
        }
        for (int k=0; k<maxSize; ++k) {
            int i = (h1 + k*h2) % maxSize;
            if (i<0) {
                i = -i;
            }
            if (value[i]==0.0 || (this.node[i]==node)) {
                return i;
            }
        }
        assert false;
        return -1;
    }

    /*
     * Increases the capacity of the internal hash table.
     */
    private void rehash() {
        assert this.size>=this.rehashThreshold;
        int newMaxSize = 2*maxSize;
        if (newMaxSize<0) {
            throw new IllegalStateException("hash table overflow");
        }
        int[] oldIndices = index;
        int[] oldNode = node;
        double[] oldValue = value;

        maxSize = newMaxSize;
        index = new int[newMaxSize];
        node = new int[newMaxSize];
        value = new double[newMaxSize];

        for (int j=0; j<size; ++j) {
            int oldIndex = oldIndices[j];
            int newIndex = index(oldNode[oldIndex]);
            index[j] = newIndex;
            node[newIndex] = oldNode[oldIndex];
            value[newIndex] = oldValue[oldIndex];
        }
        rehashThreshold = (int) (loadFactor * maxSize);
    }

    /**
     * Sets the value of the specified node to the maximum
     * of the node value immediately prior to method invocation
     * and the specified value.
     *
     * @param node the node.
     * @param value the value.
     *
     * @throws IllegalArgumentException if
     * {@code value<0.0 || Double.isNaN(value)}
     */
    public void maxUpdate(int node, double value) {
        if (value>0.0) {
            int i = index(node);
            if (this.value[i]>0.0) {
                if (value>this.value[i]) {
                    this.value[i] = value;
                }
            }
            else {
                this.index[size++] = i;
                this.node[i] = node;
                this.value[i] = value;
                if (this.size>=this.rehashThreshold) {
                    rehash();
                }
            }
        }
        else if (value>=0.0==false) {
            throw new IllegalArgumentException(String.valueOf(value));
        }
    }

    /**
     * Adds the specified value to the stored value of the specified
     * node.
     *
     * @param node the node.
     * @param value the value.
     *
     * @throws IllegalArgumentException if
     * {@code value<0.0 || Double.isNaN(value)}
     */
    public void sumUpdate(int node, double value) {
        if (value>0.0) {
            int i = index(node);
            if (this.value[i]>0.0) {
                this.value[i] += value;
            }
            else {
                this.index[size++] = i;
                this.node[i] = node;
                this.value[i] += value;
                if (this.size>=this.rehashThreshold) {
                    rehash();
                }
            }
        }
        else if (value>=0.0==false) {
            throw new IllegalArgumentException(String.valueOf(value));
        }
    }

    /**
     * Returns the number of nodes with non-zero value.
     * @return the number of nodes with non-zero value.
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
     * Returns the specified node in the list of nodes with non-zero value.
     *
     * @param index an index in the list of nodes with non-zero value.
     * @return the specified node in the list of nodes with non-zero value.
     *
     * @throws IndexOutOfBoundsException if
     * {@code index<0 || index>=this.size()}
     */
    public int enumNode(int index) {
        checkSize(index);
        return node[this.index[index]];
    }


    /**
     * Returns the value of the specified node in the list of nodes with
     * non-zero value.
     *
     * @param index an index in the list of nodes with non-zero value.
     * @return the value of the specified node in the list of nodes with
     * non-zero value.
     *
     * @throws IndexOutOfBoundsException if
     * {@code index<0 || index>=this.size()}
     */
    public double enumValue(int index) {
        checkSize(index);
        return value[this.index[index]];
    }

    /**
     * Returns the specified node value.
     *
     * @param node the first node.
     * @return the specified node value.
     */
    public double value(int node) {
        return value[index(node)];
    }

    /**
     * Sets the value of each node to 0.0.
     */
    public void clear() {
        for (int j=0; j<this.size; ++j) {
            value[index[j]] = 0.0;
        }
        size = 0;
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
        for (int j=0; j<size; ++j) {
            sb.append(" (");
            sb.append(j);
            sb.append(": node=");
            sb.append((int) enumNode(j));
            sb.append(" value=");
            sb.append(enumValue(j));
            sb.append(") ");
        }
        return sb.toString();
    }
}
