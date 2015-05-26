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
 * Class {@code SingleNodes} stores ordered node pairs and associated values.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class SingleNodes {

    private static final double loadFactor = 0.75;

    private int[] index;
    private int[] node1;
    private int[] node2;
    private double[] value;
    private int size;
    private int maxSize; // required to be a power of 2.
    private int rehashThreshold;


    /**
     * Creates a new instance of {@code SingleNodes} that has an
     * initial value of 0 for each ordered node pair.
     */
    public SingleNodes() {
        this.size = 0;
        this.maxSize = (1<<10);
        this.rehashThreshold = (int) (loadFactor * maxSize);
        this.index = new int[maxSize];
        this.node1 = new int[maxSize];
        this.node2 = new int[maxSize];
        this.value = new double[maxSize];
    }

    private static int hash1(int node1, int node2) {
        int hash = 5;
        hash = 71 * hash + node1;
        hash = 71 * hash + node2;
        return hash;
    }

    private static int hash2(int node1, int node2) {
        int hash = 7;
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
        int h1 = hash1(node1, node2);
        int h2 = hash2(node1, node2);
        if ((h2 & 1)==0) {
            // h2 must be relatively prime to maxSize, which is a power of 2
            ++h2;
        }
        for (int k=0; k<maxSize; ++k) {
            int i = (h1 + k*h2) % maxSize;
            if (i<0) {
                i = -i;
            }
            if (value[i]==0.0
                    || (this.node1[i]==node1 && this.node2[i]==node2)) {
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
        int[] oldNode1 = node1;
        int[] oldNode2 = node2;
        double[] oldValue = value;

        maxSize = newMaxSize;
        index = new int[newMaxSize];
        node1 = new int[newMaxSize];
        node2 = new int[newMaxSize];
        value = new double[newMaxSize];

        for (int j=0; j<size; ++j) {
            int oldIndex = oldIndices[j];
            int newIndex = index(oldNode1[oldIndex], oldNode2[oldIndex]);
            index[j] = newIndex;
            node1[newIndex] = oldNode1[oldIndex];
            node2[newIndex] = oldNode2[oldIndex];
            value[newIndex] = oldValue[oldIndex];
        }
        rehashThreshold = (int) (loadFactor * maxSize);
    }

    /**
     * Sets the value of the specified node pair to the maximum
     * of the node pair value immediately prior to method invocation
     * and the specified value.
     *
     * @param node1 the first node.
     * @param node2 the second node.
     * @param value the value.
     *
     * @throws IllegalArgumentException if
     * {@code value<0.0 || Double.isNaN(value)}
     */
    public void maxUpdate(int node1, int node2, double value) {
        if (value>0.0) {
            int i = index(node1, node2);
            if (this.value[i]>0.0) {
                if (value>this.value[i]) {
                    this.value[i] = value;
                }
            }
            else {
                this.index[size++] = i;
                this.node1[i] = node1;
                this.node2[i] = node2;
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
     * node pair.
     *
     * @param node1 the first node.
     * @param node2 the second node.
     * @param value the value.
     *
     * @throws IllegalArgumentException if
     * {@code value<0.0 || Double.isNaN(value)}
     */
    public void sumUpdate(int node1, int node2, double value) {
        if (value>0.0) {
            int i = index(node1, node2);
            if (this.value[i]>0.0) {
                this.value[i] += value;
            }
            else {
                this.index[size++] = i;
                this.node1[i] = node1;
                this.node2[i] = node2;
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
     * Returns the number of node pairs with non-zero value.
     * @return the number of node pairs with non-zero value.
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
     * node pairs with non-zero value.
     *
     * @param index an index in the list of node pairs with non-zero value.
     * @return the first node of the specified node pair in the list of
     * node pairs with non-zero value.
     *
     * @throws IndexOutOfBoundsException if
     * {@code index<0 || index>=this.size()}
     */
    public int enumNode1(int index) {
        checkSize(index);
        return node1[this.index[index]];
    }

    /**
     * Returns the second node of the specified node pair in the list of
     * node pairs with non-zero value.
     *
     * @param index an index in the list of node pairs with non-zero value.
     * @return the second node of the specified node pair in the list of
     * node pairs with non-zero value.
     *
     * @throws IndexOutOfBoundsException if
     * {@code index<0 || index>=this.size()}
     */
    public int enumNode2(int index) {
        checkSize(index);
        return node2[this.index[index]];
    }

    /**
     * Returns the value of the specified ordered node pair in the list of
     * node pairs with non-zero value.
     *
     * @param index an index in the list of node pairs with non-zero value.
     * @return the value of the specified ordered node pair in the list of
     * node pairs with non-zero value.
     *
     * @throws IndexOutOfBoundsException if
     * {@code index<0 || index>=this.size()}
     */
    public double enumValue(int index) {
        checkSize(index);
        return value[this.index[index]];
    }

    /**
     * Returns the specified node pair value.
     *
     * @param node1 the first node.
     * @param node2 the second node.
     * @return the specified node pair value.
     */
    public double value(int node1, int node2) {
        return value[index(node1, node2)];
    }

    /**
     * Sets the value of each node pair to 0.0.
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
            sb.append(": node1=");
            sb.append((int) enumNode1(j));
            sb.append(" node2=");
            sb.append((int) enumNode2(j));
            sb.append(" value=");
            sb.append(enumValue(j));
            sb.append(") ");
        }
        return sb.toString();
    }
}
