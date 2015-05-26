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
 * Class {@code TrioNodes} stores ordered node quartets and associated values.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class TrioNodes {

    private final double loadFactor = 0.75;

    private int size;
    private int maxSize; // required to be a power of 2.
    private int rehashThreshold;

    private int[] index;
    private int[] nodeA1;
    private int[] nodeA2;
    private int[] nodeB1;
    private int[] nodeB2;
    private double[] value;

    /**
     * Creates a new instance of {@code TrioNodes} that has an
     * initial value of 0 for each ordered node quartet.
     */
    public TrioNodes() {
        this.size = 0;
        this.maxSize = (1<<10);
        this.rehashThreshold = (int) (loadFactor * maxSize);
        this.index = new int[maxSize];
        this.nodeA1 = new int[maxSize];
        this.nodeA2 = new int[maxSize];
        this.nodeB1 = new int[maxSize];
        this.nodeB2 = new int[maxSize];
        this.value = new double[maxSize];
    }

    private static int hash1(int nodeA1, int nodeA2, int nodeB1, int nodeB2) {
        int hash = 5;
        hash = 71 * hash + nodeA1;
        hash = 71 * hash + nodeA2;
        hash = 71 * hash + nodeB1;
        hash = 71 * hash + nodeB2;
        return hash;
    }

    private static int hash2(int nodeA1, int nodeA2, int nodeB1, int nodeB2) {
        int hash = 7;
        hash = 97 * hash + nodeA1;
        hash = 97 * hash + nodeA2;
        hash = 97 * hash + nodeB1;
        hash = 97 * hash + nodeB2;
        return hash;
    }

    /*
     * Return the storage index for specified node trio.  If the key is not
     * currently stored in the hash table, the index at which the value
     * should be stored is returned.
     */
    private int index(int a1, int a2, int b1, int b2) {
        int h1 = hash1(a1, a2, b1, b2);
        int h2 = hash2(a1, a2, b1, b2);
        if ((h2 & 1)==0) {
            // h2 must be relatively prime to maxSize, which is a power of 2
            ++h2;
        }
        for (int k=0; k<maxSize; ++k) {
            int i = (h1 + k*h2) % maxSize;
            if (i<0) {
                i = -i;
            }
            if (value[i]==0.0 ||
                    (nodeA1[i]==a1 && nodeA2[i]==a2
                    && nodeB1[i]==b1 && nodeB2[i]==b2)) {
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
        assert size>=this.rehashThreshold;
        int newMaxSize = 2*maxSize;
        if (newMaxSize<0) {
            throw new IllegalStateException("hash table overflow");
        }
        int[] oldIndices = index;
        int[] oldNodeA1 = nodeA1;
        int[] oldNodeA2 = nodeA2;
        int[] oldNodeB1 = nodeB1;
        int[] oldNodeB2 = nodeB2;
        double[] oldValue = value;

        maxSize = newMaxSize;
        index = new int[newMaxSize];
        nodeA1 = new int[newMaxSize];
        nodeA2 = new int[newMaxSize];
        nodeB1 = new int[newMaxSize];
        nodeB2 = new int[newMaxSize];
        value = new double[newMaxSize];

        for (int j=0; j<size; ++j) {
            int oldIndex = oldIndices[j];
            int newIndex = index(oldNodeA1[oldIndex], oldNodeA2[oldIndex],
                    oldNodeB1[oldIndex], oldNodeB2[oldIndex]);
            index[j] = newIndex;
            nodeA1[newIndex] = oldNodeA1[oldIndex];
            nodeA2[newIndex] = oldNodeA2[oldIndex];
            nodeB1[newIndex] = oldNodeB1[oldIndex];
            nodeB2[newIndex] = oldNodeB2[oldIndex];
            value[newIndex] = oldValue[oldIndex];
        }
        rehashThreshold = (int) (loadFactor * maxSize);
    }

    /**
     * Sets the value of the specified node quarter to the maximum
     * of the node quartet value immediately prior to method invocation
     * and the specified value.
     *
     * @param nodeA1 the first node.
     * @param nodeA2 the second node.
     * @param nodeB1 the third node.
     * @param nodeB2 the fourth node.
     * @param value the value.
     *
     * @throws IllegalArgumentException if
     * {@code value<0.0 || Double.isNaN(value)}
    */
    public void maxUpdate(int nodeA1, int nodeA2, int nodeB1, int nodeB2,
            double value) {
        if (value>0.0) {
            int i = index(nodeA1, nodeA2, nodeB1, nodeB2);
            if (this.value[i]>0.0) {
                if (value>this.value[i]) {
                    this.value[i] = value;
                }
            }
            else {
                this.index[size++] = i;
                this.nodeA1[i] = nodeA1;
                this.nodeA2[i] = nodeA2;
                this.nodeB1[i] = nodeB1;
                this.nodeB2[i] = nodeB2;
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
     * node quartet.
     *
     * @param nodeA1 the first node.
     * @param nodeA2 the second node.
     * @param nodeB1 the third node.
     * @param nodeB2 the fourth node.
     * @param value the value.
     *
     * @throws IllegalArgumentException if
     * {@code value<0.0 || Double.isNaN(value)}
     */
    public void sumUpdate(int nodeA1, int nodeA2, int nodeB1, int nodeB2,
            double value) {
        if (value>0.0) {
            int i = index(nodeA1, nodeA2, nodeB1, nodeB2);
            if (this.value[i]>0.0) {
                this.value[i] += value;
            }
            else {
                this.index[size++] = i;
                this.nodeA1[i] = nodeA1;
                this.nodeA2[i] = nodeA2;
                this.nodeB1[i] = nodeB1;
                this.nodeB2[i] = nodeB2;
                this.value[i] += value;
                if (this.size>this.rehashThreshold) {
                    rehash();
                }
            }
        }
        else if (value>=0.0==false) {
            throw new IllegalArgumentException(String.valueOf(value));
        }
    }

    /**
     * Returns the number of node quartets with non-zero value.
     * @return the number of node quartets with non-zero value.
     */
    public int size() {
        return size;
    }

    private void checkSize(int index) {
        if (index >= size()) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
    }

    /**
     * Returns the first node of the specified node quartet in the list of
     * node quartets with non-zero value.
     *
     * @param index an index in the list of node quartets with non-zero value.
     * @return the first node of the specified node quartet in the list of
     * node quartets with non-zero value.
     *
     * @throws IndexOutOfBoundsException if
     * {@code index<0 || index>=this.size()}
     */
    public int enumNodeA1(int index) {
        checkSize(index);
        return nodeA1[this.index[index]];
    }

    /**
     * Returns the second node of the specified node quartet in the list of
     * node quartets with non-zero value.
     *
     * @param index an index in the list of node quartets with non-zero value.
     * @return the second node of the specified node quartet in the list of
     * node quartets with non-zero value.
     *
     * @throws IndexOutOfBoundsException if
     * {@code index<0 || index>=this.size()}
     */
    public int enumNodeA2(int index) {
        checkSize(index);
        return nodeA2[this.index[index]];
    }

    /**
     * Returns the third node of the specified node quartet in the list of
     * node quartets with non-zero value.
     *
     * @param index an index in the list of node quartets with non-zero value.
     * @return the third node of the specified node quartet in the list of
     * node quartets with non-zero value.
     *
     * @throws IndexOutOfBoundsException if
     * {@code index<0 || index>=this.size()}
     */
    public int enumNodeB1(int index) {
        checkSize(index);
        return nodeB1[this.index[index]];
    }

    /**
     * Returns the fourth node of the specified node quartet in the list of
     * node quartets with non-zero value.
     *
     * @param index an index in the list of node quartets with non-zero value.
     * @return the fourth node of the specified node quartet in the list of
     * node quartets with non-zero value.
     *
     * @throws IndexOutOfBoundsException if
     * {@code index<0 || index>=this.size()}
     */
    public int enumNodeB2(int index) {
        checkSize(index);
        return nodeB2[this.index[index]];
    }

    /**
     * Returns the value of the specified ordered node quartet in the list of
     * ordered node quartets with non-zero value.
     *
     * @param index an index in the list of node quartets with non-zero value.
     * @return the value of the specified ordered node quartet in the list of
     * ordered node quartets with non-zero value.
     *
     * @throws IndexOutOfBoundsException if
     * {@code index<0 || index>=this.size()}
     */
    public double enumValue(int index) {
        checkSize(index);
        return value[this.index[index]];
    }

    /**
     * Returns the specified ordered node quartet value.
     *
     * @param nodeA1 the first node index.
     * @param nodeA2 the second node index.
     * @param nodeB1 the third node index.
     * @param nodeB2 the fourth node index.
     * @return the specified ordered node quartet value.
     */
    public double value(int nodeA1, int nodeA2, int nodeB1, int nodeB2) {
        return value[index(nodeA1, nodeA2, nodeB1, nodeB2)];
    }

    /**
     * Sets the value of each node quartet to 0.0.
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
            sb.append(": nodeA1=");
            sb.append((int) enumNodeA1(j));
            sb.append(" nodeA2=");
            sb.append((int) enumNodeA2(j));
            sb.append(" nodeB1=");
            sb.append((int) enumNodeB1(j));
            sb.append(" nodeB2=");
            sb.append((int) enumNodeB2(j));
            sb.append(" value=");
            sb.append(enumValue(j));
            sb.append(") ");
        }
        return sb.toString();
    }
}
