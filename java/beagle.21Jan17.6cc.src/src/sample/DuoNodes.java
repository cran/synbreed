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
 * <p>Class {@code DuoNodes} stores ordered node trios and associated values.
 * </p>
 * <p>Instances of class {@code DuoNodes} are not thread safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class DuoNodes {

    private static final float loadFactor = 0.75f;

    private int size;
    private int capacity; // required to be a power of 2
    private int rehashThreshold;

    private int[] index;
    private int[] nodeAB1;
    private int[] nodeA2;
    private int[] nodeB2;
    private float[] value;

    /**
     * Creates a new instance of {@code DuoNodes} that has an
     * initial value of 0 for each ordered node trio.
     */
    public DuoNodes() {
        this.size = 0;
        this.capacity = (1<<10);
        this.rehashThreshold = (int) (loadFactor * capacity);
        this.index = new int[capacity];
        this.nodeAB1 = new int[capacity];
        this.nodeA2 = new int[capacity];
        this.nodeB2 = new int[capacity];
        this.value = new float[capacity];
    }

    private static long hash1(int nodeAB1, int nodeA2, int nodeB2) {
        long hash = 5;
        hash = 71 * hash + nodeAB1;
        hash = 71 * hash + nodeA2;
        hash = 71 * hash + nodeB2;
        return hash;
    }

    private static long hash2(int nodeAB1, int nodeA2, int nodeB2) {
        long hash = 7;
        hash = 97 * hash + nodeAB1;
        hash = 97 * hash + nodeA2;
        hash = 97 * hash + nodeB2;
        return hash;
    }

    /*
     * Return the storage index for specified node trio.  If the key is not
     * currently stored in the hash table, the index at which the value
     * should be stored is returned.
     */
    private int index(int ab1, int a2, int b2) {
        long h1 = hash1(ab1, a2, b2);
        long h2 = hash2(ab1, a2, b2);
        if ((h2 & 1)==0) {
            // h2 must be relatively prime to maxSize, which is a power of 2
            ++h2;
        }
        long l = h1;
        for (int k=0; k<capacity; ++k) {
            int i = (int) (l % capacity);
            if (value[i]==0.0 ||
                    (nodeAB1[i]==ab1 && nodeA2[i]==a2 && nodeB2[i]==b2)) {
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
        int[] oldNodeAB1 = nodeAB1;
        int[] oldNodeA2 = nodeA2;
        int[] oldNodeB2 = nodeB2;
        float[] oldValue = value;

        capacity = newMaxSize;
        index = new int[newMaxSize];
        nodeAB1 = new int[newMaxSize];
        nodeA2 = new int[newMaxSize];
        nodeB2 = new int[newMaxSize];
        value = new float[newMaxSize];

        for (int j=0; j<size; ++j) {
            int oldInd = oldIndex[j];
            int newIndex = index(oldNodeAB1[oldInd], oldNodeA2[oldInd],
                    oldNodeB2[oldInd]);
            index[j] = newIndex;
            nodeAB1[newIndex] = oldNodeAB1[oldInd];
            nodeA2[newIndex] = oldNodeA2[oldInd];
            nodeB2[newIndex] = oldNodeB2[oldInd];
            value[newIndex] = oldValue[oldInd];
        }
        rehashThreshold = (int) (loadFactor * capacity);
    }

    /**
     * Adds the specified value to the stored value of the specified
     * node trio.
     *
     * @param nodeAB1 the first node
     * @param nodeA2 the second node
     * @param nodeB2 the third node
     * @param value the value
     *
     * @throws IllegalArgumentException if
     * {@code (nodeAB1 < 0 || nodeA2 < 0 || nodeB2 < 0)}
     * @throws IllegalArgumentException if
     * {@code value <= 0 || (Double.isFinite(value) == false)}
     */
    public void sumUpdate(int nodeAB1, int nodeA2, int nodeB2, float value) {
        if (nodeAB1 < 0) {
            throw new IllegalArgumentException(String.valueOf(nodeAB1));
        }
        if (nodeA2 < 0) {
            throw new IllegalArgumentException(String.valueOf(nodeA2));
        }
        if (nodeB2 < 0) {
            throw new IllegalArgumentException(String.valueOf(nodeB2));
        }
        if (value <= 0 || (Double.isFinite(value)==false) ) {
            throw new IllegalArgumentException(String.valueOf(value));
        }
        if (value>0.0) {
            int i = index(nodeAB1, nodeA2, nodeB2);
            boolean addNode = (this.value[i]==0f);
            this.value[i] += value;
            if (addNode) {
                this.index[size++] = i;
                this.nodeAB1[i] = nodeAB1;
                this.nodeA2[i] = nodeA2;
                this.nodeB2[i] = nodeB2;
                if (this.size>=this.rehashThreshold) {
                    rehash();
                }
            }
        }
    }

    /**
     * Returns the number of node trios with non-zero value.
     * @return the number of node trios with non-zero value
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
     * Returns the first node of the specified node trio in a list of
     * node trios with non-zero value. Repeated invocations of this
     * method with the same parameter will return the same value if
     * node values are not modified between invocations. If
     * {@code (index >= 0 && index < this.size())}, then the following
     * expression will always evaluate to {@code true}:<br>
     * {@code (this.value(this.enumNodeAB1(index), this.enumNodeA2(index),
     * this.enumNodeB2(index)) == this.enumValue(index))}.
     *
     * @param index an index in a list of node trios with non-zero value.
     * @return the first node of the specified node trio in a list of
     * node trios with non-zero value
     *
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public int enumNodeAB1(int index) {
        checkSize(index);
        return nodeAB1[this.index[index]];
    }

    /**
     * Returns the second node of the specified node trio in a list of
     * node trios with non-zero value. Repeated invocations of this
     * method with the same parameter will return the same value if
     * node values are not modified between invocations. If
     * {@code (index >= 0 && index < this.size())}, then the following
     * expression will always evaluate to {@code true}:<br>
     * {@code (this.value(this.enumNodeAB1(index), this.enumNodeA2(index),
     * this.enumNodeB2(index)) == this.enumValue(index))}.
     *
     * @param index an index in a list of node trios with non-zero value
     * @return the second node of the specified node trio in a list of
     * node trios with non-zero value
     *
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public int enumNodeA2(int index) {
        checkSize(index);
        return nodeA2[this.index[index]];
    }

    /**
     * Returns the third node of the specified node trio in a list of
     * node trios with non-zero value. Repeated invocations of this
     * method with the same parameter will return the same value if
     * node values are not modified between invocations. If
     * {@code (index >= 0 && index < this.size())}, then the following
     * expression will always evaluate to {@code true}:<br>
     * {@code (this.value(this.enumNodeAB1(index), this.enumNodeA2(index),
     * this.enumNodeB2(index)) == this.enumValue(index))}.
     *
     * @param index an index in a list of node trios with non-zero value
     * @return the third node of the specified node trio in a list of
     * node trios with non-zero value
     *
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public int enumNodeB2(int index) {
        checkSize(index);
        return nodeB2[this.index[index]];
    }

    /**
     * Returns the value of the specified ordered node trio in a list of
     * node trios with non-zero value. Repeated invocations of this
     * method with the same parameter will return the same value if
     * node values are not modified between invocations. If
     * {@code (index >= 0 && index < this.size())}, then the following
     * expression will always evaluate to {@code true}:<br>
     * {@code (this.value(this.enumNodeAB1(index), this.enumNodeA2(index),
     * this.enumNodeB2(index)) == this.enumValue(index))}.
     *
     * @param index an index in a list of node trios with non-zero value
     * @return the value of the specified ordered node trio in a list of
     * node trios with non-zero value
     *
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public float enumValue(int index) {
        checkSize(index);
        return value[this.index[index]];
    }

    /**
     * Returns the specified node trio value.
     *
     * @param nodeAB1 the first node
     * @param nodeA2 the second node
     * @param nodeB2 the third node
     * @return the specified node trio value
     * @throws IllegalArgumentException if
     * {@code (nodeAB1 < 0 || nodeA2 < 0 || nodeB2 < 0)}
     */
    public float value(int nodeAB1, int nodeA2, int nodeB2) {
        if (nodeAB1 < 0) {
            throw new IllegalArgumentException(String.valueOf(nodeAB1));
        }
        if (nodeA2 < 0) {
            throw new IllegalArgumentException(String.valueOf(nodeA2));
        }
        if (nodeB2 < 0) {
            throw new IllegalArgumentException(String.valueOf(nodeB2));
        }
        return value[index(nodeAB1, nodeA2, nodeB2)];
    }

    /**
     * Sets the value of each node trio to 0.
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
            sb.append(": nodeAB1=");
            sb.append(enumNodeAB1(j));
            sb.append(" nodeA2=");
            sb.append((int) enumNodeA2(j));
            sb.append(" nodeB2=");
            sb.append(enumNodeB2(j));
            sb.append(" value=");
            sb.append(enumValue(j));
            sb.append(") ");
        }
        return sb.toString();
    }
}
