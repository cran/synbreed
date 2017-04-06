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
package dag;

import blbutil.Const;
import blbutil.IntSet;
import java.util.Arrays;

/**
 * <p>Class {@code HighCapacityDagLevel} represents a level of a leveled
 * directed acyclic graph (DAG) that can contain up to
 * {@code Integer.MAX_VALUE} edges.
 * </p>
 * <p>Instances of {@code HighCapacityDagLevel} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class HighCapacityDagLevel implements DagLevel {

/*
 * The k-th edge parent node index is stored in {@code this.parentNodes[k]}.
 * The k-th edge child node index is stored in {@code this.childNodes[k]}.
 * The k-th edge symbol is stored in {@code this.symbols[k]}.
 * The k-th edge count is stored in {@code this.edgeCounts[k]}.
 * The k-th edge conditional edge probability is stored in
 * {@code this.condEdgeProbs[k]}, and is defined to be the
 * k-th edge count divided by the k-th edge's parent node count.
 * The k-th node count is stored in {@code this.nodeCounts[k]}.
 *
 * The outgoing edges indices of the k-th parent node are stored in consecutive
 * entries of {@code this.parents} beginning with
 * {@code this.parentIndices[k]} (inclusive) and ending with
 * {@code this.parentIndices[k+1]} (exclusive).
 *
 * The ingoing edges indices of the k-th child node are stored in consecutive
 * entries of {@code this.children} beginning with
 * {@code this.childIndices[k]} (inclusive) and ending with
 * {@code this.childIndices[k+1]} (exclusive).
 */
    private final float count;
    private final int[] parentNodes;
    private final int[] childNodes;
    private final int[] parentIndices;
    private final int[] parents;
    private final int[] childIndices;
    private final int[] children;
    private final int[] symbols;
    private final float[] edgeCounts;
    private final float[] condEdgeProbs;
    private final float[] parentCounts;

    /**
     * Constructs a new {@code HighCapacityDagLevel} instance from the
     * specified data.
     *
     * @param parentNodes an array mapping edge index to parent node index
     * @param childNodes an array mapping edge index to child node index
     * @param symbols an array mapping edge index to the symbol labeling the
     * edge
     * @param counts an array mapping edge index to edge count
     *
     * @throws IllegalArgumentException if the specified arrays do not all
     * have the same length
     * @throws IllegalArgumentException if any element of the symbols array
     * is negative
     * @throws IllegalArgumentException if any two edges have the same
     * parent node and are both labeled with the same symbol
     * @throws IllegalArgumentException if the set of values of the
     * {@code parentNodes} array is not equal to {@code {0, 1, 2, ..., k}} for
     * some {@code k}
     * @throws IllegalArgumentException if the set of values of the
     * {@code childNodes} array is not equal to {@code {0, 1, 2, ..., k}}
     * for some {@code k}
     *
     * @throws NullPointerException if any parameter is {@code null}
     */
    public HighCapacityDagLevel(int[] parentNodes, int[] childNodes,
            int[] symbols, float[] counts) {
        int nEdges = checkLengths(parentNodes, childNodes, symbols, counts);
        this.parentIndices = getIndicesArray(parentNodes);
        this.childIndices = getIndicesArray(childNodes);
        this.parentNodes = parentNodes.clone();
        this.childNodes = childNodes.clone();
        this.symbols = symbols.clone();
        this.edgeCounts = counts.clone();
        this.condEdgeProbs = new float[nEdges];
        this.parents = new int[nEdges];
        this.children = new int[nEdges];

        int[] pIndices = Arrays.copyOfRange(parentIndices, 0,
                parentIndices.length-1);
        int[] cIndices = Arrays.copyOfRange(childIndices, 0,
                childIndices.length-1);
        this.parentCounts = parentCnts(parentNodes, counts, pIndices.length);
        this.count = sum(this.parentCounts);

        for (int j=0; j<nEdges; ++j) {
            int p = parentNodes[j];
            int c = childNodes[j];
            this.parents[pIndices[p]++] = j;
            this.children[cIndices[c]++] = j;
            this.condEdgeProbs[j] = counts[j] / parentCounts[p];
        }
        checkForDuplicateOutEdges(parentIndices, parents, symbols);
    }

    private static int checkLengths(int[] parentNodes, int[] childNodes,
            int[] symbols, float[] counts) {
        if ( ((parentNodes.length != childNodes.length)
                || (parentNodes.length != symbols.length))
                || (parentNodes.length != counts.length) ) {
            throw new IllegalArgumentException("inconsistent arrays");
        }
        return parentNodes.length;
    }

    private static void checkForDuplicateOutEdges(int[] parentIndices,
            int[] parents, int[] symbols) {
        IntSet indexSet = new IntSet(symbols.length);
        for (int j=1; j<parentIndices.length; ++j) {
            indexSet.clear();
            for (int k=parentIndices[j-1], n=parentIndices[j]; k<n; ++k) {
                int edge = parents[k];
                if (indexSet.add(symbols[edge])==false) {
                    throw new IllegalArgumentException("duplicate edge");
                }
            }
        }
    }

    private static int[] getIndicesArray(int[] nodes) {
        int[] countArray = elementCounts(nodes);
        int[] indicesArray = new int[countArray.length + 1];
        for (int j=1; j<indicesArray.length; ++j) {
            assert countArray[j-1]>0;
            indicesArray[j] = indicesArray[j-1] + countArray[j-1];
        }
        return indicesArray;
    }

    /*
     * Returns an array of length {@code max(nodes) + 1}
     * whose {@code j}-th element is the number of
     * elements of the specified array that have value {@code j}.
     *
     * @param nodes an array of non-negative values.
     * @return an array whose {@code j}-th element is the number of
     * elements of the specified array that have value {@code j}.
     * @throws IllegalArgumenException if set of elements of the
     * specified array is not equal to {@code {0, 1, 2, ..., k}} for some
     * {@code k}.
     */
    private static int[] elementCounts(int[] array) {
        int maxNode = max(array);
        int[] nodeCounts = new int[maxNode + 1];
        for (int i : array) {
            ++nodeCounts[i];
        }
        for (int j=0; j<nodeCounts.length; ++j) {
            if (nodeCounts[j]==0) {
                throw new IllegalArgumentException("no element with value " + j);
            }
        }
        return nodeCounts;
    }

    private static int max(int[] ia) {
        int max = 0;
        for (int i : ia) {
            if (i>max) {
                max=i;
            }
        }
        return max;
    }

    private float sum(float[] fa) {
        float sum = 0.0f;
        for (float f : fa) {
            sum += f;
        }
        return sum;
    }

    private float[] parentCnts(int[] parentNodes, float[] counts, int nNodes) {
        float[] parentCnts = new float[nNodes];
        for (int j=0; j<condEdgeProbs.length; ++j) {
            int p = parentNodes[j];
            parentCnts[p] += counts[j];
        }
        return parentCnts;
    }

    @Override
    public int nEdges() {
        return condEdgeProbs.length;
    }

    @Override
    public int nParentNodes() {
        return parentIndices.length - 1;
    }

    @Override
    public int nChildNodes() {
        return childIndices.length - 1;
    }

    @Override
    public int parentNode(int edge) {
        return parentNodes[edge];
    }

    @Override
    public int childNode(int edge) {
        return childNodes[edge];
    }

    @Override
    public int symbol(int edge) {
        return symbols[edge];
    }

    @Override
    public float edgeWeight(int edge) {
        return edgeCounts[edge];
    }

    @Override
    public float parentWeight(int parentNode) {
        return parentCounts[parentNode];
    }

    @Override
    public float condEdgeProb(int edge) {
        return condEdgeProbs[edge];
    }

    @Override
    public float edgeProb(int edge) {
        return (edgeCounts[edge] / count);
    }

    @Override
    public float parentProb(int node) {
        return (parentCounts[node] / count);
    }

    @Override
    public int nOutEdges(int parentNode) {
        return (parentIndices[parentNode+1] - parentIndices[parentNode]);
    }

    @Override
    public int outEdge(int parentNode, int outEdgeIndex) {
        if (outEdgeIndex<0 || outEdgeIndex>=nOutEdges(parentNode)) {
            throw new IndexOutOfBoundsException(String.valueOf(outEdgeIndex));
        }
        return (parents[parentIndices[parentNode] + outEdgeIndex]);
    }

    @Override
    public int outEdgeBySymbol(int parentNode, int symbol) {
        int start = parentIndices[parentNode];
        int end = parentIndices[parentNode+1];
        for (int j=start; j<end; ++j) {
            int edgeIndex = parents[j];
            if (symbols[edgeIndex]==symbol) {
                return edgeIndex;
            }
        }
        return -1;
    }

    @Override
    public int nInEdges(int childNode) {
        return (childIndices[childNode+1] - childIndices[childNode]);
    }

    @Override
    public int inEdge(int childNode, int inEdgeIndex) {
        if (inEdgeIndex<0 || inEdgeIndex>=nInEdges(childNode)) {
            throw new IndexOutOfBoundsException(String.valueOf(inEdgeIndex));
        }
        return this.children[childIndices[childNode] + inEdgeIndex];
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(500);
        sb.append(Const.nl);
        sb.append("parentNodes=");
        sb.append(Arrays.toString(parentNodes));
        sb.append(Const.nl);
        sb.append("childNodes=");
        sb.append(Arrays.toString(childNodes));
        sb.append(Const.nl);
        sb.append("symbols=");
        sb.append(Arrays.toString(symbols));
        sb.append(Const.nl);
        sb.append("condEdgeProbs=");
        sb.append(Arrays.toString(condEdgeProbs));
        sb.append(Const.nl);
        sb.append("parentCounts=");
        sb.append(Arrays.toString(parentCounts));
        sb.append(Const.nl);
        sb.append("edgeCounts=");
        sb.append(Arrays.toString(edgeCounts));
        sb.append(Const.nl);
        sb.append("parentIndices=");
        sb.append(Arrays.toString(parentIndices));
        sb.append(Const.nl);
        sb.append("parents=");
        sb.append(Arrays.toString(parents));
        sb.append(Const.nl);
        sb.append("childIndices=");
        sb.append(Arrays.toString(childIndices));
        sb.append(Const.nl);
        sb.append("children=");
        sb.append(Arrays.toString(children));
        return sb.toString();
    }
}
