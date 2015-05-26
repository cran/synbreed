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
import haplotype.HapsMarker;
import java.util.Arrays;
import vcf.Marker;

/**
 * Class {@code MergeableDagLevel} represents a level of a leveled
 * directed acyclic graph (DAG).  The class includes a public method for
 * merging parent nodes.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class MergeableDagLevel {

    private MergeableDagLevel nextLevel = null;
    private MergeableDagLevel prevLevel = null;

    private final Marker marker;
    private final int markerIndex;
    private final int nAlleles;
    private final int nHaps;

    private int[][] outEdges;  // [allele][parent node]
    private int[] child2FirstInEdge;
    private int[] inEdge2NextInEdge;

    private int[] parentNodes;    // edge -> parent node
    private int[] childNodes;     // edge -> child node
    private byte[] symbols;       // edge -> symbol
    private float[] counts;       // edge -> weight

    private int[] child2FirstHap; // child node -> first hap index
    private int[] hap2NextHap;   // current hap index -> next hap index

    /**
     * Constructs a new {@code MergeableDagLevel} instance from the specified
     * phased genotype data and haplotype weights.  The {@code previous()}
     * method of the constructed instance will return {@code null}.
     * @param data the phased genotype data.
     * @param weights an array mapping haplotype indices to non-negative
     * weights.
     *
     * @throws IllegalArgumentException if {@code weights.length!=data.nHaps()}
     * @throws NullPointerException if {@code data==null || weights==null}
     */
    public MergeableDagLevel(HapsMarker data, float[] weights) {
        checkParameters(data, weights);
        boolean isRootLevel = true;
        this.prevLevel = null;
        this.nextLevel = null;
        this.marker = data.marker();
        this.markerIndex = 0;
        this.nAlleles = data.marker().nAlleles();
        this.nHaps = data.nHaps();
        allocateAndInitializeArrays(isRootLevel, nAlleles, nHaps);
        fillArrays(data, weights);
    }

    /**
     * Constructs a new {@code MergeableDagLevel} instance with the
     * specified previous {@code MergeableDagLevel}, phased genotype data,
     * and haplotype weights.  This constructor does not alter any field
     * of the specified {@code prevLevel} object.
     * @param prevLevel the previous {@code MergeableDagLevel}.
     * @param data the phased genotype data.
     * @param weights an array mapping haplotype indices to non-negative
     * weights.
     *
     * @throws IllegalArgumentException
     *   if {@code weights.length!= data.nHaps()}
     * @throws IllegalArgumentException if {@code prevLevel.nextLevel()!=null}
     * @throws IllegalArgumentException if {@code parent.nHaps()!=data.nHaps()}
     * @throws NullPointerException if
     * {@code parent==null || data==null || weights==null}
     */
    public MergeableDagLevel(MergeableDagLevel prevLevel, HapsMarker data,
            float[] weights) {
        checkParameters(prevLevel, data, weights);
        boolean isRootLevel = false;
        this.prevLevel = prevLevel;
        this.nextLevel = null;
        this.marker = data.marker();
        this.markerIndex = prevLevel.markerIndex() + 1;
        this.nAlleles = data.marker().nAlleles();
        this.nHaps = data.nHaps();
        allocateAndInitializeArrays(isRootLevel, nAlleles, nHaps);
        fillArrays(prevLevel, data, weights);
    }

    private void checkParameters(HapsMarker data, float[] weights) {
        if (weights.length != data.nHaps()) {
             String s = "data.nHaps()=" + data.nHaps()
                    + " != weights.length=" + weights.length;
             throw new IllegalArgumentException(s);
        }
    }

    private void checkParameters(MergeableDagLevel parent, HapsMarker data,
            float[] weights) {
        checkParameters(data, weights);
        if (parent.nextLevel!=null) {
            throw new IllegalArgumentException("parent.nextLevel!=null");
        }
        if (parent.nHaps()!=data.nHaps()) {
            throw new IllegalArgumentException("inconsistent samples");
        }
        // NB: the sequences of sample ID indices are not checked
    }

    private void allocateAndInitializeArrays(boolean isRootLevel, int nAlleles,
            int nHaps) {
        int size = isRootLevel ? nAlleles : nHaps;
        this.outEdges = new int[nAlleles][size];
        this.child2FirstInEdge = new int[size];
        this.inEdge2NextInEdge = new int[size];
        this.parentNodes = new int[size];
        this.childNodes = new int[size];
        this.symbols = new byte[size];
        this.counts = new float[size];
        this.child2FirstHap = new int[nHaps];
        this.hap2NextHap = new int[nHaps];

        for (int[] oe : outEdges) {
            Arrays.fill(oe, -1);
        }
        Arrays.fill(child2FirstInEdge, -1);
        Arrays.fill(inEdge2NextInEdge, -1);
        Arrays.fill(parentNodes, -1);
        Arrays.fill(childNodes, -1);
        Arrays.fill(symbols, (byte) -1);
        Arrays.fill(child2FirstHap, -1);
        Arrays.fill(hap2NextHap, -1);
    }

    private void fillArrays(HapsMarker data, float[] weights) {
        int parentNode = 0;
        for (int hap=0, n=data.nHaps(); hap<n; ++hap) {
            byte symbol = data.allele(hap);
            float count = weights[hap];
            int edge = this.outEdges[symbol][parentNode];
            if (edge == -1) {
                edge = symbol;
                addEdge(parentNode, symbol, count, edge, hap);
            }
            else {
                assert edge == symbol;
                assert edge==childNodes[edge];
                int child = childNodes[edge];
                this.counts[edge] += count;
                this.hap2NextHap[hap] = this.child2FirstHap[child];
                this.child2FirstHap[child] = hap;
            }
        }
    }

    private void fillArrays(MergeableDagLevel prevLevel, HapsMarker data,
            float[] weights) {
        int nEdges = 0;
        for (int node=0, n=prevLevel.child2FirstHap.length; node<n; ++node) {
            if (prevLevel.child2FirstHap[node] >= 0) {
                int hap = prevLevel.child2FirstHap[node];
                while (hap != -1) {
                    byte symbol = data.allele(hap);
                    float count = weights[hap];
                    int edge = this.outEdges[symbol][node];
                    if (edge == -1) {
                        addEdge(node, symbol, count, nEdges++, hap);
                    }
                    else {
                        assert edge==childNodes[edge];
                        int child = childNodes[edge];
                        this.counts[edge] += count;
                        this.hap2NextHap[hap] = this.child2FirstHap[child];
                        this.child2FirstHap[child] = hap;
                    }
                    hap = prevLevel.hap2NextHap[hap];
                }
            }
        }
        if (nEdges < 0.75*nHaps) {
            reduceEdgeArrayLengths(nEdges);
        }
        prevLevel.removeHaplotypeIndices();
    }

    private void addEdge(int parentNode, byte symbol, float weight,
            int edge, int haplotype) {
        int childNode = edge;
        outEdges[symbol][parentNode] = edge;
        child2FirstInEdge[childNode] = edge;
        parentNodes[edge] = parentNode;
        childNodes[edge] = childNode;
        symbols[edge] = symbol;
        counts[edge] = weight;
        child2FirstHap[childNode] = haplotype;
    }

    private void reduceEdgeArrayLengths(int newLength) {
        child2FirstInEdge = Arrays.copyOf(child2FirstInEdge, newLength);
        inEdge2NextInEdge = Arrays.copyOf(inEdge2NextInEdge, newLength);
        parentNodes = Arrays.copyOf(parentNodes, newLength);
        childNodes = Arrays.copyOf(childNodes, newLength);
        symbols = Arrays.copyOf(symbols, newLength);
        counts = Arrays.copyOf(counts, newLength);
    }

    /**
     * Removes haplotype index data from {@code this}.
     */
    private void removeHaplotypeIndices() {
        this.child2FirstHap = null;
        this.hap2NextHap = null;
    }

   /**
     * Sets the previous DAG level to {@code null}, and returns
     * the previous DAG level that existed immediately prior to the invocation
     * of this method.
     * @return the previous DAG level that existed immediately prior to the
     * invocation of this method.
     */
    public MergeableDagLevel setPreviousToNull() {
        MergeableDagLevel prev = this.prevLevel;
        this.prevLevel = null;
        return prev;
    }

    /**
     * Sets the next level to the specified {@code MergeableDagLevel}.
     * @param nextLevel the next level.
     * @throws IllegalArgumentException if
     * {@code nextLevel.previousLevel()!=this}
     */
    public void setNextLevel(MergeableDagLevel nextLevel) {
        if (nextLevel.prevLevel != this) {
            throw new IllegalArgumentException("nextLevel.previousLevel!=this");
        }
        this.nextLevel = nextLevel;
    }

    /**
     * Returns the previous DAG level.
     * @return the previous DAG level.
     */
    public MergeableDagLevel previous() {
        return prevLevel;
    }

    /**
     * Returns the next DAG level.
     * @return the next DAG level.
     */
    public MergeableDagLevel next() {
        return nextLevel;
    }

    /**
     * Returns {@code true} if the specified parent node has a
     * sibling  and returns {@code false} otherwise.
     * Two parent nodes are siblings if they are connected by an
     * edge to the same parent node at the previous level of the DAG.
     *
     * @param parentNode a parent node index.
     * @return {@code true} if the specified parent node has a
     * sibling  and returns {@code false} otherwise.
     */
    public boolean hasSibling(int parentNode) {
        int edge = prevLevel.child2FirstInEdge[parentNode];
        while (edge>=0) {
            int pn = prevLevel.parentNodes[edge];
            int cnt = 0;
            for (int allele=0, n=prevLevel.nAlleles; allele<n; ++allele) {
                if (prevLevel.outEdges[allele][pn]>=0) {
                    ++cnt;
                }
            }
            if (cnt>1) {
                return true;
            }
            edge = prevLevel.inEdge2NextInEdge[edge];
        }
        return false;
    }

    /**
     * Returns an immutable {@code DagLevel} corresponding to
     * {@code this}. The parent node, edge, and child node indices
     * in the returned {@code DagLevel} are the ranks of the
     * parent node, edge, and child node indices for {@code this},
     * with rank 0 corresponding to the smallest index.
     * @return an immutable {@code DagLevel} corresponding to
     * {@code this}.
     */
    public DagLevel toDagLevel() {
         char[] modParentNodes = rankValues(
                DagUtils.removeValues(parentNodes, -1));
         char[] modChildNodes = rankValues(
                DagUtils.removeValues(childNodes, -1));
         byte[] modSymbols = DagUtils.removeValues(symbols, (byte) -1);
         float[] modCounts = DagUtils.removeValues(counts, 0.0f);
         return new ImmutableDagLevel(marker, modParentNodes, modChildNodes,
                 modSymbols, modCounts);
    }

    /*
     * Returns an array obtained by replacing each array value with it's
     * rank when the set of array values is ordered: the smallest value
     * is replaced by 0, the next smallest value is replaced by 1, etc.
     *
     * @throws NullPointerException if {@code array==null}
     * @throws IllegalArgumentException if {@code array==null}
     * @throws IllegalArgumentException if any element of array
     * is negative.
     * @throws NegativeArrayException if any element of array equals
     * {@code Integer.MAX_VALUE}.
     */
    private static char[] rankValues(int[] array) {
        if (array.length==0) {
            throw new IllegalArgumentException("array.length==0");
        }
        assert array.length < Character.MAX_VALUE;
        int[] sortedCopy = array.clone();
        Arrays.sort(sortedCopy);
        if (sortedCopy[0] < 0) {
            String s = "element<0: " + sortedCopy[0];
            throw new IllegalArgumentException(s);
        }
        int n = sortedCopy[sortedCopy.length - 1] + 1;
        char[] indexMap = new char[n];
        Arrays.fill(indexMap, Character.MAX_VALUE);
        char index = 0;
        indexMap[sortedCopy[0]] = index++;
        for (int j=1; j<sortedCopy.length; ++j) {
            if (sortedCopy[j] != sortedCopy[j-1]) {
                indexMap[sortedCopy[j]] = index++;
            }
        }
        char[] transformedArray = new char[array.length];
        for (int j=0; j<transformedArray.length; ++j) {
            transformedArray[j] = indexMap[array[j]];
        }
        return transformedArray;
    }

    /**
     * Merges the two specified parent nodes and assigns the merged
     * node to the specified {@code retainedNode} index.
     *
     * @param retainedNode a parent node which will receive ingoing and
     * outgoing edges of {@code removedNode}.
     * @param removedNode a parent node that will be deleted after merging.
     *
     * @throws IllegalArgumentException if {@code retainedNode}
     * or {@code returnedNode} is not a valid parent node index.
     */
    public void mergeParentNodes(int retainedNode, int removedNode) {
        if (isParentNode(retainedNode)==false) {
            String s = "invalid parent node: " + retainedNode;
            throw new IllegalArgumentException(s);
        }
        if (isParentNode(removedNode)==false) {
            String s = "invalid parent node: " + removedNode;
            throw new IllegalArgumentException(s);
        }
        prevLevel.mergeChildNodes(retainedNode, removedNode);
        mergeParentNodes2(retainedNode, removedNode);
    }

    private void mergeParentNodes2(int retainedNode, int removedNode) {
        for (byte j=0; j<nAlleles; ++j) {
            int retainedEdge = outEdges[j][retainedNode];
            int removedEdge = outEdges[j][removedNode];
            if (removedEdge >= 0) {
                if (retainedEdge == -1) {
                    changeParent(removedEdge, retainedNode);
                }
                else {
                    int retainedChild = childNode(retainedEdge);
                    int removedChild = childNode(removedEdge);
                    mergeEdges(retainedEdge, removedEdge);
                    if (nextLevel != null) {
                        nextLevel.mergeParentNodes2(retainedChild, removedChild);
                    }
                }
            }
        }
    }

    /*
     * Merges the two specified child nodes and assigns the merged
     * node to the specified {@code retainedNode} index.  Ingoing edges
     * to {@code removedNode} are redirected to be ingoing edges
     * to {@code retainedNode}.
     *
     * @param retainedNode a child node which will receive ingoing edges of
     * {@code removedNode}.
     * @param removedNode a child node that will be deleted after merging.
     */
    private void mergeChildNodes(int retainedNode, int removedNode) {
        int lastEdge = -1;
        int edge = child2FirstInEdge[removedNode];
        while (edge != -1) {
            assert childNodes[edge] == removedNode;
            childNodes[edge] = retainedNode;
            lastEdge = edge;
            edge = inEdge2NextInEdge[edge];
        }
        if (lastEdge != -1) {
            inEdge2NextInEdge[lastEdge] = child2FirstInEdge[retainedNode];
            child2FirstInEdge[retainedNode] = child2FirstInEdge[removedNode];
            child2FirstInEdge[removedNode] = -1;
        }
    }

    private void changeParent(int edge, int newParent) {
        int oldParent = parentNodes[edge];
        byte symbol = symbols[edge];
        assert (outEdges[symbol][oldParent] == edge);
        assert (outEdges[symbol][newParent] == -1);
        outEdges[symbol][oldParent] = -1;
        outEdges[symbol][newParent] = edge;
        parentNodes[edge] = newParent;
    }

    private void mergeEdges(int retainedEdge, int removedEdge) {
        assert symbols[retainedEdge] == symbols[removedEdge];
        assert counts[removedEdge] > 0.0f;
        counts[retainedEdge] += counts[removedEdge];
        if (nextLevel==null) {
            mergeHaplotypes(childNodes[retainedEdge], childNodes[removedEdge]);
        }
        int parentNode = parentNodes[removedEdge];
        int childNode = childNodes[removedEdge];
        byte symbol = symbols[removedEdge];
        assert inEdge2NextInEdge[child2FirstInEdge[childNode]] == -1;
        outEdges[symbol][parentNode] = -1;
        child2FirstInEdge[childNode] = -1;
        counts[removedEdge] = 0.0f;
        parentNodes[removedEdge] = -1;
        childNodes[removedEdge] = -1;
        symbols[removedEdge] = -1;
    }

    private void mergeHaplotypes(int retainedChild, int removedChild) {
        int hap = child2FirstHap[removedChild];
        while (hap2NextHap[hap] != -1) {
            hap = hap2NextHap[hap];
        }
        hap2NextHap[hap] = child2FirstHap[retainedChild];
        child2FirstHap[retainedChild] = child2FirstHap[removedChild];
        child2FirstHap[removedChild] = -1;
    }

    /**
     * Returns the marker.
     * @return the marker.
     */
    public Marker marker() {
        return marker;
    }

    /**
     * Returns the marker index.
     * @return the marker index.
     */
    public int markerIndex() {
        return this.markerIndex;
    }

    /**
     * Returns the number of sequences used to construct the DAG.
     * @return the number of sequences used to construct the DAG.
     */
    public int nHaps() {
        return this.nHaps;
    }

    /**
     * Returns the number of alleles.
     *
     * @return the number of alleles.
     */
    public int nAlleles() {
        return this.nAlleles;
    }

   /**
    * Returns the sum of weights for the sequences that pass
    * through the specified edge or 0.0f if the edge does not exist.
    *
    * @param edge index of the edge.
    * @return sum of weights for the sequences that pass
    * through the specified edge or 0.0f if the edge does not exist.
    *
    * @throws IndexOutOfBoundsException if
    * {@code edge<0 || edge>=this.nHaps()}
    */
    public float edgeCount(int edge) {
        return counts[edge];
    }

   /**
    * Returns the sum of weights for the sequences that pass
    * through the specified parent node or 0.0f if the parent node
    * does not exist.
    *
    * @param parentNode index of the parent node.
    * @return sum of weights for the sequences  that pass
    * through the specified parent node or 0.0f if the parent node
    * does not exist.
    *
    * @throws IndexOutOfBoundsException if
    * {@code parentNode<0 || parentNode>=this.nHaps()}    *
    */
    public float nodeCount(int parentNode) {
        float sum = 0.0f;
        for (int symbol=0; symbol<nAlleles; ++symbol) {
            if (outEdges[symbol][parentNode] >= 0) {
                sum += edgeCount(outEdges[symbol][parentNode]);
            }
        }
        return sum;
    }

    /**
     * Returns an array of parent node indices.
     * @return an array of parent node indices.
     */
    public int[] parentNodeArray() {
        int[] sortedReducedArray = DagUtils.removeValues(parentNodes, -1);
        Arrays.sort(sortedReducedArray);
        assert sortedReducedArray.length > 0;
        int cnt = 1;
        for (int j=1; j<sortedReducedArray.length; ++j) {
            if (sortedReducedArray[j] != sortedReducedArray[j-1]) {
                ++cnt;
            }
        }
        int[] parentNodeArray = new int[cnt];
        int index = 0;
        parentNodeArray[index++] = sortedReducedArray[0];
        for (int j=1; j<sortedReducedArray.length; ++j) {
            if (sortedReducedArray[j] != sortedReducedArray[j-1]) {
                parentNodeArray[index++] = sortedReducedArray[j];
            }
        }
        assert index==parentNodeArray.length;
        return parentNodeArray;
    }

   /**
     * Returns the parent node of the specified
     * edge or -1 if edge does not exist.
     *
     * @param edge index of the edge
     * @return the parent node of the specified
     * edge or -1 if edge does not exist.
     *
     * @throws IndexOutOfBoundsException if
     * {@code edge<0 || edge>=this.nHaps()}
     */
    public int parentNode(int edge) {
        return parentNodes[edge];
    }

    /**
     * Returns the child node of the specified
     * edge or -1 if the edge does not exist.
     *
     * @param edge the edge
     * @return the child node of the specified
     * edge or -1 if the edge does not exist.
     *
     * @throws IndexOutOfBoundsException if
     * {@code edge<0 || edge>=this.Haplotypes()}
     */
    public int childNode(int edge) {
        return childNodes[edge];
    }

    /**
     * Returns the outgoing edge of the specified
     * parent parent node that has the specified symbol, or returns
     * -1 if no such edge exists.
     *
     * @param parentNode the parent node.
     * @param symbol symbol labeling the outgoing edge.
     * @return the outgoing edge of the specified
     * parent parent node that has the specified symbol, or
     * -1 if no such edge exists.
     *
     * @throws IndexOutOfBoundsException if
     * {@code parentNode<0 || parentNode>=this.nHaps()}
     * or if {@code symbol<0 || symbol>=this.nAlleles()}.
     */
    public int outEdge(int parentNode, int symbol) {
        return outEdges[symbol][parentNode];
    }

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     *
     * @return a string representation of {@code this}.
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(1000);
        sb.append(Const.nl);
        sb.append("[ MergeableDagLevel: marker=");
        sb.append(markerIndex);
        sb.append(Const.nl);
        for (int j=0, n=nHaps(); j<n; ++j) {
            if (parentNodes[j] != -1) {
                sb.append("edge=");
                sb.append(j);
                sb.append(" parent=");
                sb.append(parentNodes[j]);
                sb.append(" child=");
                sb.append(childNodes[j]);
                sb.append(" symbol=");
                sb.append(symbols[j]);
                sb.append(" count=");
                sb.append(counts[j]);
                sb.append(Const.nl);
            }
        }
        sb.append("previous=");
        sb.append(prevLevel!=null);
        sb.append(" next=");
        sb.append(nextLevel!=null);
        sb.append(Const.nl);
        sb.append(" ]");
        return sb.toString();
    }

    private boolean isParentNode(int node) {
        if (prevLevel!=null) {
            return prevLevel.child2FirstInEdge[node]>=0;
        }
        else {
            for (int j=0; j<nAlleles; ++j) {
                if (outEdges[j][node] != -1) {
                    return true;
                }
            }
            return false;
        }
    }
}