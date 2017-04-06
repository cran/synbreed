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

import vcf.Markers;

/**
 * <p>Interface {@code Dag} represents a leveled directed acyclic graph (DAG).
 * </p>
 * All instances of {@code DAG} are required to be immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface Dag {

    /**
     * Returns the number of edges at the specified level of the DAG.
     *
     * @param level a level of the DAG
     * @return the number of edges at the specified level of the DAG.
     * @throws IndexOutOfBoundsException if
     * {@code level < 0 || level >= this.nMarkers()}
     */
    public int nEdges(int level);

    /**
     * Returns the number of parent nodes at the specified level of the DAG.
     *
     * @param level a level of the DAG
     * @return the number of parent nodes at the specified level of the DAG
     *
     * @throws IndexOutOfBoundsException if
     * {@code level < 0 || level >= this.nMarkers()}
     */
    public int nParentNodes(int level);

    /**
     * Returns the number of child nodes at the specified level of the DAG.
     *
     * @param level a level of the DAG
     * @return the number of child nodes at the specified level of the DAG
     *
     * @throws IndexOutOfBoundsException if
     * {@code level < 0 || level >= this.nMarkers()}
     */
    public int nChildNodes(int level);

    /**
     * Returns the index of the specified parent node in the DAG.
     *
     * @param level a level of the DAG.
     * @param edge the index of an edge at the specified level of the DAG
     * @return the index of the specified parent node in the DAG
     *
     * @throws IndexOutOfBoundsException if
     * {@code level < 0 || level >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code edge < 0 || edge >= this.nEdges(level)}
     */
    public int parentNode(int level, int edge);

    /**
     * Returns the index of the specified child node in the DAG.
     *
     * @param level a level of the DAG
     * @param edge the index of an edge at the specified level of the DAG
     * @return the index of the specified child node in the DAG
     *
     * @throws IndexOutOfBoundsException if
     * {@code level < 0 || level >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code edge < 0 || edge >= this.nEdges(level)}
     */
    public int childNode(int level, int edge);

    /**
     * Returns the symbol labeling the specified edge of the DAG.
     *
     * @param level a level of the DAG
     * @param edge the index of an edge at the specified level of the DAG
     * @return the symbol labeling the specified edge of the DAG
     *
     * @throws IndexOutOfBoundsException if
     * {@code level < 0 || level >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code edge < 0 || edge >= this.nEdges(level)}
     */
    public int symbol(int level, int edge);

    /**
     * Returns the sum of the weights of the sequences that pass
     * through the specified edge of the DAG.
     *
     * @param level a level of the DAG
     * @param edge the index of an edge at the specified level of the DAG
     * @return the sum of the weights of the sequences that pass
     * through the specified edge of the DAG
     *
     * @throws IndexOutOfBoundsException if
     * {@code level < 0 || level >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code edge < 0 || edge >= this.nEdges(level)}
     */
    public float edgeWeight(int level, int edge);

    /**
     * Returns the sum of the weights of the sequences that pass
     * through the specified node of the DAG.
     *
     * @param level a level of the DAG
     * @param parentNode the index of a parent node at the specified level
     * of the DAG
     * @return the sum of the weights of the sequences that pass
     * through the specified node of the DAG
     *
     * @throws IndexOutOfBoundsException if
     * {@code level < 0 || level >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code parentNode < 0 || node >= this.nParentNodes(level)}
     */
    public float parentWeight(int level, int parentNode);

    /**
     * Returns the ratio of the sum of the weights of the sequences that pass
     * through the specified edge of the DAG and
     * the sum of the weights of the sequences that pass through the parent
     * node of the specified edge of the DAG.
     *
     * @param level a level of the DAG
     * @param edge the index of an edge at the specified level of the DAG
     * @return the ratio of the sum of the weights of the sequences that pass
     * through the specified edge of the DAG and
     * the sum of the weights of the sequences that pass through the parent
     * node of the specified edge of the DAG
     *
     * @throws IndexOutOfBoundsException if
     * {@code level < 0 || level >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code edge < 0 || edge >= this.nEdges(level)}
     */
    public float condEdgeProb(int level, int edge);

    /**
     * Returns the ratio of the sum of the weights of the sequences that pass
     * through the specified edge of the DAG and the sum of the weights of all
     * sequences.
     *
     * @param level a level of the DAG
     * @param edge the index of an edge at the specified level of the DAG
     * @return the ratio of the sum of the weights of the sequences that pass
     * through the specified edge of the DAG and the sum of the weights of all
     * sequences
     *
     * @throws IndexOutOfBoundsException if
     * {@code level < 0 || level >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code edge < 0 || edge >= this.nEdges(level)}
     */
    public float edgeProb(int level, int edge);

    /**
     * Returns the ratio of the sum of the weights of the sequences that pass
     * through the specified parent node of the DAG and the sum of the weights
     * of all sequences.
     *
     * @param level a level of the DAG
     * @param parentNode the index of a parent node at the specified level
     * of the DAG
     * @return the ratio of the sum of the weights of the sequences that pass
     * through the specified parent node of the DAG and the sum of the weights
     * of all sequences
     *
     * @throws IndexOutOfBoundsException if
     * {@code level < 0 || level >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code parentNode < 0 || parentNode >= this.nParentNodes(level)}
     */
    public float parentProb(int level, int parentNode);

    /**
     * Returns the number of markers.
     *
     * @return the number of markers
     */
    public int nLevels();

    /**
     * Returns the markers represented by this DAG.
     * @return the markers represented by this DAG
     */
    public Markers markers();

    /**
     * Returns the number of nodes in the DAG.
     *
     * @return the number of node in the DAG
     */
    public long nNodes();

    /**
     * Returns the number of edges in the DAG.
     *
     * @return the number of edges in the DAG
     */
    public long nEdges();

    /**
     * Returns the maximum number of parent nodes at any level of the DAG.
     *
     * @return the maximum number of parent nodes at any level of the DAG
     */
    public int maxNodes();

    /**
     * Returns the maximum number of edges at any level of the DAG.
     *
     * @return the maximum number of edges at any level of the DAG
     */
    public int maxEdges();

    /**
     * Returns the number of outgoing edges for the specified node of the DAG.
     *
     * @param level a level of the DAG
     * @param parentNode the index of a parent node at the specified
     * level of the DAG
     * @return the number of outgoing edges for the specified node of the DAG
     *
     * @throws IndexOutOfBoundsException if
     * {@code level < 0 || level >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code parentNode < 0 || parentNode >= this.nParentNodes(level)}
     */
    public int nOutEdges(int level, int parentNode);

    /**
     * Returns the index of the specified edge in the DAG.
     *
     * @param level a level of the DAG
     * @param parentNode the index of a parent node at the specified
     * level of the DAG
     * @param outEdge the index of an outgoing edge of the specified
     * parent node
     * @return the index of the specified edge in the DAG
     *
     * @throws IndexOutOfBoundsException if
     * {@code level < 0 || level >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code parentNode < 0 || parentNode >= nParentNodes(level)}
     * @throws IndexOutOfBoundsException if
     * {@code outEdge < 0 || outEdge >= this.nOutEdges(level, parentNode)}
     */
    public int outEdge(int level, int parentNode, int outEdge);

    /**
     * Returns the index of the specified edge at the specified level of the
     * DAG or {@code -1} if no such edge exists.
     *
     * @param level a level of the DAG
     * @param parentNode the index of a parent node at the specified
     * level of the DAG
     * @param symbol a symbol labeling an outgoing edge of the specified
     * parent node of the DAG
     * @return the index of the specified edge at the specified level of the
     * DAG or {@code -1} if no such edge exists
     *
     * @throws IndexOutOfBoundsException if
     * {@code level < 0 || level >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code parentNode < 0 || parentNode >= this.nParentNodes(level)}
     * @throws IndexOutOfBoundsException if
     * {@code symbol < 0 || symbol >= this.marker(level).nAlleles()}}
     */
    public int outEdgeBySymbol(int level, int parentNode, int symbol);

    /**
     * Returns the number of ingoing edges for the specified node of the DAG.
     *
     * @param level a level of the DAG
     * @param childNode the index of a child node at the specified
     * level of the DAG
     * @return the number of ingoing edges for the specified node of the DAG
     *
     * @throws IndexOutOfBoundsException if
     * {@code level < 0 || level >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code childNode < 0 || childNode >= this.nChildNodes(level)}
     */
    public int nInEdges(int level, int childNode);

    /**
     * Returns the index of the specified edge in the DAG.
     *
     * @param level a level of the DAG
     * @param childNode the index of a child node at the specified
     * level of the DAG
     * @param inEdge the index of an ingoing edge of the specified
     * child node in the DAG
     * @return the index of the specified edge in the DAG
     *
     * @throws IndexOutOfBoundsException if
     * {@code level < 0 || level >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code childNode < 0 || childNode >= nChildNodes(level)}
     * @throws IndexOutOfBoundsException if
     * {@code inEdge < 0 || inEdge >= this.nInEdges(level, childNode)}
     */
    public int inEdge(int level, int childNode, int inEdge);

    /**
     * Returns {@code true} if the child node of the specified parent
     * edge equals the parent node of the specified child edge and
     * returns {@code false} otherwise.
     *
     * @param parentLevel a level of the DAG
     * @param parentEdge the index of an edge at the specified level
     * of the DAG
     * @param childEdge the index of an edge at level {@code (parentLevel + 1)}
     * of the DAG
     * @return {@code true} if the child node of the specified parent
     * edge equals the parent node of the specified child edge
     *
     * @throws IndexOutOfBoundsException if
     * {@code parentLevel < 0 || parentLevel >= (this.nMarkers() - 1)}
     * @throws IndexOutOfBoundsException if
     * {@code parentEdge < 0 || parentEdge >= this.nEdges(level)}
     * @throws IndexOutOfBoundsException if
     * {@code childEdge < 0 || childEdge >= this.nEdges(level + 1)}
     */
    public boolean isChildOf(int parentLevel, int parentEdge, int childEdge);

    /**
     * Returns an array of length {@code this.nMarkers()} whose {@code j}-th
     * element is the distance from the root node to
     * the child node at level {@code j} of the DAG.
     * The distance from parent node to child node at level {@code lev}
     * equals {@code -Math.log10(P)} where {@code P} is the weighted conditional
     * edge probability at level {@code lev}, when each edge {@code e} is
     * weighted by {@code this.counts(lev, e)}.
     *
     * @return an array of length {@code this.nMarkers()} whose {@code j}-th
     * element is the distance from the root node to
     * the child node at level {@code j} of the DAG
     */
    public double[] posArray();

    /**
     * Returns a description of the specified levels of the DAG.  The
     * exact details of the description are unspecified and subject to change.
     *
     * @param start the first level (inclusive)
     * @param end the last level (exclusive)
     * @return a description of the specified levels of the DAG
     *
     * @throws IllegalArgumentException if
     * {@code start < 0 || start > end || end >= this.nMarkers()}
     */
    public String toString(int start, int end);

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecifed and subject to change.
     *
     * @return a string representation of {@code this}
     */
    @Override
    String toString();
}
