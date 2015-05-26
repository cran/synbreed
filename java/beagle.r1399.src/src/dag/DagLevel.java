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

import vcf.Marker;

/**
 * <p>Interface {@code DagLevel} represents a level of a leveled directed
 * acyclic graph (DAG).
 * </p>
 * All instances of {@code DagLevel} are required to be immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface DagLevel {

    /**
     * Returns the marker corresponding to this level of the DAG.
     *
     * @return the marker corresponding to this level of the DAG.
     */
    public Marker marker();

    /**
     * Returns the number of edges at this level of the DAG.
     *
     * @return the number of edges at this level of the DAG.
     */
    public int nEdges();

    /**
     * Returns the number of parent nodes at this level of the DAG.
     *
     * @return the number of parent nodes at this level of the DAG.
     */
    public char nParentNodes();

    /**
     * Returns the number of child nodes at this level of the DAG.
     *
     * @return the number of child nodes at this level of the DAG.
     */
    public int nChildNodes();

    /**
     * Returns the index of the parent node of the specified edge
     * at this level of the DAG.
     *
     * @param edge an edge index.
     * @return the index of the parent node of the specified edge
     * at this level of the DAG.
     *
     * @throws IndexOutOfBoundsException if
     * {@code edge<0 || edge>=this.nEdges()}
     */
    public char parentNode(int edge);

    /**
     * Returns the index of the child node of the specified edge
     * at this level of the DAG.
     *
     * @param edge an edge index.
     * @return the index of the child node of the specified edge
     * at this level of the DAG.
     *
     * @throws IndexOutOfBoundsException if
     * {@code edge<0 || edge>=this.nEdges()}
     */
    public char childNode(int edge);

    /**
     * Returns the symbol labeling the specified edge at this level
     * of the DAG.
     *
     * @param edge an edge index.
     * @return the symbol labeling the specified edge at this level
     * of the DAG.
     *
     * @throws IndexOutOfBoundsException if
     * {@code edge<0 || edge>=this.nEdges()}
     */
    public byte symbol(int edge);

    /**
     * Returns the sum of weights for the sequences that pass
     * through the specified edge at this level of the DAG.
     *
     * @param edge an edge index.
     * @return the sum of weights for the sequences that pass
     * through the specified edge at this level of the DAG.
     *
     * @throws IndexOutOfBoundsException if
     * {@code edge<0 || edge>=this.nEdges()}
     */
    public float edgeCnt(int edge);

    /**
     * Returns the sum of weights for the sequences that pass
     * through the specified node at this level of the DAG.
     *
     * @param parentNode a parent node index.
     * @return the sum of weights for the sequences that pass
     * through the specified node at this level of the DAG.
     *
     * @throws IndexOutOfBoundsException if
     * {@code parentNode<0 || parentNode>=this.nParentNodes()}
     */
    public float nodeCnt(int parentNode);

    /**
     * Returns the ratio of the sum of the weights of the sequences that pass
     * through the specified edge at this level of the DAG and
     * the sum of the weights of the sequences that pass through the parent
     * node of the specified edge.
     *
     * @param edge an edge index.
     * @return the ratio of the sum of the weights of the sequences that pass
     * through the specified edge at this level of the DAG and
     * the sum of the weights of the sequences that pass through the parent
     * node of the specified edge.
     *
     * @throws IndexOutOfBoundsException if
     * {@code edge<0 || edge>=this.nEdges()}
     */
    public float condEdgeProb(int edge);

    /**
     * Returns the ratio of the sum of the weights of the sequences that pass
     * through the specified edge at this level of the DAG and
     * the sum of the weights of the sequences that pass through
     * any edge at this level of the DAG.
     *
     * @param edge an edge index.
     * @return the ratio of the sum of the weights of the sequences that pass
     * through the specified edge at this level of the DAG and
     * the sum of the weights of the sequences that pass through
     * any edge at this level of the DAG.
     *
     * @throws IndexOutOfBoundsException if
     * {@code edge<0 || edge>=this.nEdges()}
     */
    public float edgeProb(int edge);

    /**
     * Returns the ratio of the sum of the weights of the sequences that pass
     * through the specified parent node at this level of the DAG and the sum of
     * the weights of the sequences that pass through any parent node at this
     * level of the DAG.
     *
     * @param parentNode a parent node index.
     * @return the ratio of the sum of the weights of the sequences that pass
     * through the specified parent node at this level of the DAG and the sum of
     * the weights of the sequences that pass through any parent node at this
     * level of the DAG.
     *
     * @throws IndexOutOfBoundsException if
     * {@code parentNode<0 || parentNode>=this.nParentNodes()}
     */
    public float nodeProb(int parentNode);

    /**
     * Returns the number of outgoing edges of the specified parent node
     * at this level of the DAG.
     *
     * @param parentNode a parent node index.
     * @return the number of outgoing edges of the specified parent node
     * at this level of the DAG.
     *
     * @throws IndexOutOfBoundsException if
     * {@code parentNode<0 || parentNode>=nParentNodes()}
     */
    public char nOutEdges(int parentNode);

    /**
     * Returns the index of the specified edge at this level of the DAG.
     *
     * @param parentNode a parent node index.
     * @param outEdge the index of the outgoing edge of the specified
     * parent node.
     * @return the index of the specified edge at this level of the DAG.
     *
     * @throws IndexOutOfBoundsException if
     * {@code parentNode<0 || parentNode>=nParentNodes()}
     * @throws IndexOutOfBoundsException if
     * {@code outEdge<0 || outEdge>=this.nOutEdges(parentNode)}
     */
    public char outEdge(int parentNode, int outEdge);

    /**
     * Returns the index of the specified edge at this level of the
     * DAG or {@code Character.MAX_VALUE} if no such edge exists.
     *
     * @param parentNode a parent node index.
     * @param symbol a symbol labeling an outgoing edge of the specified
     * parent node.
     * @return the index of the specified edge at this level of the
     * DAG or {@code Character.MAX_VALUE} if no such edge exists.
     *
     * @throws IndexOutOfBoundsException if
     * {@code parentNode<0 || parentNode>=nParentNodes()}
     * @throws IndexOutOfBoundsException if
     * {@code symbol<0 || symbol>=this.nSumbols()}
     */
    public char outEdgeBySymbol(int parentNode, byte symbol);

    /**
     * Returns the number of ingoing edges for the specified child node
     * at this level of the DAG.
     *
     * @param childNode a child node index.
     * @return the number of ingoing edges for the specified child node
     * at this level of the DAG.
     *
     * @throws IndexOutOfBoundsException if
     * {@code childNode<0 || childNode>=this.nChildNodes()}
     */
    public char nInEdges(int childNode);

    /**
     * Returns the index of the specified edge at this level of the DAG.
     *
     * @param childNode index of the child node.
     * @param inEdge index of an ingoing edge of the specified child node.
     * @return the index of the specified edge at this level of the DAG.
     *
     * @throws IndexOutOfBoundsException if
     * {@code childNode<0 || childNode>=this.nChildNodes()}
     * @throws IndexOutOfBoundsException if
     * {@code inEdge<0 || inEdge>=this.nInEdges(childNode)}.
     */
    public char inEdge(int childNode, int inEdge);

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     *
     * @return a string representation of {@code this}.
     */
    @Override
    public String toString();

}
