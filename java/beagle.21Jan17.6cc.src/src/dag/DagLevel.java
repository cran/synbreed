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

/**
 * <p>Interface {@code DagLevel} represents a level of a leveled directed
 * acyclic graph (DAG).
 * </p>
 * <p>All instances of {@code DagLevel} are required to be immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface DagLevel {

    /**
     * Returns the number of edges at this level of the DAG.
     *
     * @return the number of edges at this level of the DAG
     */
    public int nEdges();

    /**
     * Returns the number of parent nodes at this level of the DAG.
     *
     * @return the number of parent nodes at this level of the DAG
     */
    public int nParentNodes();

    /**
     * Returns the number of child nodes at this level of the DAG.
     *
     * @return the number of child nodes at this level of the DAG
     */
    public int nChildNodes();

    /**
     * Returns the index of the parent node of the specified edge
     * at this level of the DAG.
     *
     * @param edge an edge index
     * @return the index of the parent node of the specified edge
     * at this level of the DAG
     *
     * @throws IndexOutOfBoundsException if
     * {@code edge < 0 || edge >= this.nEdges()}
     */
    public int parentNode(int edge);

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
    public int childNode(int edge);

    /**
     * Returns the symbol labeling the specified edge at this level
     * of the DAG.
     *
     * @param edge an edge index
     * @return the symbol labeling the specified edge at this level
     * of the DAG
     *
     * @throws IndexOutOfBoundsException if
     * {@code edge < 0 || edge >= this.nEdges()}
     */
    public int symbol(int edge);

    /**
     * Returns the sum of weights for the sequences that pass
     * through the specified edge at this level of the DAG.
     *
     * @param edge an edge index
     * @return the sum of weights for the sequences that pass
     * through the specified edge at this level of the DAG
     *
     * @throws IndexOutOfBoundsException if
     * {@code edge < 0 || edge >= this.nEdges()}
     */
    public float edgeWeight(int edge);

    /**
     * Returns the sum of weights for the sequences that pass
     * through the specified node at this level of the DAG.
     *
     * @param parentNode a parent node index
     * @return the sum of weights for the sequences that pass
     * through the specified node at this level of the DAG
     *
     * @throws IndexOutOfBoundsException if
     * {@code parentNode < 0 || parentNode >= this.nParentNodes()}
     */
    public float parentWeight(int parentNode);

    /**
     * Returns the conditional edge probability, which is defined to be
     * the ratio of the sum of the weights of the sequences that pass
     * through the specified edge at this level of the DAG and
     * the sum of the weights of the sequences that pass through the parent
     * node of the specified edge.
     *
     * @param edge an edge index
     * @return the conditional edge probability
     *
     * @throws IndexOutOfBoundsException if
     * {@code edge < 0 || edge >= this.nEdges()}
     */
    public float condEdgeProb(int edge);

    /**
     * Returns the edge probability, which is defined to be the ratio of the
     * sum of the weights of the sequences that pass through the specified
     * edge at this level of the DAG and the sum of the weights of the
     * sequences that pass through any edge at this level of the DAG.
     *
     * @param edge an edge index
     * @return the edge probability
     *
     * @throws IndexOutOfBoundsException if
     * {@code edge < 0 || edge >= this.nEdges()}
     */
    public float edgeProb(int edge);

    /**
     * Returns the parent node probability, which is defined to be the
     * ratio of the sum of the weights of the sequences that pass through
     * the specified parent node at this level of the DAG and the sum of
     * the weights of the sequences that pass through any parent node at this
     * level of the DAG.
     *
     * @param parentNode a parent node index
     * @return the parent node probability
     *
     * @throws IndexOutOfBoundsException if
     * {@code parentNode < 0 || parentNode >= this.nParentNodes()}
     */
    public float parentProb(int parentNode);

    /**
     * Returns the number of outgoing edges of the specified parent node
     * at this level of the DAG.
     *
     * @param parentNode a parent node index
     * @return the number of outgoing edges of the specified parent node
     * at this level of the DAG
     *
     * @throws IndexOutOfBoundsException if
     * {@code parentNode < 0 || parentNode >= nParentNodes()}
     */
    public int nOutEdges(int parentNode);

    /**
     * Returns the index of the specified edge at this level of the DAG.
     *
     * @param parentNode a parent node index
     * @param outEdge the index of the outgoing edge of the specified
     * parent node
     * @return the index of the specified edge at this level of the DAG
     *
     * @throws IndexOutOfBoundsException if
     * {@code parentNode < 0 || parentNode >= nParentNodes()}
     * @throws IndexOutOfBoundsException if
     * {@code outEdge < 0 || outEdge >= this.nOutEdges(parentNode)}
     */
    public int outEdge(int parentNode, int outEdge);

    /**
     * Returns the index of the specified edge at this level of the
     * DAG or {@code -1} if no such edge exists.
     *
     * @param parentNode a parent node index
     * @param symbol a symbol labeling an outgoing edge of the specified
     * parent node
     * @return  the index of the specified edge at this level of the
     * DAG or {@code -1} if no such edge exists
     *
     * @throws IndexOutOfBoundsException if
     * {@code parentNode < 0 || parentNode >= nParentNodes()}
     */
    public int outEdgeBySymbol(int parentNode, int symbol);

    /**
     * Returns the number of ingoing edges for the specified child node
     * at this level of the DAG.
     *
     * @param childNode a child node index
     * @return the number of ingoing edges for the specified child node
     * at this level of the DAG
     *
     * @throws IndexOutOfBoundsException if
     * {@code childNode < 0 || childNode >= this.nChildNodes()}
     */
    public int nInEdges(int childNode);

    /**
     * Returns the index of the specified edge at this level of the DAG.
     *
     * @param childNode index of the child node
     * @param inEdge index of an ingoing edge of the specified child node
     * @return the index of the specified edge at this level of the DAG
     *
     * @throws IndexOutOfBoundsException if
     * {@code childNode < 0 || childNode >= this.nChildNodes()}
     * @throws IndexOutOfBoundsException if
     * {@code inEdge < 0 || inEdge >= this.nInEdges(childNode)}
     */
    public int inEdge(int childNode, int inEdge);

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     *
     * @return a string representation of {@code this}
     */
    @Override
    public String toString();

}
