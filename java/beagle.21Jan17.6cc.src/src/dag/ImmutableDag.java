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
 * <p>Class {@code ImmutableDag} represents a leveled Directed Acyclic Graph
 * (DAG).
 * </p>
 * <p>Instances of class {@code ImmutableDag} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class ImmutableDag implements Dag {

    private final Markers markers;
    private final long nNodes;         // total number of nodes
    private final long nEdges;         // total number of edges
    private final int maxNodes;        // maximum number of nodes on any level
    private final int maxEdges;        // maximum number of edges on any level

    private final DagLevel[] dagLevels;
    private final double[] posArray;

    /**
     * Constructs a new {@code ImmutableDag} instance from the specified data.
     * @param markers the markers
     * @param levels the levels of the leveled DAG
     * @throws IllegalArgumentException if {@code levels.length == 0}
     * @throws IllegalArgumentException if {@code levels[0].nParentNodes() != 1}
     * @throws IllegalArgumentException if
     * {@code levels[j-1].nChildNodes() != levels[j].nParentNodes()} for any
     * {@code j} satisfying {@code 0 < j && j < levels.length}
     * @throws NullPointerException if
     * {@code (markers == null || levels==null)}, or if
     * {@code levels[j] == null} for any
     * {@code j} satisfying {@code 0 <= j && j < levels.length}
     */
    public ImmutableDag(Markers markers, DagLevel[] levels) {
        if (levels.length==0) {
            throw new IllegalArgumentException("levels.length==0");
        }
        if (levels[0].nParentNodes()!=1) {
            throw new IllegalArgumentException("levels[0].nParentNodes()!=1");
        }
        int cumulativeNodeCnt = 1;
        int cumulativeEdgeCnt = 0;
        int maxNodesPerLevel = 0;
        int maxEdgesPerLevel = 0;
        double[] pos = new double[levels.length];
        for (int j=0; j<levels.length; ++j) {
            if (j>0 && levels[j-1].nChildNodes()!=levels[j].nParentNodes()) {
                throw new IllegalArgumentException("inconsistent levels");
            }
            cumulativeNodeCnt += levels[j].nChildNodes();
            cumulativeEdgeCnt += levels[j].nEdges();
            if (levels[j].nChildNodes() > maxNodesPerLevel) {
                maxNodesPerLevel = levels[j].nChildNodes();
            }
            if (levels[j].nEdges() > maxEdgesPerLevel) {
                maxEdgesPerLevel = levels[j].nEdges();
            }
            double d = minusLog10CondEdgeProb(levels[j]);
            pos[j] = (j==0) ? d : (pos[j-1] + d);
        }
        this.dagLevels = levels.clone();
        this.posArray = pos;
        this.markers = markers;
        this.nNodes = cumulativeNodeCnt;
        this.nEdges = cumulativeEdgeCnt;
        this.maxEdges = maxEdgesPerLevel;
        this.maxNodes = maxNodesPerLevel;
    }

    private static double minusLog10CondEdgeProb(DagLevel level) {
        float meanScore = 0.0f;
        for (int e=0, n=level.nEdges(); e<n; ++e) {
            meanScore += level.edgeProb(e)*level.condEdgeProb(e);
        }
        double d = -Math.log10(meanScore);
        return (d<0) ? 0.0 : d;
    }

    @Override
    public int nEdges(int level) {
        return dagLevels[level].nEdges();
    }

    @Override
    public int nParentNodes(int level) {
        return dagLevels[level].nParentNodes();
    }

    @Override
    public int nChildNodes(int level) {
        return dagLevels[level].nChildNodes();
    }

    @Override
    public int parentNode(int level, int edge) {
        return dagLevels[level].parentNode(edge);
    }

    @Override
    public int childNode(int level, int edge) {
        return dagLevels[level].childNode(edge);
    }

    @Override
    public int symbol(int level, int edge) {
        return dagLevels[level].symbol(edge);
    }

    @Override
    public float edgeWeight(int level, int edge) {
        return dagLevels[level].edgeWeight(edge);
    }


    @Override
    public float parentWeight(int level, int parentNode) {
        return dagLevels[level].parentWeight(parentNode);
    }

    @Override
    public float condEdgeProb(int level, int edge) {
        return dagLevels[level].condEdgeProb(edge);
    }


    @Override
    public float edgeProb(int level, int edge) {
        return dagLevels[level].edgeProb(edge);
    }

    @Override
    public float parentProb(int level, int node) {
        return dagLevels[level].parentProb(node);
    }

    @Override
    public int nLevels() {
        return dagLevels.length;
    }

    @Override
    public Markers markers() {
        return markers;
    }

    @Override
    public long nNodes() {
        return nNodes;
    }

    @Override
    public long nEdges() {
        return nEdges;
    }

    @Override
    public int maxNodes() {
        return maxNodes;
    }

    @Override
    public int maxEdges() {
        return maxEdges;
    }

    @Override
    public int nOutEdges(int level, int parentNode) {
        return dagLevels[level].nOutEdges(parentNode);
    }

    @Override
    public int outEdge(int level, int parentNode, int outEdge) {
        return dagLevels[level].outEdge(parentNode, outEdge);
    }

    @Override
    public int outEdgeBySymbol(int level, int parentNode, int symbol) {
        return dagLevels[level].outEdgeBySymbol(parentNode, symbol);
    }

    @Override
    public int nInEdges(int level, int childNode) {
        return dagLevels[level].nInEdges(childNode);
    }

    @Override
    public int inEdge(int level, int childNode, int inEdge) {
        return dagLevels[level].inEdge(childNode, inEdge);
    }

    @Override
    public boolean isChildOf(int parentLevel, int parentEdge, int childEdge) {
        int nodeA = dagLevels[parentLevel+1].parentNode(childEdge);
        int nodeB = dagLevels[parentLevel].childNode(parentEdge);
        return nodeA==nodeB;
    }

    @Override
    public double[] posArray() {
        return posArray.clone();
    }

    @Override
    public String toString(int startLevel, int endLevel) {
        String nl = System.getProperty("line.separator");
        StringBuilder sb = new StringBuilder(1000);
        for (int j=startLevel; j<endLevel; ++j) {
            sb.append("level=");
            sb.append(j);
            sb.append(": ");
            sb.append(dagLevels[j]);
            sb.append(nl);
        }
        return sb.toString();
    }

    @Override
    public String toString() {
        String nl = System.getProperty("line.separator");
        StringBuilder sb = new StringBuilder(2000);
        sb.append("[Dag: nMarkers=");
        sb.append(dagLevels.length);
        sb.append("  nodes=");
        sb.append(nNodes);
        sb.append("  edges=");
        sb.append(nEdges);
        sb.append("  maxNodes=");
        sb.append((int) maxNodes);
        sb.append("  maxEdges=");
        sb.append((int) maxEdges);
        sb.append(nl);
        for (int j=0; j<dagLevels.length; ++j) {
            sb.append(nl);
            sb.append("level=");
            sb.append(j);
            sb.append(": ");
            sb.append(dagLevels[j]);
            sb.append(nl);
        }
        return sb.toString();
    }
}
