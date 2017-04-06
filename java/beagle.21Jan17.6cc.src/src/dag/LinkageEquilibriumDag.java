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
import java.util.Arrays;
import vcf.GL;
import vcf.Markers;

/**
 * <p>Class {@code LinkageEquilibriumDag} represents a leveled DAG with one parent
 * node at each level.
 * </p>
 * <p>Instances of class {@code LinkageEquilibriumDag} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class LinkageEquilibriumDag implements Dag {

    private final Markers markers;
    private final float[][] alleleFreq;
    private final int maxAlleles;
    private final int sumAlleles;

    /**
     * Constructs a new {@code LinkageEquilibriumDag} instance that represents
     * markers in linkage equilibrium, with one level per marker,
     * one parent node per level, one edge per allele at each level,
     * and edge count equal to the estimated allele frequency.
     * @param gl the genotype emission probabilities which determine
     * the estimated allele frequencies
     * @param minFreq the minimum allele frequency that will be used
     * @throws IllegalArgumentException if
     * {@code minFreq <= 0.0f || minFreq >= 0.5f || Float.isNaN(minFreq) == true}
     * @throws NullPointerException if {@code gl == null}
     */
    public LinkageEquilibriumDag(GL gl, float minFreq) {
        if (minFreq <= 0.0f || minFreq >= 0.5f || Float.isNaN(minFreq)) {
            throw new IllegalArgumentException(String.valueOf(minFreq));
        }
        int nMarkers = gl.nMarkers();
        int localMaxAlleles = 0;
        this.markers = gl.markers();
        this.alleleFreq = new float[nMarkers][];
        for (int marker=0; marker<nMarkers; ++marker) {
            alleleFreq[marker] = alleleFrequencies(gl, marker, minFreq);
            if (alleleFreq[marker].length > localMaxAlleles) {
                localMaxAlleles = alleleFreq[marker].length;
            }
        }
        this.maxAlleles = localMaxAlleles;
        this.sumAlleles = gl.markers().sumAlleles();
    }

    private static float[] alleleFrequencies(GL gl, int marker,
             float minFreq) {
        int nSamples = gl.nSamples();
        int nAlleles = gl.marker(marker).nAlleles();
        float[] alleleFreq = new float[nAlleles];
        float[] scaledFreq = new float[nAlleles];
        for (int sample=0; sample<nSamples; ++sample) {
            for (int a1=0; a1<nAlleles; ++a1) {
                for (int a2=0; a2<nAlleles; ++a2) {
                    float likelihood = gl.gl(marker, sample, a1, a2);
                    scaledFreq[a1] += likelihood;
                    scaledFreq[a2] += likelihood;
                }
            }
            divideEntriesBySum(scaledFreq);
            for (int j=0; j<scaledFreq.length; ++j) {
                alleleFreq[j] += scaledFreq[j];
                scaledFreq[j] = 0.0f;
            }
        }
        divideEntriesBySum(alleleFreq);
        enforceMinFrequency(alleleFreq, minFreq);
        return alleleFreq;
    }

    private static void divideEntriesBySum(float[] fa) {
        float sum = 0.0f;
        for (float f : fa) {
            sum += f;
        }
        for (int j=0; j<fa.length; ++j) {
            fa[j] /= sum;
        }
    }

    private static void enforceMinFrequency(float[] alleleFreq, float minFreq) {
        boolean changedFreq = false;
        for (int j=0; j<alleleFreq.length; ++j) {
            if (alleleFreq[j] < minFreq) {
                alleleFreq[j] = minFreq;
                changedFreq = true;
            }
        }
        if (changedFreq) {
            divideEntriesBySum(alleleFreq);
        }
    }

    private void checkLevel(int level) {
        if (level<0 || level >= alleleFreq.length) {
            throw new IllegalArgumentException("level: " + level);
        }
    }

     private void checkEdge(int level, int edge) {
        if (edge<0 || edge>=alleleFreq[level].length) {
            throw new IndexOutOfBoundsException("edge: " + (int) edge);
        }
    }

    private void checkParentNode(int level, int node) {
        checkLevel(level);
        if (node!=0) {
            throw new IndexOutOfBoundsException("node: " + (int) node);
        }
    }

    @Override
    public int nEdges(int level) {
        return alleleFreq[level].length;
    }

    @Override
    public int nParentNodes(int level) {
        checkLevel(level);
        return 1;
    }

    @Override
    public int nChildNodes(int level) {
        checkLevel(level);
        return 1;
    }

    @Override
    public int parentNode(int level, int edge) {
        checkEdge(level, edge);
        return 0;
    }

    @Override
    public int childNode(int level, int edge) {
        checkEdge(level, edge);
        return 0;
    }

    @Override
    public int symbol(int level, int edge) {
        checkEdge(level, edge);
        return edge;
    }

    @Override
    public float edgeWeight(int level, int edge) {
        return alleleFreq[level][edge];
    }

    @Override
    public float parentWeight(int level, int parentNode) {
        checkParentNode(level, parentNode);
        return 1.0f;
    }

    @Override
    public float condEdgeProb(int level, int edge) {
        return alleleFreq[level][edge];
    }

    @Override
    public float edgeProb(int level, int edge) {
        return alleleFreq[level][edge];
    }

    @Override
    public float parentProb(int level, int node) {
        checkParentNode(level, node);
        return 1.0f;
    }

    @Override
    public int nLevels() {
        return alleleFreq.length;
    }

    @Override
    public Markers markers() {
        return markers;
    }

    @Override
    public long nNodes() {
        return (alleleFreq.length + 1);
    }

    @Override
    public long nEdges() {
        return sumAlleles;
    }

    @Override
    public int maxNodes() {
        return 1;
    }

    @Override
    public int maxEdges() {
        return maxAlleles;
    }

    @Override
    public int nOutEdges(int level, int parentNode) {
        return alleleFreq[level].length;
    }

    @Override
    public int outEdge(int level, int parentNode, int outEdge) {
        checkParentNode(level, parentNode);
        checkEdge(level, outEdge);
        return outEdge;
    }

    @Override
    public int outEdgeBySymbol(int level, int parentNode, int symbol) {
        return symbol;
    }

    @Override
    public int nInEdges(int level, int childNode) {
        checkLevel(level);
        return alleleFreq[level].length;
    }

    @Override
    public int inEdge(int level, int childNode, int inEdge) {
        checkEdge(level, inEdge);
        if (childNode!=0) {
            throw new IllegalArgumentException("childNode: " + (int) childNode);
        }
        return inEdge;
    }

    @Override
    public boolean isChildOf(int parentLevel, int parentEdge, int childEdge) {
        checkEdge(parentLevel, parentEdge);
        checkEdge(parentLevel+1, childEdge);
        return true;
    }

    @Override
    public double[] posArray() {
        double[] pos = new double[alleleFreq.length];
        for (int j=0; j<pos.length; ++j) {
            double condEdgeProb = 0.0;
            for (int a=0; a<alleleFreq[j].length; ++a) {
                condEdgeProb += alleleFreq[j][a]*alleleFreq[j][a];
            }
            if (j==0) {
                pos[j] = -Math.log10(condEdgeProb);
            }
            else {
                pos[j] = pos[j-1] - Math.log10(condEdgeProb);
            }
        }
        return pos;
    }

    @Override
    public String toString(int start, int end) {
        if (start<0 || start>end || end>=alleleFreq.length) {
            String s = "start=" + start + " end=" + end;
            throw new IllegalArgumentException(s);
        }
        StringBuilder sb = new StringBuilder((end-start) * 20);
        for (int level=start; level<end; ++level) {
            sb.append(Arrays.toString(alleleFreq[level]));
            sb.append(Const.nl);
        }
        return sb.toString();
    }
}
