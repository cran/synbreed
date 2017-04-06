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

import dag.Dag;
import java.util.Arrays;
import vcf.BasicGL;
import vcf.GL;

/**
 * <p>Class {@code DuoBaumLevel} computes forward and backward Baum
 * values at a level of a hidden Markov model (HMM) whose states are
 * ordered edge trios of a leveled directed acyclic graph (DAG).
 * </p>
 * <p>Instances of class {@code SingleBaumLevel} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class DuoBaumLevel {

    private static final int INITIAL_CAPACITY = 400;
    private static final float MIN_VALUE = 100*Float.MIN_VALUE;
    private final Dag dag;
    private final GL gl;

    private int marker = -1;
    private int sampleA = -1;
    private int sampleB = -1;
    private int size = 0;

    private int capacity = INITIAL_CAPACITY;
    private int[] edgesAB1 = new int[INITIAL_CAPACITY];
    private int[] edgesA2 = new int[INITIAL_CAPACITY];
    private int[] edgesB2 = new int[INITIAL_CAPACITY];
    private float[] fwdValues = new float[INITIAL_CAPACITY];
    private float[] bwdValues = new float[INITIAL_CAPACITY];
    private float fwdValueSum = 0f;
    private float bwdValueSum = 0f;

    private int nGenotypes = 0;
    private float[] gtProbsA = new float[3];
    private float[] gtProbsB = new float[3];

    /**
     * Constructs a new {@code DuoBaumLevel} instance from the specified data.
     * @param dag the directed acyclic graph that the determines transition
     * probabilities
     * @param gl the emission probabilities
     * @throws IllegalArgumentException if
     * {@code dag.markers().equals(gl.markers()) == false}
     * @throws NullPointerException if {@code dag == null || gl == null}
     */
    public DuoBaumLevel(Dag dag, GL gl) {
        if (dag.markers().equals(gl.markers())==false) {
            throw new IllegalArgumentException("marker inconsistency");
        }
        this.dag = dag;
        this.gl = gl;
    }

    /**
     * Sets the Baum forward algorithm values for this level of the HMM
     * and records the child node trio values in the specified
     * {@code nodes} parameter. When the method call returns, the {@code nodes}
     * parameter will be reset to the child node trio values for this level of
     * the HMM.
     *
     * @param nodes child node trio values at the previous level of the HMM
     * @param marker the level of the HMM at which the Baum forward algorithm
     * probabilities will be computed
     * @param sampleA the parent sample index
     * @param sampleB the offspring sample index
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.dag().nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code sampleA < 0 || sampleA >= this.gl().nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code sampleB < 0 || sampleB >= this.gl().nSamples()}
     * @throws IndexOutOfBoundsException if any node in any node trio with
     * non-zero value is not a valid parent node at the specified level of the
     * HMM
     * @throws NullPointerException if {@code nodes == null}
     */
    public void setForwardValues(DuoNodes nodes, int marker, int sampleA,
            int sampleB) {
        this.marker = marker;
        this.sampleA = sampleA;
        this.sampleB = sampleB;
        this.nGenotypes = gl.marker(marker).nGenotypes();
        this.size = 0;
        this.fwdValueSum = 0f;
        this.bwdValueSum = 0f;
        initializeGtProbs(); // initialized here due to gtProbs() contract
        setStates(nodes);
        setChildNodes(nodes);
    }

    private void initializeGtProbs() {
        if (gtProbsA.length < nGenotypes) {
            int newLength = Math.max(nGenotypes, (3*gtProbsA.length/2 + 1));
            gtProbsA = new float[newLength];
            gtProbsB = new float[newLength];
        }
        else {
            for (int j=0; j<nGenotypes; ++j) {
                gtProbsA[j] = 0f;
                gtProbsB[j] = 0f;
            }
        }
    }

    private void setStates(DuoNodes nodes) {
        float valueSum = 0f;
        for (int j=0, n=nodes.size(); j<n; ++j) {
            int nodeAB1 = nodes.enumNodeAB1(j);
            int nodeA2 = nodes.enumNodeA2(j);
            int nodeB2 = nodes.enumNodeB2(j);
            float nodeValue = nodes.enumValue(j);
            for (int ab1=0, nAB1=dag.nOutEdges(marker, nodeAB1); ab1<nAB1; ++ab1) {
                int edgeAB1 = dag.outEdge(marker, nodeAB1, ab1);
                int symbolAB1 = dag.symbol(marker, edgeAB1);
                for (int a2=0, nA2=dag.nOutEdges(marker, nodeA2); a2<nA2; ++a2) {
                    int edgeA2 = dag.outEdge(marker, nodeA2, a2);
                    int symbolA2 = dag.symbol(marker, edgeA2);
                    float epA = gl.gl(marker, sampleA, symbolAB1, symbolA2);
                    if (epA > 0.0) {
                        for (int b2=0, nB2=dag.nOutEdges(marker, nodeB2); b2<nB2; ++b2) {
                            int edgeB2 = dag.outEdge(marker, nodeB2, b2);
                            int symbolB2 = dag.symbol(marker, edgeB2);
                            float epB = gl.gl(marker, sampleB, symbolAB1, symbolB2);
                            if (epB > 0.0) {
                                if (size == capacity) {
                                    ensureCapacity(size+1);
                                }
                                float tpAB1 = dag.condEdgeProb(marker, edgeAB1);
                                float tpA2 = dag.condEdgeProb(marker, edgeA2);
                                float tpB2 = dag.condEdgeProb(marker, edgeB2);
                                float fwdValue = (epA * epB) * nodeValue
                                        * (tpAB1 * tpA2 * tpB2);
                                if (fwdValue<MIN_VALUE && nodeValue > 0.0) {
                                    fwdValue = MIN_VALUE;
                                }
                                edgesAB1[size] = edgeAB1;
                                edgesA2[size] = edgeA2;
                                edgesB2[size] = edgeB2;
                                fwdValues[size++] = fwdValue;
                                valueSum += fwdValue;
                            }
                        }
                    }
                }
            }
        }
        assert valueSum>0.0 ^ size==0;
        for (int k=0; k<size; ++k) {
            this.fwdValues[k] /= valueSum;
        }
        fwdValueSum = valueSum;
    }

    /**
     * Stores the Baum forward algorithm child node trio values for this
     * level of the HMM in the specified {@code DuoNodes} object.
     *
     * @param nodes the node trio values that will be set
     *
     * @throws NullPointerException if {@code nodes == null}
     */
    public void setChildNodes(DuoNodes nodes) {
        nodes.clear();
        for (int k=0; k<size; ++k) {
            int nodeAB1 = dag.childNode(marker, edgesAB1[k]);
            int nodeA2 = dag.childNode(marker, edgesA2[k]);
            int nodeB2 = dag.childNode(marker, edgesB2[k]);
            nodes.sumUpdate(nodeAB1, nodeA2, nodeB2, fwdValues[k]);
        }
    }

    /**
     * Sets the Baum backward algorithm values for this level of the HMM
     * and stores the parent node trio values in the specified
     * {@code nodes} parameter. When the method call returns, the
     * ${@code nodes} parameter will be reset to the parent node trio values
     * for this level of the HMM.
     *
     * @param nodes parent node trio values at the next level of HMM
     *
     * @throws IndexOutOfBoundsException if any node in any node trio with
     * non-zero value is not a valid child node at this level of the HMM
     * @throws NullPointerException if {@code nodes == null}
     */
    public void setBackwardValues(DuoNodes nodes) {
        for (int j=0; j<size; ++j) {
            int nodeAB1 = dag.childNode(marker, edgesAB1[j]);
            int nodeA2 = dag.childNode(marker, edgesA2[j]);
            int nodeB2 = dag.childNode(marker, edgesB2[j]);
            float backwardValue = nodes.value(nodeAB1, nodeA2, nodeB2);
            bwdValues[j] = backwardValue;
            bwdValueSum += backwardValue;
        }
        nodes.clear();
        float gtProbsSum = 0f;
        for (int j=0; j<size; ++j) {
            bwdValues[j] /= bwdValueSum;
            int symbolAB1 = symbolAB1(j);
            int symbolA2 = symbolA2(j);
            int symbolB2 = symbolB2(j);
            float tpAB1 = dag.condEdgeProb(marker, edgesAB1[j]);
            float tpA2 = dag.condEdgeProb(marker, edgesA2[j]);
            float tpB2 = dag.condEdgeProb(marker, edgesB2[j]);

            float stateProb = fwdValues[j] * bwdValues[j];
            int gtIndexA = BasicGL.genotype(symbolAB1, symbolA2);
            int gtIndexB = BasicGL.genotype(symbolAB1, symbolB2);
            // gtProbs[AB] assumed to be initialized in setForwardValues() method
            gtProbsA[gtIndexA] += stateProb;
            gtProbsB[gtIndexB] += stateProb;
            gtProbsSum += stateProb;

            float epA = gl.gl(marker, sampleA, symbolAB1, symbolA2);
            float epB = gl.gl(marker, sampleB, symbolAB1, symbolB2);
            float bwdValue = bwdValues[j] * (tpAB1 * tpA2 * tpB2) * (epA*epB);
            if (bwdValue < MIN_VALUE && bwdValues[j] > 0f) {
                bwdValue = MIN_VALUE;
            }
            int pnAB1 = dag.parentNode(marker, edgesAB1[j]);
            int pnA2 = dag.parentNode(marker, edgesA2[j]);
            int pnB2 = dag.parentNode(marker, edgesB2[j]);
            nodes.sumUpdate(pnAB1, pnA2, pnB2, bwdValue);
        }
        for (int j=0; j<nGenotypes; ++j) {
            gtProbsA[j] /= gtProbsSum;
            gtProbsB[j] /= gtProbsSum;
        }
    }

    /**
     * Returns the directed acyclic graph that determines the transition
     * probabilities.
     * @return the directed acyclic graph that determines the transition
     * probabilities
     */
    public Dag dag() {
        return dag;
    }

    /**
     * Returns the emission probabilities.
     * @return the emission probabilities
     */
    public GL gl() {
        return gl;
    }

    /**
     * Return the level of the HMM.
     * @return the level of the HMM
     */
    public int marker() {
        return marker;
    }

    /**
     * Return the number of possible genotypes at this level of the HMM.
     * @return the number of possible genotypes at this level of the HMM
     */
    public int nGenotypes() {
        return nGenotypes;
    }

    /**
     * Returns the specified posterior genotype probability for the parent.
     * Returns 0 if the Baum backward probabilities have not been set.
     * @param gt a genotype index
     * @return the specified posterior genotype probability for the parent
     * @throws IndexOutOfBoundsException if
     * {@code gt < 0 || gt >= this.nGenotypes()}
     */
    public float gtProbsA(int gt) {
        checkGT(gt);
        return gtProbsA[gt];
    }

    /**
     * Returns the specified posterior genotype probability for the offspring.
     * Returns 0 if the Baum backward probabilities have not been set.
     * @param gt a genotype index
     * @return the specified posterior genotype probability for the offspring
     * @throws IndexOutOfBoundsException if
     * {@code gt < 0 || gt >= this.nGenotypes()}
     */
    public float gtProbsB(int gt) {
        checkGT(gt);
        return gtProbsB[gt];
    }

    private void checkGT(int gt) {
        if (gt >= nGenotypes) {
            throw new IllegalArgumentException(String.valueOf(gt));
        }
    }

    /**
     * Return the number of states with nonzero forward probability at
     * this level of the HMM.
     *
     * @return the number of states with nonzero forward probability at
     * this level of the HMM
     */
    public int size() {
        return size;
    }

    private void checkIndex(int state) {
        if (state >= size) {
            throw new IndexOutOfBoundsException(String.valueOf(size));
        }
    }

    /**
     * Returns the first edge of the specified HMM state with nonzero forward
     * probability.
     * @param state an index of a HMM state with nonzero forward probability
     * @return the first edge of the specified HMM state with nonzero forward
     * probability
     *
     * @throws IndexOutOfBoundsException if
     * {@code state < 0 || state >= this.size()}
     */
    public int edgeAB1(int state) {
        checkIndex(state);
        return edgesAB1[state];
    }

    /**
     * Returns the second edge of the specified HMM state with nonzero forward
     * probability.
     * @param state an index of a HMM state with nonzero forward probability
     * @return the second edge of the specified HMM state with nonzero forward
     * probability
     *
     * @throws IndexOutOfBoundsException if
     * {@code state < 0 || state >= this.size()}
     */
    public int edgeA2(int state) {
        checkIndex(state);
        return edgesA2[state];
    }

    /**
     * Returns the third edge of the specified HMM state with nonzero forward
     * probability.
     * @param state an index of a HMM state with nonzero forward probability
     * @return the third edge of the specified HMM state with nonzero forward
     * probability
     *
     * @throws IndexOutOfBoundsException if
     * {@code state < 0 || state >= this.size()}
     */
    public int edgeB2(int state) {
        checkIndex(state);
        return edgesB2[state];
    }

    /**
     * Returns the parent node of the first edge of the specified HMM state
     * with nonzero forward probability.
     *
     * @param state an index of a HMM state with nonzero forward probability
     * @return the parent node of the first edge of the specified HMM state
     * with nonzero forward probability
     *
     * @throws IndexOutOfBoundsException if
     * {@code state < 0 || state >= this.size()}
     */
    public int parentNodeAB1(int state) {
        checkIndex(state);
        return dag.parentNode(marker, edgesAB1[state]);
    }

    /**
     * Returns the parent node of the second edge of the specified HMM state
     * with nonzero forward probability.
     *
     * @param state an index of a HMM state with nonzero forward probability
     * @return the parent node of the second edge of the specified HMM state
     * with nonzero forward probability
     *
     * @throws IndexOutOfBoundsException if
     * {@code state < 0 || state >= this.size()}
     */
    public int parentNodeA2(int state) {
        checkIndex(state);
        return dag.parentNode(marker, edgesA2[state]);
    }

    /**
     * Returns the parent node of the third edge of the specified HMM state
     * with nonzero forward probability.
     *
     * @param state an index of a HMM state with nonzero forward probability
     * @return the parent node of the third edge of the specified HMM state
     * with nonzero forward probability
     *
     * @throws IndexOutOfBoundsException if
     * {@code state < 0 || state >= this.size()}
     */
    public int parentNodeB2(int state) {
        checkIndex(state);
        return dag.parentNode(marker, edgesB2[state]);
    }

    /**
     * Returns the child node of the first edge of the specified HMM state
     * with nonzero forward probability.
     *
     * @param state an index of a HMM state with nonzero forward probability
     * @return the child node of the first edge of the specified HMM state
     * with nonzero forward probability
     *
     * @throws IndexOutOfBoundsException if
     * {@code state < 0 || state >= this.size()}
     */
    public int childNodeAB1(int state) {
        checkIndex(state);
        return dag.childNode(marker, edgesAB1[state]);
    }

    /**
     * Returns the child node of the second edge of the specified HMM state
     * with nonzero forward probability.
     *
     * @param state an index of a HMM state with nonzero forward probability
     * @return the child node of the second edge of the specified HMM state
     * with nonzero forward probability
     *
     * @throws IndexOutOfBoundsException if
     * {@code state < 0 || state >= this.size()}
     */
    public int childNodeA2(int state) {
        checkIndex(state);
        return dag.childNode(marker, edgesA2[state]);
    }

    /**
     * Returns the child node of the third edge of the specified HMM state
     * with nonzero forward probability.
     *
     * @param state an index of a HMM state with nonzero forward probability
     * @return the child node of the third edge of the specified HMM state
     * with nonzero forward probability
     *
     * @throws IndexOutOfBoundsException if
     * {@code state < 0 || state >= this.size()}
     */
    public int childNodeB2(int state) {
        checkIndex(state);
        return dag.childNode(marker, edgesB2[state]);
    }

    /**
     * Returns the symbol for the first edge of the specified HMM state
     * with nonzero forward probability.
     *
     * @param state an index of a HMM state with nonzero forward probability
     * @return the symbol for the first edge of the specified HMM state
     * with nonzero forward probability
     *
     * @throws IndexOutOfBoundsException if
     * {@code state < 0 || state >= this.size()}
     */
    public int symbolAB1(int state) {
        return dag.symbol(marker, edgeAB1(state));
    }

    /**
     * Returns the symbol for the second edge of the specified HMM state
     * with nonzero forward probability.
     *
     * @param state an index of a HMM state with nonzero forward probability
     * @return the symbol for the second edge of the specified HMM state
     * with nonzero forward probability
     *
     * @throws IndexOutOfBoundsException if
     * {@code state < 0 || state >= this.size()}
     */
    public int symbolA2(int state) {
        return dag.symbol(marker, edgeA2(state));
    }

    /**
     * Returns the symbol for the third edge of the specified HMM state
     * with nonzero forward probability.
     *
     * @param state an index of a HMM state with nonzero forward probability
     * @return the symbol for the third edge of the specified HMM state
     * with nonzero forward probability
     *
     * @throws IndexOutOfBoundsException if
     * {@code state < 0 || state >= this.size()}
     */
    public int symbolB2(int state) {
        return dag.symbol(marker, edgeB2(state));
    }

    /**
     * Returns the normalized forward value for the specified HMM state
     * with nonzero forward probability.
     * The normalized forward value is obtained by dividing the
     * forward value by the sum of the forward values at this level
     * of the HMM.
     *
     * @param state an index of a HMM state with nonzero forward probability
     *
     * @return the normalized forward value for the specified HMM state
     * with nonzero forward probability
     *
     * @throws IndexOutOfBoundsException if
     * {@code state < 0 || state >= this.size()}
     */
    public float forwardValue(int state) {
        checkIndex(state);
        return fwdValues[state];
    }

    /**
     * Returns the normalized backward value for the specified HMM state
     * with nonzero forward probability.
     * The normalized backward value is obtained by dividing the
     * backward value by the sum of the backward values at this level
     * of the HMM.
     *
     * @param state an index of a state with nonzero forward probability
     *
     * @return the normalized backward value for the specified HMM state
     * with nonzero forward probability
     *
     * @throws IndexOutOfBoundsException if
     * {@code state < 0 || state >= this.size()}
     */
    public float backwardValue(int state) {
        checkIndex(state);
        return bwdValues[state];
    }

    /**
     * Returns the sum of the forward values at this level of the HMM
     * when the forward values are computed using forward values
     * from the previous level that are normalized to sum to 1.
     * @return the sum of the forward values at this level of the HMM
     */
    public float forwardValuesSum() {
        return fwdValueSum;
    }

    /**
     * Returns the sum of the backward values at this level of the HMM
     * when the backward values are computed using backward
     * values from the next level that are normalized to sum to 1.
     * @return the sum of the backward values at this level of the HMM
     */
    public float backwardValuesSum() {
        return bwdValueSum;
    }

    /**
     * Returns a string description of {@code this}.  The exact details
     * of the description are unspecified and subject to change.
     *
     * @return a string description of {@code this}.
     */
    @Override
    public String toString() {
        String space = " ";
        String sep = " | ";
        StringBuilder sb = new StringBuilder(100);
        sb.append("level=");
        sb.append(marker);
        sb.append(" size=");
        sb.append(size);
        sb.append(" forwardValuesSum=");
        sb.append(fwdValueSum);
        sb.append(" backwardSum=");
        sb.append(bwdValueSum);
        for (int j=0; j<size; ++j) {
            sb.append(sep);
            sb.append("j=");
            sb.append(j);
            sb.append(": ");
            sb.append( (int) edgeAB1(j));
            sb.append(space);
            sb.append( (int) edgeA2(j));
            sb.append(space);
            sb.append( (int) edgeB2(j));
            sb.append(space);
            sb.append(forwardValue(j));
            sb.append(space);
            sb.append(backwardValue(j));
        }
        sb.append(sep);
        return sb.toString();
    }

    /*
     * Increases the state capacity of array fields as necessary
     * to be greater than or equal to the specified minimum capacity.
     *
     * @param minCapacity the desired minimum state capacity
     */
    private void ensureCapacity(int minCapacity) {
        if (minCapacity > capacity) {
            capacity = (capacity * 3)/2 + 1;
            if (capacity < minCapacity) {
                capacity = minCapacity;
            }
            edgesAB1 = Arrays.copyOf(edgesAB1, capacity);
            edgesA2 = Arrays.copyOf(edgesA2, capacity);
            edgesB2 = Arrays.copyOf(edgesB2, capacity);
            fwdValues = Arrays.copyOf(fwdValues, capacity);
            bwdValues = Arrays.copyOf(bwdValues, capacity);
        }
    }
}
