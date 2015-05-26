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
 * Class {@code SingleBaumLevel} computes forward and backward Baum
 * values at a level of a hidden Markov model (HMM) whose states are
 * ordered edge pairs of a leveled directed acyclic graph (DAG).
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class SingleBaumLevel {

    private static final int INITIAL_CAPACITY = 400;
    private static final double MIN_VALUE = 100*Double.MIN_VALUE;
    private final Dag dag;
    private final GL gl;

    private int marker = -1;
    private int sample = -1;
    private int size=0;

    private int capacity = INITIAL_CAPACITY;
    private int[] edges1 = new int[INITIAL_CAPACITY];
    private int[] edges2 = new int[INITIAL_CAPACITY];
    private double[] fwdValues = new double[INITIAL_CAPACITY];
    private double[] bwdValues = new double[INITIAL_CAPACITY];
    private double fwdValueSum = 0.0;
    private double bwdValueSum = 0.0;

    private int nGenotypes = 0;
    private double[] gtProbs = new double[3];

    /**
     * Constructs a new {@code SingleBaumLevel} instance.
     * @param dag the directed acyclic graph that the determines transition
     * probabilities.
     * @param gl the emission probabilities.
     * @throws IllegalArgumentException if
     * {@code dag.markers().equals(gl.markers())==false}
     * @throws NullPointerException if {@code dag==null || gl==null}
     */
    public SingleBaumLevel(Dag dag, GL gl) {
        if (dag.markers().equals(gl.markers())==false) {
            throw new IllegalArgumentException("marker inconsistency");
        }
        this.dag = dag;
        this.gl = gl;
    }

    /**
     * Initializes the node pair values for the Baum forward algorithm.
     *
     * @param nodes the node pair values to be initialized.
     */
    public static void initializeNodes(SingleNodes nodes) {
        nodes.clear();
        nodes.sumUpdate(0, 0, 1.0);
    }

    /**
     * Sets the Baum forward algorithm values for this level of the HMM
     * and records the child node pair values in the specified
     * {@code nodes} parameter.
     *
     * @param nodes child node pair values at the previous level of HMM.  When
     * the method call returns, this parameter will be reset to the child
     * node pair values for this level of the HMM.
     * @param marker the level of the HMM at which the Baum forward algorithm
     * values will be computed.
     * @param sample a sample index.
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker<0 || marker>=this.dag().nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code sample<0 || sample>=this.gl().nSamples()}
     * @throws IndexOutOfBoundsException if either node in any node pair with
     * non-zero value is not a valid parent node at the specified level of the
     * HMM
     * @throws NullPointerException if {@code nodes==null}
     */
    public void setForwardValues(SingleNodes nodes, int marker, int sample) {
        this.marker = marker;
        this.sample = sample;
        this.nGenotypes = gl.marker(marker).nGenotypes();
        this.size = 0;
        this.fwdValueSum = 0.0;
        this.bwdValueSum = 0.0;
        initializeGtProbs(); // called here due to grProbs() contract
        setStates(nodes);
        setChildNodes(nodes);
    }

    private void initializeGtProbs() {
        if (gtProbs.length < nGenotypes) {
            int newLength = Math.max(nGenotypes, (3*gtProbs.length/2 + 1));
            gtProbs = new double[newLength];
        }
        else {
            for (int j=0; j<nGenotypes; ++j) {
                gtProbs[j] = 0.0;
            }
        }
    }

    private void setStates(SingleNodes nodes) {
        double valueSum = 0.0;
        for (int j=0, n=nodes.size(); j<n; ++j) {
            int node1 = nodes.enumNode1(j);
            int node2 = nodes.enumNode2(j);
            for (int i1=0, nI1=dag.nOutEdges(marker, node1); i1<nI1; ++i1) {
                int edge1 = dag.outEdge(marker, node1, i1);
                byte symbol1 = dag.symbol(marker, edge1);
                for (int i2=0, nI2=dag.nOutEdges(marker, node2); i2<nI2; ++i2) {
                    int edge2 = dag.outEdge(marker, node2, i2);
                    byte symbol2 = dag.symbol(marker, edge2);
                    float ep = gl.gl(marker, sample, symbol1, symbol2);
                    if (ep > 0.0) {
                        if (size == capacity) {
                            ensureCapacity(size+1);
                        }
                        edges1[size] = edge1;
                        edges2[size] = edge2;
                        double tp1 = dag.condEdgeProb(marker, edge1);
                        double tp2 = dag.condEdgeProb(marker, edge2);
                        double fwdValue = ep * nodes.enumValue(j) * (tp1 * tp2);
                        if (fwdValue<MIN_VALUE && nodes.enumValue(j) > 0.0) {
                            fwdValue = MIN_VALUE;
                        }
                        fwdValues[size++] = fwdValue;
                        valueSum += fwdValue;
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
     * Stores the Baum forward algorithm child node pair values for this
     * level of the HMM in the specified {@code SingleNodes} object.
     *
     * @param nodes the node pair values that will be set.
     *
     * @throws NullPointerException if {@code nodes==null}
     */
    public void setChildNodes(SingleNodes nodes) {
        nodes.clear();
        for (int k=0; k<size; ++k) {
            int node1 = dag.childNode(marker, edges1[k]);
            int node2 = dag.childNode(marker, edges2[k]);
            nodes.sumUpdate(node1, node2, fwdValues[k]);
        }
    }

    /**
     * Initializes the node pair values for the Baum backward algorithm.
     *
     * @param nodes the node pair values to be initialized.
     * @throws NullPointerException if {@code nodes==null}
     */
    public void setInitialBackwardValues(SingleNodes nodes) {
        nodes.clear();
        for (int j=0; j<size; ++j) {
            int node1 = dag.childNode(marker, edges1[j]);
            int node2 = dag.childNode(marker, edges2[j]);
            nodes.maxUpdate(node1, node2, 1.0);
        }
        setBackwardValues(nodes);
    }

    /**
     * Sets the Baum backward algorithm values for this level of the HMM
     * and stores the parent node pair values in the specified
     * {@code nodes} parameter.
     *
     * @param nodes parent node pair values at the next level of HMM.  When
     * the method call returns, this parameter will be reset to the parent
     * node pair values for this level of the HMM.
     *
     * @throws IndexOutOfBoundsException if either node in any node pair with
     * non-zero value is not a valid child node at the {@code this.marker()}
     * level of the HMM
     * @throws NullPointerException if {@code nodes==null}
     */
    public void setBackwardValues(SingleNodes nodes) {
        bwdValueSum = 0.0;
        double gtProbsSum = 0.0;
        for (int j=0; j<size; ++j) {
            int node1 = dag.childNode(marker, edges1[j]);
            int node2 = dag.childNode(marker, edges2[j]);
            double backwardValue = nodes.value(node1, node2);
            bwdValues[j] = backwardValue;
            bwdValueSum += backwardValue;
        }
        nodes.clear();
        for (int j=0; j<size; ++j) {
            bwdValues[j] /= bwdValueSum;
            int edge1 = edges1[j];
            int edge2 = edges2[j];
            byte symbol1 = symbol1(j);
            byte symbol2 = symbol2(j);
            int node1 = dag.parentNode(marker, edge1);
            int node2 = dag.parentNode(marker, edge2);
            double tp1 = dag.condEdgeProb(marker, edge1);
            double tp2 = dag.condEdgeProb(marker, edge2);

            double stateProb = fwdValues[j] * bwdValues[j];
            int gtIndex = BasicGL.genotype(symbol1, symbol2);
            // gtProbs initialized in setForwardValues() method
            gtProbs[gtIndex] += stateProb;
            gtProbsSum += stateProb;

            double ep = gl.gl(marker, sample, symbol1, symbol2);
            double bwdValue = bwdValues[j] * (tp1 * tp2) * ep;
            if (bwdValue < MIN_VALUE && bwdValues[j]>0.0) {
                bwdValue = MIN_VALUE;
            }
            nodes.sumUpdate(node1, node2, bwdValue);
        }
        for (int j=0; j<nGenotypes; ++j) {
            gtProbs[j] /= gtProbsSum;
        }
    }

    /**
     * Returns the directed acyclic graph that determines the transition
     * probabilities.
     * @return the directed acyclic graph that determines the transition
     * probabilities.
     */
    public Dag dag() {
        return dag;
    }

    /**
     * Returns the emission probabilities.
     * @return the emission probabilities.
     */
    public GL gl() {
        return gl;
    }

    /**
     * Return the level of the HMM.
     * @return the level of the HMM.
     */
    public int marker() {
        return marker;
    }

    /**
     * Return the number of possible genotypes at this level of the HMM.
     * @return the number of possible genotypes at this level of the HMM.
     */
    public int nGenotypes() {
        return nGenotypes;
    }

    /**
     * Returns the specified posterior genotype probability.  Returns 0
     * if the Baum backward values have not been set.
     * @param gt a genotype index.
     * @return the specified posterior genotype probability.
     * @throws IndexOutOfBoundsException if
     * {@code gt<0 || gt>=this.nGenotypes()}
     */
    public double gtProbs(int gt) {
        if (gt >= nGenotypes) {
            throw new IllegalArgumentException(String.valueOf(gt));
        }
        return gtProbs[gt];
    }

    /**
     * Return the number of states with nonzero forward probability at
     * this level of the HMM.
     *
     * @return the number of states with nonzero forward probability at
     * this level of the HMM.
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
     * Returns the DAG level edge index for the first edge of the
     * specified HMM state with nonzero forward probability.
     * @param state an index of a HMM state at this level with nonzero
     * forward probability.
     * @return the DAG level edge index for the first edge of the
     * specified HMM state with nonzero forward probability.
     *
     * @throws IndexOutOfBoundsException if
     * {@code state<0 || state>=this.size()}
     */
    public int edge1(int state) {
        checkIndex(state);
        return edges1[state];
    }

    /**
     * Returns the DAG level edge index for the second edge of the
     * specified HMM state with nonzero forward probability.
     * @param state an index of a HMM state at this level with nonzero
     * forward probability.
     * @return the DAG level edge index for the second edge of the
     * specified HMM state with nonzero forward probability.
     *
     * @throws IndexOutOfBoundsException if
     * {@code state<0 || state>=this.size()}
     */
    public int edge2(int state) {
        checkIndex(state);
        return edges2[state];
    }

    /**
     * Returns the DAG level parent node index for the parent node of the
     * first edge of the specified HMM state with nonzero forward probability.
     *
     * @param state an index of a HMM state at this level with nonzero
     * forward probability.
     * @return the DAG level parent node index for the parent node of the
     * first edge of the specified HMM state with nonzero forward probability.
     *
     * @throws IndexOutOfBoundsException if
     * {@code state<0 || state>=this.size()}
     */
    public int parentNode1(int state) {
        checkIndex(state);
        return dag.parentNode(marker, edges1[state]);
    }

    /**
     * Returns the DAG level parent node index for the parent node of the
     * second edge of the specified HMM state with nonzero forward probability.
     *
     * @param state an index of a HMM state at this level with nonzero
     * forward probability.
     * @return the DAG level parent node index for the parent node of the
     * second edge of the specified HMM state with nonzero forward probability.
     *
     * @throws IndexOutOfBoundsException if
     * {@code state<0 || state>=this.size()}
     */
    public int parentNode2(int state) {
        checkIndex(state);
        return dag.parentNode(marker, edges2[state]);
    }

    /**
     * Returns the DAG level child node index for the child node of the
     * first edge of the specified HMM state with nonzero forward probability.
     *
     * @param state an index of a HMM state at this level with nonzero
     * forward probability.
     * @return the DAG level child node index for the child node of the
     * first edge of the specified HMM state with nonzero forward probability.
     *
     * @throws IndexOutOfBoundsException if
     * {@code state<0 || state>=this.size()}
     */
    public int childNode1(int state) {
        checkIndex(state);
        return dag.childNode(marker, edges1[state]);
    }

    /**
     * Returns the DAG level child node index for the child node of the
     * second edge of the specified HMM state with nonzero forward probability.
     *
     * @param state an index of a HMM state at this level with nonzero
     * forward probability.
     * @return the DAG level child node index for the child node of the
     * second edge of the specified HMM state with nonzero forward probability.
     *
     * @throws IndexOutOfBoundsException if
     * {@code state<0 || state>=this.size()}
     */
    public int childNode2(int state) {
        checkIndex(state);
        return dag.childNode(marker, edges2[state]);
    }

    /**
     * Returns the symbol for the first edge of the specified HMM state
     * with nonzero forward probability.
     *
     * @param state an index of a HMM state at this level with nonzero
     * forward probability.
     * @return the symbol for the first edge of the specified HMM state
     * with nonzero forward probability.
     *
     * @throws IndexOutOfBoundsException if
     * {@code state<0 || state>=this.size()}
     */
    public byte symbol1(int state) {
        return dag.symbol(marker, edge1(state));
    }

    /**
     * Returns the symbol for the second edge of the specified HMM state
     * with nonzero forward probability.
     *
     * @param state an index of a HMM state at this level with nonzero
     * forward probability.
     * @return the symbol for the second edge of the specified HMM state
     * with nonzero forward probability.
     *
     * @throws IndexOutOfBoundsException if
     * {@code state<0 || state>=this.size()}
     */
    public byte symbol2(int state) {
        return dag.symbol(marker, edge2(state));
    }

    /**
     * Returns the normalized forward value for the specified HMM state
     * with nonzero forward probability.
     * The normalized forward value is obtained by dividing the
     * forward value by the sum of the forward values at this level
     * of the HMM.
     *
     * @param state an index of a HMM state at this level with nonzero
     * forward probability.
     *
     * @return the normalized forward value for the specified HMM state
     * with nonzero forward probability.
     *
     * @throws IndexOutOfBoundsException if
     * {@code state<0 || state>=this.size()}
     */
    public double forwardValue(int state) {
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
     * @param state an index of a state with nonzero backward value.
     *
     * @return the normalized backward value for the specified HMM state
     * with nonzero forward probability.
     *
     * @throws IndexOutOfBoundsException if
     * {@code state<0 || state>=this.size()}
     */
    public double backwardValue(int state) {
        checkIndex(state);
        return bwdValues[state];
    }

    /**
     * Returns the sum of the forward values at this level of the HMM
     * when the forward values are computed using normalized forward values
     * from the previous level that are normalized to sum to 1.
     * @return the sum of the forward values at this level of the HMM.
     */
    public double forwardValuesSum() {
        return fwdValueSum;
    }

    /**
     * Returns the sum of the backward values at this level of the HMM
     * when the backward values are computed using normalized backward
     * values from the next level that are normalized to sum to 1.
     * @return the sum of the backward values at this level of the HMM.
     */
    public double backwardValuesSum() {
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
            sb.append( (int) edge1(j));
            sb.append(space);
            sb.append( (int) edge2(j));
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
     * @param minCapacity the desired minimum state capacity.
     */
    private void ensureCapacity(int minCapacity) {
        if (minCapacity >capacity) {
            capacity = (capacity * 3)/2 + 1;
            if (capacity < minCapacity) {
                capacity = minCapacity;
            }
            edges1 = Arrays.copyOf(edges1, capacity);
            edges2 = Arrays.copyOf(edges2, capacity);
            fwdValues = Arrays.copyOf(fwdValues, capacity);
            bwdValues = Arrays.copyOf(bwdValues, capacity);
        }
    }
}
