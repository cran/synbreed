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
 * Class {@code TrioBaumLevel} computes forward and backward Baum
 * values at a level of a hidden Markov model (HMM) whose states are
 * ordered edge quartets of a leveled directed acyclic graph (DAG).
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class TrioBaumLevel {

    private static final int INITIAL_CAPACITY = 400;
    private static final double MIN_VALUE = 100*Double.MIN_VALUE;
    private final Dag dag;
    private final GL gl;

    private int marker = -1;
    private int sampleA = -1;
    private int sampleB = -1;
    private int sampleC = -1;
    private int size=0;

    private int capacity = INITIAL_CAPACITY;
    private int[] edgesA1 = new int[INITIAL_CAPACITY];
    private int[] edgesA2 = new int[INITIAL_CAPACITY];
    private int[] edgesB1 = new int[INITIAL_CAPACITY];
    private int[] edgesB2 = new int[INITIAL_CAPACITY];
    private double[] fwdValues = new double[INITIAL_CAPACITY];
    private double[] bwdValues = new double[INITIAL_CAPACITY];
    private double fwdValueSum = 0.0;
    private double bwdValueSum = 0.0;

    private int nGenotypes = 0;
    private double[] gtProbsA = new double[3];
    private double[] gtProbsB = new double[3];
    private double[] gtProbsC = new double[3];


    /**
     * Constructs a new {@code TrioBaumLevel} instance.
     * @param dag the directed acyclic graph that the determines transition
     * probabilities.
     * @param gl the emission probabilities.
     * @throws IllegalArgumentException if
     * {@code dag.markers().equals(gl.markers())==false}
     * @throws NullPointerException if {@code dag==null || gl==null}
     */
    public TrioBaumLevel(Dag dag, GL gl) {
        if (dag.markers().equals(gl.markers())==false) {
            throw new IllegalArgumentException("marker inconsistency");
        }
        this.dag = dag;
        this.gl = gl;
    }

    /**
     * Initializes the node quartet values for the Baum forward algorithm.
     *
     * @param nodes the node quartet values to be initialized.
     */
    public static void initializeNodes(TrioNodes nodes) {
        nodes.clear();
        nodes.sumUpdate(0, 0, 0, 0, 1.0);
    }

    /**
     * Sets the Baum forward algorithm values for this level of the HMM
     * and records the child node quartet values in the specified
     * {@code nodes} parameter.
     *
     * @param nodes child node quartet values at the previous level of HMM.
     * When the method call returns, this parameter will be reset to the child
     * node quartet values for this level of the HMM.
     * @param marker the level of the HMM at which the Baum forward algorithm
     * probabilities will be computed.
     * @param sampleA the father's sample index.
     * @param sampleB the mother's sample index.
     * @param sampleC the offspring's sample index.
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker<0 || marker>=this.dag().nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code sampleA<0 || sampleA>=this.gl().nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code sampleB<0 || sampleB>=this.gl().nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code sampleC<0 || sampleC>=this.gl().nSamples()}
     * @throws IndexOutOfBoundsException if any node in any node quartet with
     * non-zero value is not a valid node parent at the specified level of the
     * HMM
     * @throws NullPointerException if {@code nodes==null}
     */
    public void setForwardValues(TrioNodes nodes, int marker, int sampleA,
            int sampleB, int sampleC) {
        this.marker = marker;
        this.sampleA = sampleA;
        this.sampleB = sampleB;
        this.sampleC = sampleC;
        this.nGenotypes = gl.marker(marker).nGenotypes();
        this.size = 0;
        this.fwdValueSum = 0.0;
        this.bwdValueSum = 0.0;
        initializeGtProbs(); // called here due to gtProbs[ABC]() contract
        setStates(nodes);
        setChildNodes(nodes);
    }

    private void initializeGtProbs() {
        if (gtProbsA.length < nGenotypes) {
            int newLength = Math.max(nGenotypes, (3*gtProbsA.length/2 + 1));
            gtProbsA = new double[newLength];
            gtProbsB = new double[newLength];
            gtProbsC = new double[newLength];
        }
        else {
            for (int j=0; j<nGenotypes; ++j) {
                gtProbsA[j] = 0.0;
                gtProbsB[j] = 0.0;
                gtProbsC[j] = 0.0;
            }
        }
    }

     private void setStates(TrioNodes nodes) {
        double valueSum = 0.0;
        for (int j=0, n=nodes.size(); j<n; ++j) {
            int nodeA1 = nodes.enumNodeA1(j);
            int nodeA2 = nodes.enumNodeA2(j);
            int nodeB1 = nodes.enumNodeB1(j);
            int nodeB2 = nodes.enumNodeB2(j);
            double nodeValue = nodes.enumValue(j);
            for (int a1=0, nA1=dag.nOutEdges(marker, nodeA1); a1<nA1; ++a1) {
                int edgeA1 = dag.outEdge(marker, nodeA1, a1);
                byte symbolA1 = dag.symbol(marker, edgeA1);
                for (int a2=0, nA2=dag.nOutEdges(marker, nodeA2); a2<nA2; ++a2) {
                    int edgeA2 = dag.outEdge(marker, nodeA2, a2);
                    byte symbolA2 = dag.symbol(marker, edgeA2);
                    float epA = gl.gl(marker, sampleA, symbolA1, symbolA2);
                    if (epA > 0.0f) {
                        for (int b1=0, nB1=dag.nOutEdges(marker, nodeB1); b1<nB1; ++b1) {
                            int edgeB1 = dag.outEdge(marker, nodeB1, b1);
                            byte symbolB1 = dag.symbol(marker, edgeB1);
                            for (int b2=0, nB2=dag.nOutEdges(marker, nodeB2); b2<nB2; ++b2) {
                                int edgeB2 = dag.outEdge(marker, nodeB2, b2);
                                byte symbolB2 = dag.symbol(marker, edgeB2);
                                float epB = gl.gl(marker, sampleB, symbolB1, symbolB2);
                                if (epB > 0.0f) {
                                    float epC = gl.gl(marker, sampleC, symbolA1, symbolB1);
                                    if (epC > 0.0f) {
                                        if (size == capacity) {
                                            ensureCapacity(size+1);
                                        }
                                        double tpA1 = dag.condEdgeProb(marker, edgeA1);
                                        double tpA2 = dag.condEdgeProb(marker, edgeA2);
                                        double tpB1 = dag.condEdgeProb(marker, edgeB1);
                                        double tpB2 = dag.condEdgeProb(marker, edgeB2);
                                        double fwdValue = (epA * epB * epC) * nodeValue
                                                * (tpA1 * tpA2 * tpB1 * tpB2);
                                        if (fwdValue<MIN_VALUE && nodeValue > 0.0) {
                                            fwdValue = MIN_VALUE;
                                        }
                                        edgesA1[size] = edgeA1;
                                        edgesA2[size] = edgeA2;
                                        edgesB1[size] = edgeB1;
                                        edgesB2[size] = edgeB2;
                                        fwdValues[size++] = fwdValue;
                                        valueSum += fwdValue;
                                    }
                                }
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
     * Stores the Baum forward algorithm child node quartet values for this
     * level of the HMM in the specified {@code TrioNodes} object.
     *
     * @param nodes the node quartet values that will be set.
     *
     * @throws NullPointerException if {@code nodes==null}
     */
    public void setChildNodes(TrioNodes nodes) {
        nodes.clear();
        for (int k=0; k<size; ++k) {
            int nodeA1 = dag.childNode(marker, edgesA1[k]);
            int nodeA2 = dag.childNode(marker, edgesA2[k]);
            int nodeB1 = dag.childNode(marker, edgesB1[k]);
            int nodeB2 = dag.childNode(marker, edgesB2[k]);
            nodes.sumUpdate(nodeA1, nodeA2, nodeB1, nodeB2, fwdValues[k]);
        }
    }

    /**
     * Initializes the node quartet values for the Baum backward algorithm.
     *
     * @param nodes the node quartet values to be initialized.
     * @throws NullPointerException if {@code nodes==null}
     */
    public void setInitialBackwardValues(TrioNodes nodes) {
        nodes.clear();
        for (int j=0; j<size; ++j) {
            int nodeA1 = dag.childNode(marker, edgesA1[j]);
            int nodeA2 = dag.childNode(marker, edgesA2[j]);
            int nodeB1 = dag.childNode(marker, edgesB1[j]);
            int nodeB2 = dag.childNode(marker, edgesB2[j]);
            nodes.maxUpdate(nodeA1, nodeA2, nodeB1, nodeB2, 1.0);
        }
        setBackwardValues(nodes);
    }

    /**
     * Sets the Baum backward algorithm values for this level of the HMM
     * and stores the parent node quartet values in the specified
     * {@code nodes} parameter.
     *
     * @param nodes parent node quartet values at the next level of HMM.  When
     * the method call returns, this parameter will be reset to the parent
     * node quartet values for this level of the HMM.
     *
     * @throws IndexOutOfBoundsException if any node in any node quartet with
     * non-zero value is not a valid child node at the {@code this.marker()}
     * level of the HMM
     * @throws NullPointerException if {@code nodes==null}
     */
    public void setBackwardValues(TrioNodes nodes) {
        bwdValueSum = 0.0;
        double gtProbsSum = 0.0;
        for (int j=0; j<size; ++j) {
            int nodeA1 = dag.childNode(marker, edgesA1[j]);
            int nodeA2 = dag.childNode(marker, edgesA2[j]);
            int nodeB1 = dag.childNode(marker, edgesB1[j]);
            int nodeB2 = dag.childNode(marker, edgesB2[j]);
            double backwardValue = nodes.value(nodeA1, nodeA2, nodeB1, nodeB2);
            bwdValues[j] = backwardValue;
            bwdValueSum += backwardValue;
        }
        nodes.clear();
        for (int j=0; j<size; ++j) {
            bwdValues[j] /= bwdValueSum;
            byte symbolA1 = symbolA1(j);
            byte symbolA2 = symbolA2(j);
            byte symbolB1 = symbolB1(j);
            byte symbolB2 = symbolB2(j);
            int nodeA1 = dag.parentNode(marker, edgesA1[j]);
            int nodeA2 = dag.parentNode(marker, edgesA2[j]);
            int nodeB1 = dag.parentNode(marker, edgesB1[j]);
            int nodeB2 = dag.parentNode(marker, edgesB2[j]);
            double pA1 = dag.condEdgeProb(marker, edgesA1[j]);
            double pA2 = dag.condEdgeProb(marker, edgesA2[j]);
            double pB1 = dag.condEdgeProb(marker, edgesB1[j]);
            double pB2 = dag.condEdgeProb(marker, edgesB2[j]);

            double stateProb = fwdValues[j] * bwdValues[j];
            int gtIndexA = BasicGL.genotype(symbolA1, symbolA2);
            int gtIndexB = BasicGL.genotype(symbolB1, symbolB2);
            int gtIndexC = BasicGL.genotype(symbolA1, symbolB1);
            // gtProbsA, gtProbsB, gtProbsC initialized in setForwardValues() method
            gtProbsA[gtIndexA] += stateProb;
            gtProbsB[gtIndexB] += stateProb;
            gtProbsC[gtIndexC] += stateProb;
            gtProbsSum += stateProb;

            double emA = gl.gl(marker, sampleA, symbolA1, symbolA2);
            double emB = gl.gl(marker, sampleB, symbolB1, symbolB2);
            double emC = gl.gl(marker, sampleC, symbolA1, symbolB1);
            double bwdValue = bwdValues[j] * (pA1 * pA2 * pB1 * pB2)
                    * (emA * emB * emC);
            if (bwdValue<MIN_VALUE && bwdValues[j]>0.0) {
                bwdValue = MIN_VALUE;
            }
            nodes.sumUpdate(nodeA1, nodeA2, nodeB1, nodeB2, bwdValue);
        }
        for (int j=0; j<nGenotypes; ++j) {
            gtProbsA[j] /= gtProbsSum;
            gtProbsB[j] /= gtProbsSum;
            gtProbsC[j] /= gtProbsSum;
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
     * Returns the specified posterior genotype probability for the father.
     * Returns 0 if the Baum backward probabilities have not been set.
     * @param gt a genotype index.
     * @return the specified posterior genotype probability for the father.
     * @throws IndexOutOfBoundsException if
     * {@code gt<0 || gt>=this.nGenotypes()}
     */
    public double gtProbsA(int gt) {
        checkGT(gt);
        return gtProbsA[gt];
    }

    /**
     * Returns the specified posterior genotype probability for the mother.
     * Returns 0 if the Baum backward probabilities have not been set.
     * @param gt a genotype index.
     * @return the specified posterior genotype probability for the mother.
     * @throws IndexOutOfBoundsException if
     * {@code gt<0 || gt>=this.nGenotypes()}
     */
    public double gtProbsB(int gt) {
        checkGT(gt);
        return gtProbsB[gt];
    }

    /**
     * Returns the specified posterior genotype probability for the offspring.
     * Returns 0 if the Baum backward probabilities have not been set.
     * @param gt a genotype index.
     * @return the specified posterior genotype probability for the offspring.
     * @throws IndexOutOfBoundsException if
     * {@code gt<0 || gt>=this.nGenotypes()}
     */
    public double gtProbsC(int gt) {
        checkGT(gt);
        return gtProbsC[gt];
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
    public int edgeA1(int state) {
        checkIndex(state);
        return edgesA1[state];
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
    public int edgeA2(int state) {
        checkIndex(state);
        return edgesA2[state];
    }

    /**
     * Returns the DAG level edge index for the third edge of the
     * specified HMM state with nonzero forward probability.
     * @param state an index of a HMM state at this level with nonzero
     * forward probability.
     * @return the DAG level edge index for the third edge of the
     * specified HMM state with nonzero forward probability.
     *
     * @throws IndexOutOfBoundsException if
     * {@code state<0 || state>=this.size()}
     */
    public int edgeB1(int state) {
        checkIndex(state);
        return edgesB1[state];
    }

    /**
     * Returns the DAG level edge index for the fourth edge of the
     * specified HMM state with nonzero forward probability.
     * @param state an index of a HMM state at this level with nonzero
     * forward probability.
     * @return the DAG level edge index for the fourth edge of the
     * specified HMM state with nonzero forward probability.
     *
     * @throws IndexOutOfBoundsException if
     * {@code state<0 || state>=this.size()}
     */
    public int edgeB2(int state) {
        checkIndex(state);
        return edgesB2[state];
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
    public int parentNodeA1(int state) {
        checkIndex(state);
        return dag.parentNode(marker, edgesA1[state]);
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
    public int parentNodeA2(int state) {
        checkIndex(state);
        return dag.parentNode(marker, edgesA2[state]);
    }

    /**
     * Returns the DAG level parent node index for the parent node of the
     * third edge of the specified HMM state with nonzero forward probability.
     *
     * @param state an index of a HMM state at this level with nonzero
     * forward probability.
     * @return the DAG level parent node index for the parent node of the
     * third edge of the specified HMM state with nonzero forward probability.
     *
     * @throws IndexOutOfBoundsException if
     * {@code state<0 || state>=this.size()}
     */
    public int parentNodeB1(int state) {
        checkIndex(state);
        return dag.parentNode(marker, edgesB1[state]);
    }

    /**
     * Returns the DAG level parent node index for the parent node of the
     * fourth edge of the specified HMM state with nonzero forward probability.
     *
     * @param state an index of a HMM state at this level with nonzero
     * forward probability.
     * @return the DAG level parent node index for the parent node of the
     * fourth edge of the specified HMM state with nonzero forward probability.
     *
     * @throws IndexOutOfBoundsException if
     * {@code state<0 || state>=this.size()}
     */
    public int parentNodeB2(int state) {
        checkIndex(state);
        return dag.parentNode(marker, edgesB2[state]);
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
    public int childNodeA1(int state) {
        checkIndex(state);
        return dag.childNode(marker, edgesA1[state]);
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
    public int childNodeA2(int state) {
        checkIndex(state);
        return dag.childNode(marker, edgesA2[state]);
    }

    /**
     * Returns the DAG level child node index for the child node of the
     * third edge of the specified HMM state with nonzero forward probability.
     *
     * @param state an index of a HMM state at this level with nonzero
     * forward probability.
     * @return the DAG level child node index for the child node of the
     * third edge of the specified HMM state with nonzero forward probability.
     *
     * @throws IndexOutOfBoundsException if
     * {@code state<0 || state>=this.size()}
     */
    public int childNodeB1(int state) {
        checkIndex(state);
        return dag.childNode(marker, edgesB1[state]);
    }

    /**
     * Returns the DAG level child node index for the child node of the
     * fourth edge of the specified HMM state with nonzero forward probability.
     *
     * @param state an index of a HMM state at this level with nonzero
     * forward probability.
     * @return the DAG level child node index for the child node of the
     * fourth edge of the specified HMM state with nonzero forward probability.
     *
     * @throws IndexOutOfBoundsException if
     * {@code state<0 || state>=this.size()}
     */
    public int childNodeB2(int state) {
        checkIndex(state);
        return dag.childNode(marker, edgesB2[state]);
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
    public byte symbolA1(int state) {
        return dag.symbol(marker, edgeA1(state));
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
    public byte symbolA2(int state) {
        return dag.symbol(marker, edgeA2(state));
    }

    /**
     * Returns the symbol for the third edge of the specified HMM state
     * with nonzero forward probability.
     *
     * @param state an index of a HMM state at this level with nonzero
     * forward probability.
     * @return the symbol for the third edge of the specified HMM state
     * with nonzero forward probability.
     *
     * @throws IndexOutOfBoundsException if
     * {@code state<0 || state>=this.size()}
     */
    public byte symbolB1(int state) {
        return dag.symbol(marker, edgeB1(state));
    }

    /**
     * Returns the symbol for the fourth edge of the specified HMM state
     * with nonzero forward probability.
     *
     * @param state an index of a HMM state at this level with nonzero
     * forward probability.
     * @return the symbol for the fourth edge of the specified HMM state
     * with nonzero forward probability.
     *
     * @throws IndexOutOfBoundsException if
     * {@code state<0 || state>=this.size()}
     */
    public byte symbolB2(int state) {
        return dag.symbol(marker, edgeB2(state));
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
            sb.append( (int) edgeA1(j));
            sb.append(space);
            sb.append( (int) edgeA2(j));
            sb.append(space);
            sb.append( (int) edgeB1(j));
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
     * @param minCapacity the desired minimum state capacity.
     */
    private void ensureCapacity(int minCapacity) {
        if (minCapacity > capacity) {
            capacity = (capacity * 3)/2 + 1;
            if (capacity < minCapacity) {
                capacity = minCapacity;
            }
            edgesA1 = Arrays.copyOf(edgesA1, capacity);
            edgesA2 = Arrays.copyOf(edgesA2, capacity);
            edgesB1 = Arrays.copyOf(edgesB1, capacity);
            edgesB2 = Arrays.copyOf(edgesB2, capacity);
            fwdValues = Arrays.copyOf(fwdValues, capacity);
            bwdValues = Arrays.copyOf(bwdValues, capacity);
        }
    }
}
