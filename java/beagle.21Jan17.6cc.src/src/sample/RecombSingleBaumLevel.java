/*
 * Copyright 2014 Brian L. Browning
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package sample;

import dag.Dag;
import java.util.Arrays;
import vcf.BasicGL;
import vcf.GL;

/**
 * <p>Class {@code RestrictedSingleBaumLevel}  computes forward and backward
 * Baum values at a level of a hidden Markov model (HMM) whose states are
 * ordered edge pairs of a leveled directed acyclic graph (DAG).
 * </p>
 * <p>Instances of class {@code RestrictedSingleBaumLevel} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class RecombSingleBaumLevel {

    private static final int INITIAL_CAPACITY = 400;
    private static final float MIN_VALUE = 100*Float.MIN_VALUE;
    private final SamplerData samplerData;
    private final Dag dag;
    private final GL gl;

    private int marker = -1;
    private int sample = -1;
    private int size = 0;

    private int capacity = INITIAL_CAPACITY;
    private int[] edges1 = new int[INITIAL_CAPACITY];
    private int[] edges2 = new int[INITIAL_CAPACITY];
    private float[] fwdValues = new float[INITIAL_CAPACITY];
    private float[] bwdValues = new float[INITIAL_CAPACITY];
    private float fwdValueSum = 0.0f;
    private float bwdValueSum = 0.0f;

    private int nGenotypes = 0;
    private float[] gtProbs = new float[3];

    /**
     * Constructs a new {@code RecombSingleBaumLevel} instance from the
     * specified data.
     * @param samplerData the analysis data
     * @throws NullPointerException if {@code samplerData == null}
     */
    public RecombSingleBaumLevel(SamplerData samplerData) {
        this.samplerData = samplerData;
        this.dag = samplerData.rdag().dag();
        this.gl = samplerData.gl();
    }

    /**
     * Sets the Baum forward algorithm values for this level of the HMM
     * and records the child node pair values in the specified
     * {@code nodes} parameter. When the method call returns, the {@code nodes}
     * parameter will be reset to the child node pair values for this level of
     * the HMM.
     *
     * @param nodes child node pair values at the previous level of the HMM
     * @param permittedStates the permitted diploid model states
     * @param marker the level of the HMM at which the Baum forward algorithm
     * values will be computed
     * @param sample a sample index
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.dag().nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.gl().nSamples()}
     * @throws IndexOutOfBoundsException if either node in any node pair with
     * non-zero value is not a valid parent node at the specified level of the
     * HMM
     * @throws NullPointerException if
     * {@code nodes == null || permittedStates == null}
     */
    public void setForwardValues(RecombSingleNodes nodes,
            DiploidStates permittedStates, int marker, int sample) {
        this.marker = marker;
        this.sample = sample;
        this.nGenotypes = gl.marker(marker).nGenotypes();
        this.size = 0;
        this.fwdValueSum = 0.0f;
        this.bwdValueSum = 0.0f;
        initializeGtProbs(); // initialized here due to gtProbs() contract
        setStates(nodes, permittedStates);
        setChildNodes(nodes);
    }

    private void initializeGtProbs() {
        if (gtProbs.length < nGenotypes) {
            int newLength = Math.max(nGenotypes, (3*gtProbs.length/2 + 1));
            gtProbs = new float[newLength];
        }
        else {
            Arrays.fill(gtProbs, 0, nGenotypes, 0f);

        }
    }

    private void setStates(RecombSingleNodes nodes,
            DiploidStates permittedStates) {
        float valueSum = 0.0f;
        permittedStates.setMarker(marker);
        while (permittedStates.hasNext()) {
            permittedStates.next();
            int edge1 = permittedStates.edge1();
            int edge2 = permittedStates.edge2();
            float fwdValue = fwdValue(edge1, edge2, nodes);
            if (fwdValue > 0.0) {
                if (size == capacity) {
                    ensureCapacity(size+1);
                }
                edges1[size] = edge1;
                edges2[size] = edge2;
                fwdValues[size++] = fwdValue;
                valueSum += fwdValue;
            }
        }
        if (valueSum <= 0f) {
            throw new IllegalStateException(String.valueOf(valueSum));
        }
        for (int k=0; k<size; ++k) {
            this.fwdValues[k] /= valueSum;
        }
        fwdValueSum = valueSum;
    }

    private float fwdValue(int edge1, int edge2, RecombSingleNodes nodes) {
        float fwdValue = 0.0f;
        int symbol1 = dag.symbol(marker, edge1);
        int symbol2 = dag.symbol(marker, edge2);
        float emProb = gl.gl(marker, sample, symbol1, symbol2);

        if (emProb > 0.0) {
            float pRecom = samplerData.pRecomb(marker);
            float rec0 = (1-pRecom)*(1-pRecom);
            float rec1 = pRecom*(1-pRecom);
            float rec2 = pRecom*pRecom;
            float ep1 = dag.condEdgeProb(marker, edge1);
            float ep2 = dag.condEdgeProb(marker, edge2);
            float ep = ep1*ep2;

            int pn1 = dag.parentNode(marker, edge1);
            int pn2 = dag.parentNode(marker, edge2);
            float pnp1 = dag.parentProb(marker, pn1);
            float pnp2 = dag.parentProb(marker, pn2);

            fwdValue = rec0*emProb*ep*nodes.value(pn1, pn2);
            fwdValue += rec1*emProb*ep*pnp2*nodes.sumNode1Value(pn1);
            fwdValue += rec1*emProb*ep*pnp1*nodes.sumNode2Value(pn2);
            fwdValue += rec2*emProb*ep*pnp1*pnp2*nodes.sumValue();
            if (fwdValue<MIN_VALUE) {
                fwdValue = MIN_VALUE;
            }
        }
        return fwdValue;
    }

    /**
     * Stores the Baum forward algorithm child node pair values for this
     * level of the HMM in the specified {@code SingleNodes} object.
     *
     * @param nodes the node pair values that will be set
     *
     * @throws NullPointerException if {@code nodes == null}
     */
    public void setChildNodes(RecombSingleNodes nodes) {
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
     * @param nodes the node pair values to be initialized
     * @throws NullPointerException if {@code nodes == null}
     */
    public void setInitialBackwardValues(RecombSingleNodes nodes) {
        float bwdValue = 1.0f/size;
        bwdValueSum = size;
        nodes.clear();
        for (int j=0; j<size; ++j) {
            bwdValues[j] = bwdValue;
            int gtIndex = BasicGL.genotype(symbol1(j), symbol2(j));
            gtProbs[gtIndex] += fwdValues[j];

            float nextBaseBwdValue =  nextBaseBwdValue(edges1[j], edges2[j],
                    bwdValues[j]);
            int pn1 = dag.parentNode(marker, edges1[j]);
            int pn2 = dag.parentNode(marker, edges2[j]);
            nodes.sumUpdate(pn1, pn2, nextBaseBwdValue);
        }
    }

    /**
     * Sets the Baum backward algorithm values for this level of the HMM
     * and stores the parent node pair values in the specified
     * {@code nodes} parameter.  When the method call returns, this
     * {@code nodes} parameter will be reset to the parent
     * node pair values for this level of the HMM.
     *
     * @param nodes parent node pair values at the next level of HMM
     *
     * @throws IndexOutOfBoundsException if either node in any node pair with
     * non-zero value is not a valid child node at this level of the HMM
     * @throws NullPointerException if {@code nodes == null}
     */
    public void setBackwardValues(RecombSingleNodes nodes) {
        for (int j=0; j<size; ++j) {
            int child1 = dag.childNode(marker, edges1[j]);
            int child2 = dag.childNode(marker, edges2[j]);
            float childProb1 = dag.parentProb(marker+1, child1);
            float childProb2 = dag.parentProb(marker+1, child2);

            float pRecom = samplerData.pRecomb(marker+1);
            float rec0 = (1-pRecom)*(1-pRecom);
            float rec1 = pRecom*(1-pRecom);
            float rec2 = pRecom*pRecom;
            float bwdValue = rec0*nodes.value(child1, child2)
                    / (childProb1*childProb2);
            bwdValue += rec1*nodes.sumNode1Value(child1)/childProb1;
            bwdValue += rec1*nodes.sumNode2Value(child2)/childProb2;
            bwdValue += rec2*nodes.sumValue();
            bwdValues[j] = bwdValue;
            bwdValueSum += bwdValue;
        }
        nodes.clear();
        float gtProbsSum = 0f;
        for (int j=0; j<size; ++j) {
            bwdValues[j] /= bwdValueSum;
            float stateProb = (fwdValues[j] * bwdValues[j]);
            int gtIndex = BasicGL.genotype(symbol1(j), symbol2(j));
            // gtProbs assumed to be initialized in setForwardValues() method
            gtProbs[gtIndex] += stateProb;
            gtProbsSum += stateProb;

            float nextBaseBwdValue =  nextBaseBwdValue(edges1[j], edges2[j],
                    bwdValues[j]);
            if (nextBaseBwdValue>0f) {
                int pn1 = dag.parentNode(marker, edges1[j]);
                int pn2 = dag.parentNode(marker, edges2[j]);
                nodes.sumUpdate(pn1, pn2, nextBaseBwdValue);
            }
        }
        for (int j=0; j<nGenotypes; ++j) {
            gtProbs[j] /= gtProbsSum;
        }
    }

    private float nextBaseBwdValue(int edge1, int edge2, float lastBwdValue) {
        float ep1 = dag.edgeProb(marker, edge1);
        float ep2 = dag.edgeProb(marker, edge2);
        int symbol1 = dag.symbol(marker, edge1);
        int symbol2 = dag.symbol(marker, edge2);
        float emProb = gl.gl(marker, sample, symbol1, symbol2);
        float value = emProb*lastBwdValue*ep1*ep2;
        if (value<MIN_VALUE && lastBwdValue>0.0) {
            value = MIN_VALUE;
        }
        return value;
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
     * Returns the specified posterior genotype probability.  Returns 0
     * if the Baum backward probabilities have not been set.
     * @param gt a genotype index
     * @return the specified posterior genotype probability
     * @throws IndexOutOfBoundsException if
     * {@code gt < 0 || gt >= this.nGenotypes()}
     */
    public float gprobs(int gt) {
        if (gt >= nGenotypes) {
            String s = "gt=" + gt + " >= nGenotypes()=" + nGenotypes;
            throw new IllegalArgumentException(s);
        }
        return gtProbs[gt];
    }

    /**
     * Returns the current capacity of this level.
     * @return the current capacity of this level
     */
    public int capacity() {
        return edges1.length;
    }

    /**
     * Resets the size of this level to 0 and resets the capacity of this
     * level to the specified value.
     *
     * @param newCapacity the new capacity
     * @throws IllegalArgumentException if {@code newCapacity < 0}
     */
    public void reset(int newCapacity) {
        if (newCapacity<0) {
            throw new IllegalArgumentException(String.valueOf(newCapacity));
        }
        size = 0;
        capacity = newCapacity;
        edges1 = new int[newCapacity];
        edges2 = new int[newCapacity];
        fwdValues = new float[newCapacity];
        bwdValues = new float[newCapacity];
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
            String s = "state=" + state + " size()=" + size();
            throw new IndexOutOfBoundsException(s);
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
    public int edge1(int state) {
        checkIndex(state);
        return edges1[state];
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
    public int edge2(int state) {
        checkIndex(state);
        return edges2[state];
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
    public int parentNode1(int state) {
        checkIndex(state);
        return dag.parentNode(marker, edges1[state]);
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
    public int parentNode2(int state) {
        checkIndex(state);
        return dag.parentNode(marker, edges2[state]);
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
    public int childNode1(int state) {
        checkIndex(state);
        return dag.childNode(marker, edges1[state]);
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
    public int childNode2(int state) {
        checkIndex(state);
        return dag.childNode(marker, edges2[state]);
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
    public int symbol1(int state) {
        return dag.symbol(marker, edge1(state));
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
    public int symbol2(int state) {
        return dag.symbol(marker, edge2(state));
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
     * @param minCapacity the desired minimum state capacity
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
