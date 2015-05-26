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
import vcf.AL;

/**
 * Class {@code HapBaumLevel} computes forward and backward Baum values for a
 * haploid hidden Markov model (HMM) whose states are edges of a leveled
 * directed acyclic graph (DAG).
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class HapBaumLevel {

    private static final int INITIAL_CAPACITY=100;
    private static final double MIN_VALUE=100*Double.MIN_VALUE;
    private final Dag dag;
    private final AL al;

    private int marker=-1;
    private int hap=-1;
    private int size=0;

    private int capacity=INITIAL_CAPACITY;
    private int[] edges=new int[INITIAL_CAPACITY];
    private double[] fwdValues=new double[INITIAL_CAPACITY];
    private double[] bwdValues=new double[INITIAL_CAPACITY];
    private double fwdValueSum=0.0;
    private double bwdValueSum=0.0;

    private int nAlleles=0;
    private double[] alProbs=new double[3];

    /**
     * Constructs a new {@code HapBaumLevel} instance.
     *
     * @param dag the directed acyclic graph that the determines transition
     * probabilities.
     * @param al the emission probabilities
     * @throws IllegalArgumentException if
     * {@code dag.markers().equals(al.markers())==false}
     * @throws NullPointerException if {@code dag==null || al==null}
     */
    public HapBaumLevel(Dag dag, AL al) {
        if (dag.markers().equals(al.markers())==false) {
            throw new IllegalArgumentException("marker inconsistency");
        }
        this.dag=dag;
        this.al=al;
    }

    /**
     * Initializes the node values for the Baum forward algorithm.
     *
     * @param nodes the node values to be initialized.
     */
    public static void initializeNodes(HapNodes nodes) {
        nodes.clear();
        nodes.sumUpdate(0, 1.0);
    }

    /**
     * Sets the Baum forward algorithm values for this level of the HMM and
     * records the child node values in the specified {@code nodes} parameter.
     *
     * @param nodes child node values at the previous level of HMM. When the
     * method call returns, this parameter will be reset to the child node
     * values for this level of the HMM.
     * @param marker the level of the HMM at which the Baum forward algorithm
     * values will be computed.
     * @param haplotype a haplotype index.
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker<0 || marker>=this.dag().nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code haplotype<0 || haplotype>=2*this.al().nSamples()}
     * @throws IndexOutOfBoundsException if any node with non-zero value
     * is not a valid parent node at the specified level of the HMM
     * @throws NullPointerException if {@code nodes==null}
     */
    public void setForwardValues(HapNodes nodes, int marker, int haplotype) {
        this.marker=marker;
        this.hap=haplotype;
        this.nAlleles=al.marker(marker).nAlleles();
        this.size=0;
        this.fwdValueSum=0.0;
        this.bwdValueSum=0.0;
        setStates(nodes);
        setChildNodes(nodes);
    }

    private void setStates(HapNodes nodes) {
        double valueSum=0.0;
        for (int j=0, n=nodes.size(); j<n; ++j) {
            int node=nodes.enumNode(j);
            for (int k=0, m=dag.nOutEdges(marker, node); k<m; ++k) {
                int edge=dag.outEdge(marker, node, k);
                byte symbol=dag.symbol(marker, edge);
                float ep=al.al(marker, hap, symbol);
                if (ep>0.0f) {
                    if (size==capacity) {
                        ensureCapacity(size+1);
                    }
                    edges[size]=edge;
                    double tp=dag.condEdgeProb(marker, edge);
                    double fwdValue=ep*nodes.enumValue(j)*tp;
                    if (fwdValue<MIN_VALUE) {
                        assert nodes.enumValue(j)>0.0;
                        fwdValue=MIN_VALUE;
                    }
                    fwdValues[size++]=fwdValue;
                    valueSum+=fwdValue;
                }
            }
        }
        assert valueSum>0.0 ^ size==0;
        for (int k=0; k<size; ++k) {
            this.fwdValues[k]/=valueSum;
        }
        fwdValueSum=valueSum;
    }

    /**
     * Stores the Baum forward algorithm child node values for this
     * level of the HMM in the specified {@code HapNodes} object.
     *
     * @param nodes the node values that will be set.
     *
     * @throws NullPointerException if {@code nodes==null}
     */
    public void setChildNodes(HapNodes nodes) {
        nodes.clear();
        for (int k=0; k<size; ++k) {
            int node = dag.childNode(marker, edges[k]);
            nodes.sumUpdate(node, fwdValues[k]);
        }
    }

    /**
     * Initializes the node values for the Baum backward algorithm.
     *
     * @param nodes the node values to be initialized.
     * @throws NullPointerException if {@code nodes==null}
     */
    public void setInitialBackwardValues(HapNodes nodes) {
        nodes.clear();
        for (int j=0; j<size; ++j) {
            int node = dag.childNode(marker, edges[j]);
            nodes.maxUpdate(node, 1.0);
        }
        setBackwardValues(nodes);
    }

    /**
     * Sets the Baum backward algorithm values for this level of the HMM
     * and stores the parent node values in the specified {@code nodes}
     * parameter.
     *
     * @param nodes parent node values at the next level of HMM.  When
     * the method call returns, this parameter will be reset to the parent
     * node values for this level of the HMM.
     *
     * @throws IndexOutOfBoundsException if any node with non-zero value is
     * not a valid child node at the {@code this.marker()} level of the HMM
     * @throws NullPointerException if {@code nodes==null}
     */
    public void setBackwardValues(HapNodes nodes) {
        bwdValueSum=0.0;
        if (alProbs.length<nAlleles) {
            int newLength=Math.max(nAlleles, (3*alProbs.length/2+1));
            alProbs=new double[newLength];
        }
        Arrays.fill(alProbs, 0, nAlleles, 0.0);
        double alProbsSum=0.0;
        for (int j=0; j<size; ++j) {
            int node=dag.childNode(marker, edges[j]);
            double backwardValue=nodes.value(node);
            bwdValues[j]=backwardValue;
            bwdValueSum+=backwardValue;
        }
        nodes.clear();
        for (int j=0; j<size; ++j) {
            bwdValues[j]/=bwdValueSum;
            int edge=edges[j];
            byte symbol=symbol(j);
            int node=dag.parentNode(marker, edge);
            double tp=dag.condEdgeProb(marker, edge);

            double stateProb=fwdValues[j]*bwdValues[j];
            alProbs[symbol]+=stateProb;
            alProbsSum+=stateProb;

            double bwdValue=bwdValues[j]*tp*al.al(marker, hap, symbol);
            if (bwdValue<MIN_VALUE&&bwdValues[j]>0.0) {
                bwdValue=MIN_VALUE;
            }
            nodes.sumUpdate(node, bwdValue);
        }
        for (int j=0; j<nAlleles; ++j) {
            alProbs[j]/=alProbsSum;
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
    public AL emissions() {
        return al;
    }

    /**
     * Return the level of the HMM.
     * @return the level of the HMM.
     */
    public int marker() {
        return marker;
    }

    /**
     * Return the number of possible alleles at this level of the HMM.
     * @return the number of possible alleles at this level of the HMM.
     */
    public int nAlleles() {
        return nAlleles;
    }

    /**
     * Returns the specified posterior allele probability.  Returns 0
     * if the Baum backward values have not been set.
     * @param allele an allele index.
     * @return the specified posterior genotype probability.
     * @throws IndexOutOfBoundsException if
     * {@code allele<0 || allele>=this.nAlleles()}
     */
    public double alProbs(int allele) {
        if (allele >= nAlleles) {
            throw new IllegalArgumentException(String.valueOf(allele));
        }
        return alProbs[allele];
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
        if (state>=size) {
            throw new IndexOutOfBoundsException(String.valueOf(size));
        }
    }

    /**
     * Returns the DAG level edge index for the specified HMM state with
     * nonzero forward probability.
     * @param state an index of a HMM state at this level with nonzero
     * forward probability.
     * @return the DAG level edge index for the specified HMM state with
     * nonzero forward probability.
     *
     * @throws IndexOutOfBoundsException if
     * {@code state<0 || state>=this.size()}
     */
    public int edge(int state) {
        checkIndex(state);
        return edges[state];
    }

    /**
     * Returns the DAG level parent node index for the parent node of the
     * specified HMM state with nonzero forward probability.
     *
     * @param state an index of a HMM state at this level with nonzero
     * forward probability.
     * @return the DAG level parent node index for the parent node of the
     * specified HMM state with nonzero forward probability.
     *
     * @throws IndexOutOfBoundsException if
     * {@code state<0 || state>=this.size()}
     */
    public int parentNode(int state) {
        checkIndex(state);
        return dag.parentNode(marker, edges[state]);
    }

    /**
     * Returns the DAG level child node index for the child node of the
     * specified HMM state with nonzero forward probability.
     *
     * @param state an index of a HMM state at this level with nonzero
     * forward probability.
     * @return the DAG level child node index for the child node of the
     * specified HMM state with nonzero forward probability.
     *
     * @throws IndexOutOfBoundsException if
     * {@code state<0 || state>=this.size()}
     */
    public int childNode(int state) {
        checkIndex(state);
        return dag.childNode(marker, edges[state]);
    }

    /**
     * Returns the symbol of the specified HMM state with nonzero forward
     * probability.
     *
     * @param state an index of a HMM state at this level with nonzero
     * forward probability.
     * @return the symbol of the specified HMM state with nonzero forward
     * probability.
     *
     * @throws IndexOutOfBoundsException if
     * {@code state<0 || state>=this.size()}
     */
    public byte symbol(int state) {
        return dag.symbol(marker, edge(state));
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
     * Returns a string description of {@code this}. The exact details of the
     * description are unspecified and subject to change.
     *
     * @return a string description of {@code this}.
     */
    @Override
    public String toString() {
        String space=" ";
        String sep=" | ";
        StringBuilder sb=new StringBuilder(100);
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
            sb.append((int) edge(j));
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
        if (minCapacity>capacity) {
            capacity=(capacity*3)/2+1;
            if (capacity<minCapacity) {
                capacity=minCapacity;
            }
            edges=Arrays.copyOf(edges, capacity);
            fwdValues=Arrays.copyOf(fwdValues, capacity);
            bwdValues=Arrays.copyOf(bwdValues, capacity);
        }
    }
}
