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
package ibd;

import dag.Dag;
import haplotype.HapPairs;
import sample.DuoBaumLevel;
import sample.DuoNodes;
import sample.HapBaumLevel;
import sample.HapNodes;
import sample.SingleBaumLevel;
import sample.SingleNodes;
import vcf.GL;
import vcf.HbdAL;

/**
 * <p>Class {@code IbdBaum} estimates LOD scores for an IBD versus a non-IBD
 * model, and it estimates LOD scores for an HBD versus a non-HBD model.
 * </p>
 * <p>Instances of class {@code IbdBaum} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class IbdBaum {

    private static final double log_e_10 = Math.log(10.0);
    private final Dag dag;
    private final GL gl;
    private final int nMarkers;

    private final SingleNodes fwdNodesA;
    private final SingleNodes fwdNodesB;
    private final HapNodes fwdNodesHbd;
    private final DuoNodes fwdNodesIbd;

    private final SingleBaumLevel scratchSingleLevel;
    private final HapBaumLevel scratchHapLevel;
    private final DuoBaumLevel scratchDuoLevel;

    /**
     * Creates a new {@code IbdBaum} instance from the specified data.
     *
     * @param dag the directed acyclic graph that determines the
     * transition probabilities
     * @param gl the HMM emission probabilities
     *
     * @throws IllegalArgumentException
     * if {@code dag.markers().equals(gl.markers()) == false}
     * @throws NullPointerException if {@code dag == null || gl == null}
     */
    public IbdBaum(Dag dag, GL gl) {
        if (dag.markers().equals(gl.markers())==false) {
            throw new IllegalArgumentException("inconsistent markers");
        }
        this.dag = dag;
        this.gl = gl;
        this.nMarkers = dag.nLevels();

        this.fwdNodesHbd = new HapNodes();
        this.fwdNodesA = new SingleNodes();
        this.fwdNodesB = new SingleNodes();
        this.fwdNodesIbd = new DuoNodes();

        this.scratchSingleLevel = new SingleBaumLevel(dag, gl);
        this.scratchHapLevel = new HapBaumLevel(dag, new HbdAL(gl));
        this.scratchDuoLevel = new DuoBaumLevel(dag, gl);
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
     * Returns the HMM emission probabilities.
     * @return the HMM emission probabilities
     */
    public GL gl() {
        return gl;
    }

    /**
     * Returns the homozygosity-by-descent (HBD) LOD score.
     * @param sample the sample index
     * @param start the start marker index (inclusive)
     * @param end the end marker index (exclusive)
     * @return the HBD LOD score
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.gl().nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code start < 0 || start > end || end > this.dag.nMarkers()}
     */
    public double hbdLod(int sample, int start, int end) {
        checkStartAndEnd(start, end, nMarkers);
        if (start==end) {
            return 0.0f;
        }
        setInitialNodes(dag, start, fwdNodesHbd);
        setInitialNodes(dag, start, fwdNodesA);
        double altHbdLogLike = logLikelihood(fwdNodesHbd, scratchHapLevel,
                sample, start, end);
        double nullHbdLogLike = logLikelihood(fwdNodesA, scratchSingleLevel,
                sample, start, end);
        return lod(altHbdLogLike - nullHbdLogLike);
    }

    /**
     * Returns the identity-by-descent (IBD) LOD score.
     * @param sampleA the first sample index
     * @param sampleB the second sample index
     * @param start the start marker index (inclusive)
     * @param end the end marker index (exclusive)
     * @return the IBD LOD score
     * @throws IndexOutOfBoundsException if
     * {@code sampleA < 0 || sampleA >= this.gl().nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code sampleB < 0 || sampleB >= this.gl().nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code start < 0 || start > end || end > this.dag.nMarkers()}
     */
    public double ibdLod(int sampleA, int sampleB, int start, int end) {
        checkStartAndEnd(start, end, nMarkers);
        if (start==end) {
            return 0.0f;
        }
        setInitialNodes(dag, start, fwdNodesA);
        setInitialNodes(dag, start, fwdNodesB);
        setInitialNodes(dag, start, fwdNodesIbd);

        double nullLogLike = 0.0;
        nullLogLike += logLikelihood(fwdNodesA, scratchSingleLevel, sampleA,
                start, end);
        nullLogLike += logLikelihood(fwdNodesB, scratchSingleLevel, sampleB,
                start, end);

        double altLogLike = logLikelihood(fwdNodesIbd, scratchDuoLevel,
                sampleA, sampleB, start, end);
        return lod(altLogLike - nullLogLike);
    }

    private static void checkStartAndEnd(int start, int end, int nMarkers) {
        if (start<0 || start>end || end>nMarkers) {
            String s = "start=" + start + " end=" + end + " nMarkers="
                    + nMarkers;
            throw new IllegalArgumentException(s);
        }
    }

    private static void setInitialNodes(Dag dag, int marker, HapNodes nodes) {
        nodes.clear();
        int nNodes = dag.nParentNodes(marker);
        for (int node=0; node<nNodes; ++node) {
            float p = dag.parentProb(marker, node);
            nodes.sumUpdate(node, p);
        }
    }

    private static void setInitialNodes(Dag dag, int marker, SingleNodes nodes) {
        nodes.clear();
        int n = dag.nParentNodes(marker);
        for (int n1=0; n1<n; ++n1) {
            float p1 = dag.parentProb(marker, n1);
            for (int n2=0; n2<n; ++n2) {
                float p2 = dag.parentProb(marker, n2);
                nodes.sumUpdate(n1, n2, (p1*p2));
            }
        }
    }

    private static void setInitialNodes(Dag dag, int marker, DuoNodes nodes) {
        nodes.clear();
        int n = dag.nParentNodes(marker);
        for (int n1=0; n1<n; ++n1) {
            float p1 = dag.parentProb(marker, n1);
            for (int n2=0; n2<n; ++n2) {
                float p2 = dag.parentProb(marker, n2);
                for (int n3=0; n3<n; ++n3) {
                    float p3 = dag.parentProb(marker, n3);
                    nodes.sumUpdate(n1, n2, n3, (p1*p2*p3));
                }
            }
        }
    }

    private static double logLikelihood(HapNodes nodes,
            HapBaumLevel level, int sample, int start, int end) {
        int hap = 2*sample;
        double sum = 0.0;
        for (int j=start; j<end; ++j) {
            level.setForwardValues(nodes, j, hap);
            sum += Math.log(level.forwardValuesSum());
        }
        return sum;
    }

    private static double logLikelihood(SingleNodes nodes,
            SingleBaumLevel level, int sample, int start, int end) {
        double sum = 0.0;
        for (int j=start; j<end; ++j) {
            level.setForwardValues(nodes, j, sample);
            sum += Math.log(level.forwardValuesSum());
        }
        return sum;
    }

    private static double logLikelihood(DuoNodes nodes,
            DuoBaumLevel ibdLevel, int sampleA, int sampleB,
            int start, int end) {
        double sum = 0.0;
        for (int level=start; level<end; ++level) {
            ibdLevel.setForwardValues(nodes, level, sampleA, sampleB);
            sum += Math.log(ibdLevel.forwardValuesSum());
        }
        return sum;
    }

    private static double lod(double logLR) {
        if (Double.isNaN(logLR)) {
            return 0.0;
        }
        else {
            return (logLR / log_e_10);
        }
    }

    /**
     * Returns the estimated frequency of the haplotype segment on the LOD
     * {@code (-Math.log10)} scale.
     * @param hap a haplotype index
     * @param start the start marker index (inclusive)
     * @param end the end marker index (exclusive)
     * @param hapPairs the list of haplotype pairs
     * @param dag the directed acyclic graph that determines the HMM
     * transition probabilities
     * @return the estimated frequency of the haplotype segment on the LOD
     * {@code (-Math.log10)} scale
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap >= haps.nHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code start < 0 || start > end || end > dag.nMarkers()}
     * @throws NullPointerException if {@code dag == null || haps == null}
     */
    public static double freqLod(int hap, int start, int end, HapPairs hapPairs,
            Dag dag) {
        double minSumProbs = 1e-100;
        checkStartAndEnd(start, end, dag.nLevels());
        double sumProbs = 0.0;
        for (int node=0, n=dag.nParentNodes(start); node<n; ++node) {
            int lastNode = node;
            double p = dag.parentProb(start, node);
            for (int m=start; m<end && p>0.0; ++m) {
                int allele = hapPairs.allele(m, hap);
                int e = dag.outEdgeBySymbol(m, lastNode, allele);
                if (e == -1) {
                    p = 0.0;
                    break;
                }
                else {
                    p *= dag.condEdgeProb(m, e);
                    lastNode = dag.childNode(m, e);
                }
            }
            sumProbs += p;
        }
        if (sumProbs < minSumProbs) {
            sumProbs = minSumProbs;
        }
        return lod(-Math.log(sumProbs));
    }
}
