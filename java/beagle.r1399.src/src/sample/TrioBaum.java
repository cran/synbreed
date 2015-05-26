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
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import haplotype.HapPair;
import haplotype.BitHapPair;
import vcf.GL;

/**
 * Class {@code Trio Baum} implements the Baum forward and backward
 * algorithms for a hidden Markov model (HMM) of a parent-offspring trio's
 * genotype data.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class TrioBaum {

    private final Dag dag;
    private final GL gl;
    private final int nMarkers;
    private final int nCopies;
    private final long seed;
    private final Random random;

    private final int[] nodeA1;
    private final int[] nodeA2;
    private final int[] nodeB1;
    private final int[] nodeB2;
    private final double[] nodeValue;

    private final byte[][] allelesA1;
    private final byte[][] allelesA2;
    private final byte[][] allelesB1;
    private final byte[][] allelesB2;

    private final TrioBaumLevel[] levels;
    private final TrioNodes fwdNodes;
    private final TrioNodes bwdNodes;

    private int windowIndex = -9999;
    private int arrayIndex = -9999;

    /**
     * Creates a new {@code trioBaum} instance.
     *
     * @param dag the directed acyclic graph that determines the
     * transition probabilities.
     * @param gl the emission probabilities.
     * @param seed the initial random seed.
     * @param nCopies the number of haplotype pairs that will be sampled for
     * each individual in a parent-offspring duo.
     *
     * @throws IllegalArgumentException
     * if {@code nCopies<1 || dag.markers().equals(gl.markers())==false}
     * @throws NullPointerException if {@code dag==null || gl==null}
     */
    public TrioBaum(Dag dag, GL gl, long seed, int nCopies) {
        if (dag.markers().equals(gl.markers())==false) {
            throw new IllegalArgumentException("inconsistent markers");
        }
        if (nCopies < 1) {
            throw new IllegalArgumentException("nCopies<1: " + nCopies);
        }
        this.dag = dag;
        this.gl = gl;
        this.nMarkers = dag.nMarkers();
        this.nCopies = nCopies;
        this.seed = seed;
        this.random = new Random(seed);

        this.nodeA1 = new int[nCopies];
        this.nodeA2 = new int[nCopies];
        this.nodeB1 = new int[nCopies];
        this.nodeB2 = new int[nCopies];
        this.nodeValue = new double[nCopies];
        this.allelesA1 = new byte[nCopies][gl.nMarkers()];
        this.allelesA2 = new byte[nCopies][gl.nMarkers()];
        this.allelesB1 = new byte[nCopies][gl.nMarkers()];
        this.allelesB2 = new byte[nCopies][gl.nMarkers()];

        int size = (int) Math.ceil(Math.sqrt(1 + 8*dag.nMarkers())/2.0) + 1;
        this.levels = new TrioBaumLevel[size];
        for (int j=0; j<levels.length; ++j) {
            levels[j] = new TrioBaumLevel(dag, gl);
        }
        this.fwdNodes = new TrioNodes();
        this.bwdNodes = new TrioNodes();
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
     * Returns the number of haplotype pairs that are sampled for each
     * individual in a parent-offspring trio.
     * @return the number of haplotype pairs that are sampled for each
     * individual in a parent-offspring trio.
     */
    public int nCopies() {
        return nCopies;
    }

    /**
     * Returns the initial random seed.
     * @return the initial random seed.
     */
    public long seed() {
        return seed;
    }

    /**
     * <p>Returns a list of {@code this.nCopies()} sampled haplotype pairs for
     * the specified father ({@code sampleA}), mother ({@code sampleB}), and
     * offspring ({@code sampleC}).
     * Haplotype pairs are sampled conditional on the HMM with transition
     * probabilities determined by {@code this.dag()} and
     * emission probabilities determined by {@code this.gl()}.
     * </p>
     * The contract for this method is unspecified if no father haplotype
     * pairs, no mother haplotype pairs, or no offspring haplotype pairs
     * are consistent with the HMM.
     *
     * @param sampleA the sample index of the father.
     * @param sampleB the sample index of the mother.
     * @param sampleC the sample index of the offspring.
     * @return a list of {@code this.nCopies()} sampled haplotype pairs for the
     * specified individuals.
     *
     * @throws IndexOutOfBoundsException if
     * {@code sampleA<0 || sampleA>=this.gl().nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code sampleB<0 || sampleB>=this.gl().nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code sampleC<0 || sampleC>=this.gl().nSamples()}
     */
    public List<HapPair> sample(int sampleA, int sampleB,
            int sampleC) {
        forwardAlgorithm(sampleA, sampleB, sampleC);
        initSampleAlleles(currentLevel(), sampleA, sampleB, sampleC);
        for (int j=nMarkers-2; j>=0; --j) {
            TrioBaumLevel level = previousLevel(sampleA, sampleB, sampleC);
            sampleAlleles(level, sampleA, sampleB, sampleC);
        }
        return hapList(sampleA, sampleB, sampleC);
    }

    /**
     * <p>Returns a list of {@code this.nCopies()} sampled haplotype pairs for
     * the specified father ({@code sampleA}), mother ({@code sampleB}), and
     * offspring ({@code sampleC}).
     * Haplotype pairs are sampled conditional on the HMM with transition
     * probabilities determined by {@code this.dag()} and
     * emission probabilities determined by {@code this.gl()}.
     * Posterior genotype probabilities for the father, mother, and offspring
     * are written to {@code gtProbsA}, {@code gtProbsB}, and
     * {@code gtProbsC} respectively.  The posterior probability of the
     * {@code j}-th genotype for the {@code k}-th marker is stored
     * in element {@code gl.markers().sumGenotypes(k) + j} of the
     * {@code gtProbsA}, {@code gtProbsB}, and {@code gtProbsC}
     * arrays.
     * </p>
     * The contract for this method is unspecified if no father haplotype
     * pairs, no mother haplotype pairs, or no offspring haplotype pairs
     * are consistent with the HMM.
     *
     * @param sampleA the sample index of the father.
     * @param sampleB the sample index of the mother.
     * @param sampleC the sample index of the offspring.
     * @param gtProbsA an array to which posterior genotype probabilities
     * for the father will be written.
     * @param gtProbsB an array to which posterior genotype probabilities
     * for the mother will be written.
     * @param gtProbsC an array to which posterior genotype probabilities
     * for the offspring will be written.
     * @return a list of {@code this.nCopies()} sampled haplotype pairs for the
     * specified individuals.
     *
     * @throws IndexOutOfBoundsException if
     * {@code sampleA<0 || sampleA>=this.gl().nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code sampleB<0 || sampleB>=this.gl().nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code sampleC<0 || sampleC>=this.gl().nSamples()}
     *
     * @throws IllegalArgumentException if
     * {@code gtProbsA.length!=this.gl().markers().sumGenotypes()}
     * @throws IllegalArgumentException
     * {@code gtProbsB.length!=this.gl().markers().sumGenotypes()}
     * @throws IllegalArgumentException
     * {@code gtProbsC.length!=this.gl().markers().sumGenotypes()}
     * @throws NullPointerException if
     * {@code gtProbsA==null || gtProbsB==null || gtProbsC==null}
     */
    public List<HapPair> sample(int sampleA, int sampleB, int sampleC,
            double[] gtProbsA, double[] gtProbsB, double[] gtProbsC) {
        checkGtProbs(gtProbsA, gtProbsB, gtProbsC);
        forwardAlgorithm(sampleA, sampleB, sampleC);

        initSampleAlleles(currentLevel(), sampleA, sampleB, sampleC);
        currentLevel().setInitialBackwardValues(bwdNodes);
        setGtProbs(currentLevel(), gtProbsA, gtProbsB, gtProbsC);

        for (int j=nMarkers-2; j>=0; --j) {
            TrioBaumLevel level = previousLevel(sampleA, sampleB, sampleC);
            sampleAlleles(level, sampleA, sampleB, sampleC);
            level.setBackwardValues(bwdNodes);
            setGtProbs(level, gtProbsA, gtProbsB, gtProbsC);
        }
        return hapList(sampleA, sampleB, sampleC);
    }

    private void checkGtProbs(double[] gtProbsA, double[] gtProbsB,
            double[] gtProbsC) {
        int n = gl.markers().sumGenotypes();
        if (gtProbsA.length!=n || gtProbsB.length!=n || gtProbsC.length!=n) {
            String s = "array length error";
            throw new IllegalArgumentException(s);
        }
    }

    private void setGtProbs(TrioBaumLevel level, double[] gtProbsA,
            double[] gtProbsB, double[] gtProbsC) {
        int m = level.marker();
        int nGenotypes = gl.marker(m).nGenotypes();
        int base = gl.markers().sumGenotypes(m);
        for (int j=0; j<nGenotypes; ++j) {
            gtProbsA[base + j] = level.gtProbsA(j);
            gtProbsB[base + j] = level.gtProbsB(j);
            gtProbsC[base + j] = level.gtProbsC(j);
        }
    }

    private List<HapPair> hapList(int sampleA, int sampleB, int sampleC) {
        List<HapPair> hapList = new ArrayList<>(2*nCopies);
        for (int copy=0; copy<nCopies; ++copy) {
            HapPair hapsA = new BitHapPair(gl.markers(),
                    gl.samples().idIndex(sampleA), allelesA1[copy], allelesA2[copy]);
            HapPair hapsB = new BitHapPair(gl.markers(),
                    gl.samples().idIndex(sampleB), allelesB1[copy], allelesB2[copy]);
            HapPair hapsC = new BitHapPair(gl.markers(),
                    gl.samples().idIndex(sampleC), allelesA1[copy], allelesB1[copy]);
            hapList.add(hapsA);
            hapList.add(hapsB);
            hapList.add(hapsC);
        }
        return hapList;
    }

    private void initSampleAlleles(TrioBaumLevel level, int sampleA,
            int sampleB, int sampleC) {
        int m = level.marker();
        for (int copy=0; copy<nCopies; ++copy) {
            int state = initialRandomState(level);
            nodeA1[copy] = level.parentNodeA1(state);
            nodeA2[copy] = level.parentNodeA2(state);
            nodeB1[copy] = level.parentNodeB1(state);
            nodeB2[copy] = level.parentNodeB2(state);
            nodeValue[copy] =  parentSum(level, sampleA, sampleB, sampleC, state);
            allelesA1[copy][m] = level.symbolA1(state);
            allelesA2[copy][m] = level.symbolA2(state);
            allelesB1[copy][m] = level.symbolB1(state);
            allelesB2[copy][m] = level.symbolB2(state);
        }
    }

    private int initialRandomState(TrioBaumLevel level) {
        double d = random.nextDouble();
        double sum = 0.0;
        for (int j=0, n=level.size(); j<n; ++j) {
            sum += level.forwardValue(j);
            if (d <= sum) {
                return j;
            }
        }
        assert false;
        return level.size()-1; // error in finite bit arithmetic encountered
    }

    private double parentSum(TrioBaumLevel level, int sampleA, int sampleB,
            int sampleC, int state) {
        int marker = level.marker();
        double fwdValue = level.forwardValuesSum()*level.forwardValue(state);
        int edgeA1 = level.edgeA1(state);
        int edgeA2 = level.edgeA2(state);
        int edgeB1 = level.edgeB1(state);
        int edgeB2 = level.edgeB2(state);
        double tpA1 = dag.condEdgeProb(marker, edgeA1);
        double tpA2 = dag.condEdgeProb(marker, edgeA2);
        double tpB1 = dag.condEdgeProb(marker, edgeB1);
        double tpB2 = dag.condEdgeProb(marker, edgeB2);
        byte symbolA1 = dag.symbol(marker, edgeA1);
        byte symbolA2 = dag.symbol(marker, edgeA2);
        byte symbolB1 = dag.symbol(marker, edgeB1);
        byte symbolB2 = dag.symbol(marker, edgeB2);
        double epA = gl.gl(marker, sampleA, symbolA1, symbolA2);
        double epB = gl.gl(marker, sampleB, symbolB1, symbolB2);
        double epC = gl.gl(marker, sampleC, symbolA1, symbolB1);
        return fwdValue / ( (epA*epB*epC) * (tpA1*tpA2*tpB1*tpB2) );
    }

    private void sampleAlleles(TrioBaumLevel level, int sampleA, int sampleB,
            int sampleC) {
        int m = level.marker();
        for (int copy=0; copy<nCopies; ++copy) {
            int state = randomPreviousState(level, nodeA1[copy], nodeA2[copy],
                    nodeB1[copy], nodeB2[copy], nodeValue[copy]);
            nodeA1[copy] = level.parentNodeA1(state);
            nodeA2[copy] = level.parentNodeA2(state);
            nodeB1[copy] = level.parentNodeB1(state);
            nodeB2[copy] = level.parentNodeB2(state);
            nodeValue[copy] =  parentSum(level, sampleA, sampleB, sampleC, state);
            allelesA1[copy][m] = level.symbolA1(state);
            allelesA2[copy][m] = level.symbolA2(state);
            allelesB1[copy][m] = level.symbolB1(state);
            allelesB2[copy][m] = level.symbolB2(state);
        }
    }

    private int randomPreviousState(TrioBaumLevel level, int nodeA1,
            int nodeA2, int nodeB1, int nodeB2, double nodeValue) {
        double d = random.nextDouble() * nodeValue;
        double sum = 0.0;
        for (int j=0, n=level.size(); j<n; ++j) {
            if ( nodeA1==level.childNodeA1(j)
                    && nodeA2==level.childNodeA2(j)
                    && nodeB1==level.childNodeB1(j)
                    && nodeB2==level.childNodeB2(j) ) {
                sum += level.forwardValue(j);
                if (d <= sum) {
                    return j;
                }
            }
        }
        return level.size()-1; // error in finite bit arithmetic encountered
    }

    private TrioBaumLevel nextLevel() {
        ++arrayIndex;
        if (arrayIndex == levels.length) {
            ++windowIndex;
            arrayIndex = windowIndex;
        }
        return levels[arrayIndex];
    }

    private TrioBaumLevel currentLevel() {
        return levels[arrayIndex];
    }

    private TrioBaumLevel previousLevel(int sampleA, int sampleB, int sampleC) {
        if (arrayIndex == windowIndex) {
            --windowIndex;
            arrayIndex = windowIndex;
            levels[arrayIndex].setChildNodes(fwdNodes);
            int startLevel = levels[windowIndex].marker() + 1;
            int endLevel = startLevel + (levels.length - (windowIndex + 1) );
            for (int marker=startLevel; marker<endLevel; ++marker) {
                nextLevel().setForwardValues(fwdNodes, marker,
                        sampleA, sampleB, sampleC);
            }
            return currentLevel();
        }
        else {
            return levels[--arrayIndex];
        }
    }

    private void forwardAlgorithm(int sampleA, int sampleB, int sampleC) {
        TrioBaumLevel.initializeNodes(fwdNodes);
        this.windowIndex = -1;
        this.arrayIndex = levels.length - 1;
        for (int marker=0; marker<nMarkers; ++marker) {
            nextLevel().setForwardValues(fwdNodes, marker, sampleA, sampleB,
                    sampleC);
        }
    }
}
