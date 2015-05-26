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
 * Class {@code DuoBaum} implements the Baum forward and backward
 * algorithms for a hidden Markov model (HMM) of a parent-offspring duo's
 * genotype data.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class DuoBaum {

    private final Dag dag;
    private final GL gl;
    private final int nMarkers;
    private final int nCopies;
    private final long seed;
    private final Random random;

    private final int[] nodeAB1;
    private final int[] nodeA2;
    private final int[] nodeB2;
    private final double[] nodeValue;

    private final byte[][] allelesAB1;
    private final byte[][] allelesA2;
    private final byte[][] allelesB2;

    private final DuoBaumLevel[] levels;
    private final DuoNodes fwdNodes;
    private final DuoNodes bwdNodes;

    private int windowIndex = -9999;
    private int arrayIndex = -9999;

    /**
     * Creates a new {@code DuoBaum} instance.
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
    public DuoBaum(Dag dag, GL gl, long seed, int nCopies) {
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

        this.nodeAB1 = new int[nCopies];
        this.nodeA2 = new int[nCopies];
        this.nodeB2 = new int[nCopies];
        this.nodeValue = new double[nCopies];
        this.allelesAB1 = new byte[nCopies][gl.nMarkers()];
        this.allelesA2 = new byte[nCopies][gl.nMarkers()];
        this.allelesB2 = new byte[nCopies][gl.nMarkers()];

        int size = (int) Math.ceil(Math.sqrt(1 + 8*dag.nMarkers())/2.0) + 1;
        this.levels = new DuoBaumLevel[size];
        for (int j=0; j<levels.length; ++j) {
            levels[j] = new DuoBaumLevel(dag, gl);
        }
        this.fwdNodes = new DuoNodes();
        this.bwdNodes = new DuoNodes();
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
     * individual in a parent-offspring duo.
     * @return the number of haplotype pairs that are sampled for each
     * individual in a parent-offspring duo.
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
     * the specified parent ({@code sampleA}) and offspring ({@code sampleB}).
     * Haplotype pairs are sampled conditional on the HMM with transition
     * probabilities determined by {@code this.dag()} and
     * emission probabilities determined by {@code this.gl()}.
     * </p>
     * The contract for this method is unspecified if no parent haplotype
     * pairs or no offspring haplotype pairs are consistent with the HMM.
     *
     * @param sampleA the sample index of the parent.
     * @param sampleB the sample index of the offspring.
     * @return a list of {@code this.nCopies()} sampled haplotype pairs for the
     * specified individuals.
     *
     * @throws IndexOutOfBoundsException if
     * {@code sampleA<0 || sampleA>=this.gl().nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code sampleB<0 || sampleB>=this.gl().nSamples()}
     */
    public List<HapPair> sample(int sampleA, int sampleB) {
        forwardAlgorithm(sampleA, sampleB);
        initSampleAlleles(currentLevel(), sampleA, sampleB);
        for (int j=nMarkers-2; j>=0; --j) {
            DuoBaumLevel level = previousLevel(sampleA, sampleB);
            sampleAlleles(level, sampleA, sampleB);
        }
        return hapList(sampleA, sampleB);
    }

    /**
     * <p>Returns a list of {@code this.nCopies()} sampled haplotype pairs for
     * the specified parent ({@code sampleA}) and offspring ({@code sampleB}).
     * Haplotype pairs are sampled conditional on the HMM with transition
     * probabilities determined by {@code this.dag()} and
     * emission probabilities determined by {@code this.gl()}.
     * Posterior genotype probabilities for the parent and offspring are
     * written to {@code gtProbsA} and {@code gtProbsB} respectively.
     * The posterior probability of the {@code j}-th genotype for
     * the {@code k}-th marker is stored at index
     * {@code gl.markers().sumGenotypes(k) + j} in the {@code gtProbsA}
     * and {@code gtProbsB} arrays.
     * </p>
     * The contract for this method is unspecified if no parent haplotype
     * pairs or no offspring haplotype pairs are consistent with the HMM.
     *
     * @param sampleA the sample index of the parent.
     * @param sampleB the sample index of the offspring.
     * @param gtProbsA an array to which posterior genotype probabilities
     * for the parent will be written.
     * @param gtProbsB an array to which posterior genotype probabilities
     * for the offspring will be written.
     * @return a list of {@code this.nCopies()} sampled haplotype pairs for the
     * specified individuals.
     *
     * @throws IndexOutOfBoundsException if
     * {@code sampleA<0 || sampleA>=this.gl().nSamples()}
     * @throws IndexOutOfBoundsException if
     * {@code sampleB<0 || sampleB>=this.gl().nSamples()}
     * @throws IllegalArgumentException if
     * {@code gtProbsA.length!=this.gl().markers().sumGenotypes()}
     * @throws IllegalArgumentException if
     * {@code gtProbsB.length!=this.gl().markers().sumGenotypes()}
     * @throws NullPointerException if
     * {@code gtProbsA==null || gtProbsB==null}
     */
    public List<HapPair> sample(int sampleA, int sampleB,
            double[] gtProbsA, double[] gtProbsB) {
        checkGprobs(gtProbsA, gtProbsB);
        forwardAlgorithm(sampleA, sampleB);

        initSampleAlleles(currentLevel(), sampleA, sampleB);
        currentLevel().setInitialBackwardValues(bwdNodes);
        setGprobs(currentLevel(), gtProbsA, gtProbsB);

        for (int j=nMarkers-2; j>=0; --j) {
            DuoBaumLevel level = previousLevel(sampleA, sampleB);
            sampleAlleles(level, sampleA, sampleB);
            level.setBackwardValues(bwdNodes);
            setGprobs(level, gtProbsA, gtProbsB);
        }
        return hapList(sampleA, sampleB);
    }

    private void checkGprobs(double[] gtProbsA, double[] gtProbsB) {
        int n = gl.markers().sumGenotypes();
        if (gtProbsA.length!=n || gtProbsB.length!=n) {
            String s = "gtProbs.length!=gl.markers().sumGenotypes()";
            throw new IllegalArgumentException(s);
        }
    }

    private void setGprobs(DuoBaumLevel level, double[] gtProbsA, double[] gtProbsB) {
        int m = level.marker();
        int nGenotypes = gl.marker(m).nGenotypes();
        int base = gl.markers().sumGenotypes(m);
        for (int j=0; j<nGenotypes; ++j) {
            gtProbsA[base + j] = level.gtProbsA(j);
            gtProbsB[base + j] = level.gtProbsB(j);
        }
    }

    private List<HapPair> hapList(int sampleA, int sampleB) {
        List<HapPair> hapList = new ArrayList<>(2*nCopies);
        for (int copy=0; copy<nCopies; ++copy) {
            HapPair hapsA = new BitHapPair(gl.markers(),
                    gl.samples().idIndex(sampleA), allelesAB1[copy], allelesA2[copy]);
            HapPair hapsB = new BitHapPair(gl.markers(),
                    gl.samples().idIndex(sampleB), allelesAB1[copy], allelesB2[copy]);
            hapList.add(hapsA);
            hapList.add(hapsB);
        }
        return hapList;
    }

    private void initSampleAlleles(DuoBaumLevel level, int sampleA, int sampleB) {
        int m = level.marker();
        for (int copy=0; copy<nCopies; ++copy) {
            int state = initialRandomState(level);
            nodeAB1[copy] = level.parentNodeAB1(state);
            nodeA2[copy] = level.parentNodeA2(state);
            nodeB2[copy] = level.parentNodeB2(state);
            nodeValue[copy] =  parentSum(level, sampleA, sampleB, state);
            allelesAB1[copy][m] = level.symbolAB1(state);
            allelesA2[copy][m] = level.symbolA2(state);
            allelesB2[copy][m] = level.symbolB2(state);
        }
    }

    private int initialRandomState(DuoBaumLevel level) {
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

    private double parentSum(DuoBaumLevel level, int sampleA, int sampleB,
            int state) {
        int marker = level.marker();
        double fwdValue = level.forwardValuesSum()*level.forwardValue(state);
        int edgeAB1 = level.edgeAB1(state);
        int edgeA2 = level.edgeA2(state);
        int edgeB2 = level.edgeB2(state);
        double tpAB1 = dag.condEdgeProb(marker, edgeAB1);
        double tpA2 = dag.condEdgeProb(marker, edgeA2);
        double tpB2 = dag.condEdgeProb(marker, edgeB2);
        byte symbolAB1 = dag.symbol(marker, edgeAB1);
        byte symbolA2 = dag.symbol(marker, edgeA2);
        byte symbolB2 = dag.symbol(marker, edgeB2);
        double epA = gl.gl(marker, sampleA, symbolAB1, symbolA2);
        double epB = gl.gl(marker, sampleB, symbolAB1, symbolB2);
        return fwdValue / ( (epA*epB) * (tpAB1*tpA2*tpB2) );
    }

    private void sampleAlleles(DuoBaumLevel level, int sampleA, int sampleB) {
        int m = level.marker();
        for (int copy=0; copy<nCopies; ++copy) {
            int state = randomPreviousState(level, nodeAB1[copy], nodeA2[copy],
                    nodeB2[copy], nodeValue[copy]);
            nodeAB1[copy] = level.parentNodeAB1(state);
            nodeA2[copy] = level.parentNodeA2(state);
            nodeB2[copy] = level.parentNodeB2(state);
            nodeValue[copy] =  parentSum(level, sampleA, sampleB, state);
            allelesAB1[copy][m] = level.symbolAB1(state);
            allelesA2[copy][m] = level.symbolA2(state);
            allelesB2[copy][m] = level.symbolB2(state);
        }
    }

    private int randomPreviousState(DuoBaumLevel level, int nodeAB1,
            int nodeA2, int nodeB2, double nodeValue) {
        double d = random.nextDouble() * nodeValue;
        double sum = 0.0;
        for (int j=0, n=level.size(); j<n; ++j) {
            if ( nodeAB1==level.childNodeAB1(j)
                    && nodeA2==level.childNodeA2(j)
                    && nodeB2==level.childNodeB2(j) ) {
                sum += level.forwardValue(j);
                if (d <= sum) {
                    return j;
                }
            }
        }
        return level.size()-1; // error in finite bit arithmetic encountered
    }

    private DuoBaumLevel nextLevel() {
        ++arrayIndex;
        if (arrayIndex == levels.length) {
            ++windowIndex;
            arrayIndex = windowIndex;
        }
        return levels[arrayIndex];
    }

    private DuoBaumLevel currentLevel() {
        return levels[arrayIndex];
    }

    private DuoBaumLevel previousLevel(int sampleA, int sampleB) {
        if (arrayIndex == windowIndex) {
            --windowIndex;
            arrayIndex = windowIndex;
            levels[arrayIndex].setChildNodes(fwdNodes);
            int startLevel = levels[windowIndex].marker() + 1;
            int endLevel = startLevel + (levels.length - (windowIndex + 1) );
            for (int marker=startLevel; marker<endLevel; ++marker) {
                nextLevel().setForwardValues(fwdNodes, marker,
                        sampleA, sampleB);
            }
            return currentLevel();
        }
        else {
            return levels[--arrayIndex];
        }
    }

    private void forwardAlgorithm(int sampleA, int sampleB) {
        DuoBaumLevel.initializeNodes(fwdNodes);
        this.windowIndex = -1;
        this.arrayIndex = levels.length - 1;
        for (int marker=0; marker<nMarkers; ++marker) {
            nextLevel().setForwardValues(fwdNodes, marker, sampleA, sampleB);
        }
    }
}
