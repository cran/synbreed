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
import haplotype.HapPair;
import haplotype.BitHapPair;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import vcf.AL;

/**
 * Class {@code HapBaum}  implements the Baum forward and backward
 * algorithms for a hidden Markov model (HMM) for an individual's genotype data.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class HapBaum {

    private final Dag dag;
    private final AL al;
    private final int nMarkers;
    private final int nCopies;
    private final long seed;
    private final Random random;

    private final int[] node;
    private final double[] nodeValue;

    private final byte[][] alleles1;
    private final byte[][] alleles2;

    private final double[] alProbs1;
    private final double[] alProbs2;

    private final HapBaumLevel[] levels;
    private final HapNodes fwdNodes;
    private final HapNodes bwdNodes;

    private int windowIndex = -9999;
    private int arrayIndex = -9999;

    /**
     * Creates a {@code HaplotypeBaum}  instance.
     *
     * @param dag the directed acyclic graph that determines the
     * transition probabilities.
     * @param al the emission probabilities.
     * @param seed the initial random seed.
     * @param nCopies the number of haplotype pairs that will be sampled for
     * each individual.
     *
     * @throws IllegalArgumentException
     * if {@code nCopies<1 || dag.markers().equals(al.markers())==false}
     * @throws NullPointerException if {@code dag==null || al==null}
     */
    public HapBaum(Dag dag, AL al, long seed, int nCopies) {
        if (dag.markers().equals(al.markers())==false) {
            throw new IllegalArgumentException("inconsistent markers");
        }
        if (nCopies < 1) {
            throw new IllegalArgumentException("nCopies<1: " + nCopies);
        }
        this.dag = dag;
        this.al = al;
        this.nMarkers = dag.nMarkers();
        this.nCopies = nCopies;
        this.seed = seed;
        this.random = new Random(seed);

        this.node = new int[nCopies];
        this.nodeValue = new double[nCopies];
        this.alleles1 = new byte[nCopies][al.nMarkers()];
        this.alleles2 = new byte[nCopies][al.nMarkers()];
        this.alProbs1 = new double[al.markers().sumAlleles()];
        this.alProbs2 = new double[al.markers().sumAlleles()];

        int size = (int) Math.ceil(Math.sqrt(1 + 8*dag.nMarkers())/2.0) + 1;
        this.levels = new HapBaumLevel[size];
        for (int j=0; j<levels.length; ++j) {
            levels[j] = new HapBaumLevel(dag, al);
        }
        this.fwdNodes = new HapNodes();
        this.bwdNodes = new HapNodes();
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
    public AL al() {
        return al;
    }

    /**
     * Returns the number of haplotype pairs that are sampled for each
     * individual.
     * @return the number of haplotype pairs that are sampled for each
     * individual.
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
     * <p>Returns a list of {@code this.nCopies()} sampled haplotype pairs for the
     * specified individual. Haplotype pairs are sampled conditional on the
     * HMM with transition probabilities determined by {@code this.dag()} and
     * emission probabilities determined by {@code this.al()}.
     * </p>
     * The contract for this method is unspecified if no haplotype pair
     * is consistent with the HMM.
     *
     * @param sample a sample index.
     * @return a list of {@code this.nCopies()} sampled haplotype pairs for the
     * specified individual.
     *
     * @throws IndexOutOfBoundsException if
     * {@code sample<0 || sample>=this.al().nSamples()}
     */
    public List<HapPair> randomSample(int sample) {
        int hap = 2*sample;
        randomSample(hap, alleles1);
        randomSample(++hap, alleles2);
        return hapList(sample);
    }

    private void randomSample(int hap, byte[][] alleles) {
        forwardAlgorithm(hap);
        initSampleAlleles(currentLevel(), hap, alleles);
        for (int j=nMarkers-2; j>=0; --j) {
            HapBaumLevel level = previousLevel(hap);
            sampleAlleles(level, hap, alleles);
        }
    }

    /**
     * <p>Returns a list of {@code this.nCopies()} sampled haplotype pairs for the
     * specified individual. Haplotype pairs are sampled conditional on the
     * HMM with transition probabilities determined by {@code this.dag()} and
     * emission probabilities determined by {@code this.al()}.
     * Posterior genotype probabilities  are written to the specified array.
     * The posterior probability of the {@code j}-th genotype for
     * the {@code k}-th marker is stored at index
     * {@code gl.markers().sumGenotypes(k) + j} in the {@code gprobs}
     * array.
     * </p>
     * The contract for this method is unspecified if no haplotype pair
     * is consistent with the HMM.
     *
     * @param sample the sample index.
     * @param gtProbs a array to which posterior genotype probabilities
     * for the sample will be written.
     * @return a list of {@code this.nCopies()} sampled haplotype pairs for the
     * specified individual.
     *
     * @throws IndexOutOfBoundsException if
     * {@code sample<0 || sample>=this.al().nSamples()}
     * @throws IllegalArgumentException if
     * {@code gprobs.length!=this.al().markers().sumGenotypes()}
     * @throws NullPointerException if {@code gprobs==null}
     */
    public List<HapPair> randomSample(int sample, double[] gtProbs) {
        checkGprobs(gtProbs);
        int hap = 2*sample;
        randomSample(hap, alleles1, alProbs1);
        randomSample(++hap, alleles2, alProbs2);
        setGprobs(gtProbs);
        return hapList(sample);
    }

    private void checkGprobs(double[] gtProbs) {
        if (gtProbs.length != al.markers().sumGenotypes()) {
            String s = "gtProbs.length!=al.markers().sumGenotypes()";
            throw new IllegalArgumentException(s);
        }
    }

    private void randomSample(int hap, byte[][] alleles, double[] alProbs) {
        forwardAlgorithm(hap);
        initSampleAlleles(currentLevel(), hap, alleles);
        currentLevel().setInitialBackwardValues(bwdNodes);
        setAlProbs(currentLevel(), alProbs);
        for (int j=nMarkers-2; j>=0; --j) {
            HapBaumLevel level = previousLevel(hap);
            sampleAlleles(level, hap, alleles);
            level.setBackwardValues(bwdNodes);
            setAlProbs(level, alProbs);
        }
    }

    private void setGprobs(double[] gtProbs) {
        int index = 0;
        int alEnd = 0;
        for (int m=0; m<nMarkers; ++m) {
            int alStart = alEnd;
            alEnd = al.markers().sumAlleles(m+1);
            for (int a2=alStart; a2<alEnd; ++a2) {
                for (int a1=alStart; a1<a2; ++a1) {
                    gtProbs[index++] = (alProbs1[a1]*alProbs2[a2]
                            + alProbs1[a2]*alProbs2[a1]);
                }
                gtProbs[index++] = alProbs1[a2]*alProbs2[a2];
            }
        }
        assert index==gtProbs.length;
    }

    private void setAlProbs(HapBaumLevel level, double[] alProbs) {
        if (alProbs != null) {
            int m = level.marker();
            int nAlleles = al.marker(m).nAlleles();
            int base = al.markers().sumAlleles(m);
            for (int j=0; j<nAlleles; ++j) {
                alProbs[base + j] = level.alProbs(j);
            }
        }
    }

    private List<HapPair> hapList(int sample) {
        List<HapPair> hapList = new ArrayList<>(2*nCopies);
        for (int copy=0; copy<nCopies; ++copy) {
            HapPair haps = new BitHapPair(al.markers(),
                    al.samples().idIndex(sample),
                    alleles1[copy], alleles2[copy]);
            hapList.add(haps);
        }
        return hapList;
    }

    private void initSampleAlleles(HapBaumLevel level, int hap, byte[][] alleles) {
        int m = level.marker();
        for (int copy=0; copy<nCopies; ++copy) {
            int state = initialRandomState(level);
            node[copy] = level.parentNode(state);
            nodeValue[copy] =  parentSum(level, hap, state);
            alleles[copy][m] = level.symbol(state);
        }
    }

    private int initialRandomState(HapBaumLevel level) {
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

    private double parentSum(HapBaumLevel level, int hap, int state) {
        int marker = level.marker();
        double fwdValue = level.forwardValuesSum()*level.forwardValue(state);
        int edge = level.edge(state);
        double tp = dag.condEdgeProb(marker, edge);
        byte symbol = dag.symbol(marker, edge);
        double ep = al.al(marker, hap, symbol);
        return fwdValue / ( ep*tp );
    }

    private void sampleAlleles(HapBaumLevel level, int hap, byte[][] alleles) {
        int m = level.marker();
        for (int copy=0; copy<nCopies; ++copy) {
            int state = randomPreviousState(level, node[copy], nodeValue[copy]);
            node[copy] = level.parentNode(state);
            nodeValue[copy] =  parentSum(level, hap, state);
            alleles[copy][m] = level.symbol(state);
        }
    }

    private int randomPreviousState(HapBaumLevel level, int node,
            double nodeValue) {
        double d = random.nextDouble() * nodeValue;
        double sum = 0.0;
        for (int j=0, n=level.size(); j<n; ++j) {
            if ( (node==level.childNode(j)) ) {
                sum += level.forwardValue(j);
                if (d <= sum) {
                    return j;
                }
            }
        }
        return level.size()-1; // error in finite bit arithmetic encountered
    }

    private HapBaumLevel nextLevel() {
        ++arrayIndex;
        if (arrayIndex == levels.length) {
            ++windowIndex;
            arrayIndex = windowIndex;
        }
        return levels[arrayIndex];
    }

    private HapBaumLevel currentLevel() {
        return levels[arrayIndex];
    }

    private HapBaumLevel previousLevel(int sample) {
        if (arrayIndex == windowIndex) {
            --windowIndex;
            arrayIndex = windowIndex;
            levels[arrayIndex].setChildNodes(fwdNodes);
            int startLevel = levels[windowIndex].marker() + 1;
            int endLevel = startLevel + (levels.length - (windowIndex + 1) );
            for (int marker=startLevel; marker<endLevel; ++marker) {
                nextLevel().setForwardValues(fwdNodes, marker, sample);
            }
            return currentLevel();
        }
        else {
            return levels[--arrayIndex];
        }
    }

    private void forwardAlgorithm(int hap) {
        HapBaumLevel.initializeNodes(fwdNodes);
        this.windowIndex = -1;
        this.arrayIndex = levels.length - 1;
        for (int marker=0; marker<nMarkers; ++marker) {
            nextLevel().setForwardValues(fwdNodes, marker, hap);
        }
    }
}
