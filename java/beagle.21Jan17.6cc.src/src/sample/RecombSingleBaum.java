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
import haplotype.HapPair;
import haplotype.BitHapPair;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import vcf.GL;

/**
 * <p>Class {@code RestrictedSingleBaum} implements the Baum forward and
 * backward algorithms for a hidden Markov model (HMM) of an individual's
 * genotype data. The HMM transition probabilities model recent
 * genetic recombination by allowing jumps between states that are not
 * connected by a node.
 * </p>
 * <p>Instances of class {@code RestrictedSingleBaum} are not thread-safe.
 * </p>
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class RecombSingleBaum implements SingleBaumInterface {

    private final SamplerData samplerData;
    private final Dag dag;
    private final RestrictedDag rdag;
    private final GL gl;
    private final int nMarkers;
    private final int nSamplesPerIndividual;
    private final long seed;
    private final Random random;

    private final int[] node1;
    private final int[] node2;
    private final float[] baseTrProb;
    private final float[] maxSum;

    private final int[][] alleles1;
    private final int[][] alleles2;

    private final RecombSingleBaumLevel[] levels;
    private final RecombSingleNodes fwdNodes;
    private final RecombSingleNodes bwdNodes;

    private int windowIndex = -9999;
    private int arrayIndex = -9999;

    /**
     * Creates a new {@code RestrictedSingleBaum} instance from the specified
     * data.
     *
     * @param samplerData the analysis data
     * @param seed the random seed
     * @param nSamplesPerIndividual the number of haplotype pairs that
     * will be sampled for each individual
     * @param lowMem {@code true} if a low memory algorithm should be used, and
     * {@code false} otherwise
     *
     * @throws IllegalArgumentException if {@code nSamplesPerIndividual < 1}
     * @throws NullPointerException if {@code samplerData == null}
     */
    public RecombSingleBaum(SamplerData samplerData, long seed,
            int nSamplesPerIndividual, boolean lowMem) {
        if (nSamplesPerIndividual < 1) {
            throw new IllegalArgumentException(
                    String.valueOf(nSamplesPerIndividual));
        }
        this.samplerData = samplerData;
        this.dag = samplerData.rdag().dag();
        this.rdag = samplerData.rdag();
        this.gl = samplerData.gl();
        this.nMarkers = samplerData.nMarkers();
        this.nSamplesPerIndividual = nSamplesPerIndividual;
        this.seed = seed;
        this.random = new Random(seed);

        this.node1 = new int[nSamplesPerIndividual];
        this.node2 = new int[nSamplesPerIndividual];
        this.baseTrProb = new float[nSamplesPerIndividual];
        this.maxSum = new float[nSamplesPerIndividual];

        this.alleles1 = new int[nSamplesPerIndividual][nMarkers];
        this.alleles2 = new int[nSamplesPerIndividual][nMarkers];

        int size = dag.nLevels();
        if (lowMem) {
            size = (int) Math.ceil(Math.sqrt(1 + 8*dag.nLevels())/2.0) + 1;
        }
        this.levels = new RecombSingleBaumLevel[size];
        for (int j=0; j<levels.length; ++j) {
            levels[j] = new RecombSingleBaumLevel(samplerData);
        }
        this.fwdNodes = new RecombSingleNodes(dag.maxNodes());
        this.bwdNodes = new RecombSingleNodes(dag.maxNodes());
    }

    @Override
    public Dag dag() {
        return rdag.dag();
    }

    @Override
    public GL gl() {
        return gl;
    }

    @Override
    public int nSamplesPerIndividual() {
        return nSamplesPerIndividual;
    }

    @Override
    public long seed() {
        return seed;
    }

    @Override
    public List<HapPair> randomSample(int sample) {
        DiploidStates permittedStates = rdag.singleStates(sample);
        forwardAlgorithm(sample, permittedStates);
        initSampleAlleles(currentLevel(), sample);
        for (int j=nMarkers-2; j>=0; --j) {
            RecombSingleBaumLevel level
                    = previousLevel(sample, permittedStates);
            sampleAlleles(level, sample);
        }
        pruneLevels();
        return hapList(sample);
    }

    @Override
    public List<HapPair> randomSample(int sample, double[] gprobs) {
        checkGprobs(gprobs);
        DiploidStates permittedStates = rdag.singleStates(sample);
        forwardAlgorithm(sample, permittedStates);
        initSampleAlleles(currentLevel(), sample);
        currentLevel().setInitialBackwardValues(bwdNodes);
        setGprobs(currentLevel(), gprobs);
        for (int j=nMarkers-2; j>=0; --j) {
            RecombSingleBaumLevel level = previousLevel(sample, permittedStates);
            sampleAlleles(level, sample);
            level.setBackwardValues(bwdNodes);
            setGprobs(level, gprobs);
        }
        pruneLevels();
        return hapList(sample);
    }

    private void pruneLevels() {
        int meanSize = estMeanSize();
        int capacityThreshold = 3*meanSize;
        int newCapacity = 3*meanSize/2 + 1;
        for (int j=0; j<levels.length; ++j) {
            if (levels[j].capacity() > capacityThreshold) {
                levels[j].reset(newCapacity);
            }
        }
    }

    private int estMeanSize() {
        int nLevelsToSample = 20;
        long sizeSum = 0;
        for (int j=0; j<nLevelsToSample; ++j) {
            sizeSum += levels[random.nextInt(levels.length)].size();
        }
        return (int) (sizeSum / nLevelsToSample);
    }

    private void checkGprobs(double[] gprobs) {
        int n = gl.markers().sumGenotypes();
        if (gprobs.length != n) {
            throw new IllegalArgumentException(String.valueOf(n));
        }
    }

    private void setGprobs(RecombSingleBaumLevel level, double[] gprobs) {
        if (gprobs != null) {
            int m = level.marker();
            int nGenotypes = gl.marker(m).nGenotypes();
            int base = gl.markers().sumGenotypes(m);
            for (int j=0; j<nGenotypes; ++j) {
                gprobs[base + j] = level.gprobs(j);
            }
        }
    }

    private List<HapPair> hapList(int sample) {
        List<HapPair> hapList = new ArrayList<>(2*nSamplesPerIndividual);
        for (int j=0; j<nSamplesPerIndividual; ++j) {
            HapPair haps = new BitHapPair(gl.markers(), gl.samples(), sample,
                    alleles1[j], alleles2[j]);
            hapList.add(haps);
        }
        return hapList;
    }

    private void initSampleAlleles(RecombSingleBaumLevel level, int sample) {
        for (int j=0; j<nSamplesPerIndividual; ++j) {
            saveCurrentData(level, sample, j, initialRandomState(level));
        }
    }

    private int initialRandomState(RecombSingleBaumLevel level) {
        float d = random.nextFloat();
        float sum = 0.0f;
        for (int j=0, n=level.size(); j<n; ++j) {
            sum += level.forwardValue(j);
            if (d <= sum) {
                return j;
            }
        }
        return level.size()-1; // error in finite bit arithmetic encountered
    }

    private void saveCurrentData(RecombSingleBaumLevel level, int sample,
            int copy, int stateIndex) {
        int m = level.marker();
        int e1 = level.edge1(stateIndex);
        int e2 = level.edge2(stateIndex);
        int s1 = level.symbol1(stateIndex);
        int s2 = level.symbol2(stateIndex);
        node1[copy] = level.parentNode1(stateIndex);
        node2[copy] = level.parentNode2(stateIndex);
        float p1 = dag.edgeProb(m, e1);
        float p2 = dag.edgeProb(m, e2);
        baseTrProb[copy] = p1*p2;

        maxSum[copy] = level.forwardValue(stateIndex) * level.forwardValuesSum()
                / gl.gl(m, sample, s1, s2);
        alleles1[copy][m] = s1;
        alleles2[copy][m] = s2;
    }

    private void sampleAlleles(RecombSingleBaumLevel level, int sample) {
        for (int j=0; j<nSamplesPerIndividual; ++j) {
            saveCurrentData(level, sample, j, randomPreviousState(level, j));
        }
    }

    private int randomPreviousState(RecombSingleBaumLevel level, int copy) {
        int m = level.marker();
        float np1 = dag.parentProb(m+1, node1[copy]);
        float np2 = dag.parentProb(m+1, node2[copy]);
        float pRecomb = samplerData.pRecomb(m+1);
        float d = random.nextFloat() * maxSum[copy];
        float sum = 0.0f;
        for (int j=0, n=level.size(); j<n; ++j) {
            float tp = 0.0f;
            boolean noJump1 = level.childNode1(j)==node1[copy];
            boolean noJump2 = level.childNode2(j)==node2[copy];
            if (noJump1 && noJump2) {
                tp += (1-pRecomb)*(1-pRecomb)*baseTrProb[copy]/ (np1*np2);
            }
            if (noJump1) {
                tp += (1-pRecomb)*pRecomb*baseTrProb[copy] / np1;
            }
            if (noJump2) {
                tp += pRecomb*(1-pRecomb)*baseTrProb[copy] / np2;
            }
            tp += pRecomb*pRecomb*baseTrProb[copy];

            sum += (level.forwardValue(j)*tp);
            if (d <= sum) {
                return j;
            }
        }
        return level.size()-1; // if reached due to rounding
    }

    private RecombSingleBaumLevel nextLevel() {
        ++arrayIndex;
        if (arrayIndex == levels.length) {
            ++windowIndex;
            arrayIndex = windowIndex;
        }
        return levels[arrayIndex];
    }

    private RecombSingleBaumLevel currentLevel() {
        return levels[arrayIndex];
    }

    private RecombSingleBaumLevel previousLevel(int sample,
            DiploidStates permittedStates) {
        if (arrayIndex == windowIndex) {
            --windowIndex;
            arrayIndex = windowIndex;
            levels[arrayIndex].setChildNodes(fwdNodes);
            int startLevel = levels[windowIndex].marker() + 1;
            int endLevel = startLevel + (levels.length - (windowIndex + 1) );
            for (int marker=startLevel; marker<endLevel; ++marker) {
                nextLevel().setForwardValues(fwdNodes, permittedStates, marker, sample);
            }
            return currentLevel();
        }
        else {
            return levels[--arrayIndex];
        }
    }

    private void forwardAlgorithm(int sample, DiploidStates permittedStates) {
        fwdNodes.clear();
        fwdNodes.sumUpdate(0, 0, 1.0f);
        this.windowIndex = -1;
        this.arrayIndex = levels.length - 1;
        for (int marker=0; marker<nMarkers; ++marker) {
            nextLevel().setForwardValues(fwdNodes, permittedStates, marker,
                    sample);
        }
    }
}
