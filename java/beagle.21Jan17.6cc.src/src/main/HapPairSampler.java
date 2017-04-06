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
package main;

import blbutil.Utilities;
import dag.Dag;
import dag.LinkageEquilibriumDag;
import dag.MergeableDag;
import haplotype.BasicHapPairs;
import haplotype.HapPair;
import haplotype.HapPairs;
import haplotype.RevHapPairs;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import sample.ConsumeSingleSamples;
import sample.SingleBaum;
import vcf.GL;
import vcf.RevGL;

/**
 * <p>Class {@code HapPairSampler} samples haplotype pairs and
 * estimates posterior genotype probabilities.
 * </p>
 * <p>Instances of class {@code HapPairSampler} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class HapPairSampler {

    private final Par par;
    private final RunStats runStats;

    /**
     * Constructs a new {@code HapPairSampler} instance from the specified data.
     * @param par the analysis parameters
     * @param runStats the object to which run-time statistics will be written
     * @throws NullPointerException if {@code par == null || runStats == null}
     */
    public HapPairSampler(Par par, RunStats runStats) {
        if (par==null) {
            throw new IllegalArgumentException("par==null");
        }
        if (runStats==null) {
            throw new IllegalArgumentException("runStats==null");
        }
        this.par = par;
        this.runStats = runStats;
    }

    /**
     * Returns a list of sampled haplotype pairs.  Haplotype pairs are
     * sampled conditional on the observed genotype data and a haplotype
     * frequency model in which all markers are in linkage equilibrium.
     *
     * @param cd the input data for the current marker window
     * @return the a list of sampled haplotype pairs
     *
     * @throws NullPointerException if {@code cd == null}
     */
    public List<HapPair> initialHaps(CurrentData cd) {
        GL freqGL = cd.targetGL();
        GL emitGL = cd.targetGL();
        boolean useRevDag = false;
        float minAlleleFreq = 0.0001f;
        Dag dag = new LinkageEquilibriumDag(freqGL, minAlleleFreq);
        List<HapPair> sampledHaps = new ArrayList<>();
        sampledHaps = Collections.synchronizedList(sampledHaps);
        sample(dag, emitGL, useRevDag, par.nsamples(),
                sampledHaps, par.nthreads());
        return new ArrayList<>(sampledHaps);
    }

    /**
     * Returns a list of sampled haplotype pairs. Haplotype pairs are
     * sampled conditional on the observed genotype and a haplotype
     * frequency model constructed from the specified {@code hapPairs}.
     * The contract for this method is undefined if the specified
     * {@code hapPairs} and {@code gv} are inconsistent with the input data
     * contained in the {@code cd} parameter.
     *
     * @param cd the input data for the current marker window
     * @param hapPairs the haplotype pairs used to build the haplotype
     * frequency model
     * @param useRevDag {@code true} if the order of markers should
     * be reversed when building the haplotype frequency model, and
     * {@code false} otherwise
     * @param gv the current scaled genotype probabilities for the target
     * samples or {@code null} if genotype probabilities are not to be estimated
     * @return the sampled haplotype pairs
     *
     * @throws IllegalArgumentException if {@code haps.isEmpty() == true}
     * @throws NullPointerException if {@code cd == null || hapPairs == null}
     */
    public List<HapPair> sample(CurrentData cd, List<HapPair> hapPairs,
            boolean useRevDag, GenotypeValues gv) {
        if (hapPairs.isEmpty()) {
            throw new IllegalArgumentException("hapPairs.isEmpty()");
        }
        int nThreads =  par.nthreads();
        int nSampledHaps = par.nsamples()*cd.nTargetSamples();
        GL gl = gl(cd, useRevDag);
        Dag dag = getDagsAndUpdatePos(cd, hapPairs, useRevDag);
        List<HapPair> sampledHaps = synchronizedEmptyList(nSampledHaps);
        if (gv!=null) {
            if (useRevDag) {
                gv = new RevGenotypeValues(gv);
            }
            sample(dag, gl, useRevDag,par.nsamples(), sampledHaps, gv, nThreads);
        }
        else {
            sample(dag, gl, useRevDag, par.nsamples(), sampledHaps, nThreads);
        }
        return new ArrayList<>(sampledHaps);
    }

    private GL gl(CurrentData cd, boolean useRevDag) {
        GL gl = cd.targetGL();
        if (useRevDag) {
            gl = new RevGL(gl);
        }
        return gl;
    }

    private static List<HapPair> synchronizedEmptyList(int capacity) {
        List<HapPair> sampledHaps = new ArrayList<>(capacity);
        sampledHaps = Collections.synchronizedList(sampledHaps);
        return sampledHaps;
    }

    private Dag getDagsAndUpdatePos(CurrentData cd, List<HapPair> hapPairs,
            boolean useRevDag) {
        cd.addRestrictedRefHapPairs(hapPairs);
        HapPairs dagHaps = new BasicHapPairs(hapPairs);
        if (useRevDag) {
            dagHaps = new RevHapPairs(dagHaps);
        }
        float[] wts = cd.weights().get(dagHaps);
        Dag dag = makeDag(dagHaps, wts, par.modelscale());
        runStats.setDagStats(dag);
        return dag;
    }

    private Dag makeDag(HapPairs hapPairs, float[] weights, float scale) {
        long t0 = System.nanoTime();
        int nInitLevels = 500;
        Dag dag = MergeableDag.dag(hapPairs, weights, scale, nInitLevels);
        runStats.buildNanos(System.nanoTime() - t0);
        return dag;
    }

    @SuppressWarnings({"BroadCatchBlock", "TooBroadCatch"})
    private void sample(Dag dag, GL gl, boolean markersAreReversed,
            int nSamples, List<HapPair> sampledHaps, int nThreads) {
        long t0 = System.nanoTime();
        Random rand = new Random(par.seed());
        final BlockingQueue<Integer> qIn = new ArrayBlockingQueue<>(3*nThreads);
        ExecutorService es = Executors.newFixedThreadPool(nThreads);
        for (int j=0; j<nThreads; ++j) {
            SingleBaum sb = new SingleBaum(dag, gl, rand.nextLong(),
                    nSamples, par.lowmem());
            es.submit(new ConsumeSingleSamples(markersAreReversed, sb, qIn,
                    sampledHaps));
        }
        try {
            for (int j=0, n=gl.nSamples(); j<n; ++j) {
                qIn.put(j);
            }
            for (int j=0; j<nThreads; ++j) {
                qIn.put(ConsumeSingleSamples.POISON);
            }
            es.shutdown();
            es.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
        }
        catch (Throwable e) {
            Utilities.exit("ERROR", e);
        }
        runStats.sampleNanos(System.nanoTime() - t0);
    }

    @SuppressWarnings({"BroadCatchBlock", "TooBroadCatch"})
    private void sample(Dag dag, GL gl, boolean markersAreReversed, int nCopies,
            List<HapPair> sampledHaps, GenotypeValues gv,
            int nThreads) {
        long t0 = System.nanoTime();
        Random rand = new Random(par.seed());
        final BlockingQueue<Integer> qIn = new ArrayBlockingQueue<>(3*nThreads);
        ExecutorService es = Executors.newFixedThreadPool(nThreads);
        for (int j=0; j<nThreads; ++j) {
            SingleBaum sb = new SingleBaum(dag, gl, rand.nextLong(),
                    nCopies, par.lowmem());
            es.submit(new ConsumeSingleSamples(markersAreReversed, sb, qIn,
                    sampledHaps, gv));
        }
        try {
            for (int j=0, n=gl.nSamples(); j<n; ++j) {
                qIn.put(j);
            }
            for (int j=0; j<nThreads; ++j) {
                qIn.put(ConsumeSingleSamples.POISON);
            }
            es.shutdown();
            es.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
        }
        catch (Throwable e) {
            Utilities.exit("ERROR", e);
        }
        runStats.sampleNanos(System.nanoTime() - t0);
    }
}
