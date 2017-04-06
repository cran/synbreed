/*
 * Copyright (C) 2014 Brian L. Browning
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package main;

import blbutil.Utilities;
import haplotype.HapPair;
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
import sample.RecombSingleBaum;
import sample.SamplerData;
import sample.SingleBaumInterface;

/**
 * <p>Class {@code RecombHapPairSamples} samples haplotype pairs and
 * estimates posterior genotype probabilities using a haplotype frequency
 * model that permits transitions between any two states at adjacent markers.
 * </p>
 * <p>Instances of class {@code RecombHapPairSampler} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class RecombHapPairSampler {

    private static final int nCopies = 4;

    private final Par par;
    private final RunStats runStats;
    private double edgePairsPerMarker;

    /**
     * Constructs a new {@code RecombHapPairSampler} instance from the
     * specified data.
     * @param par the analysis parameters
     * @param runStats the object to which run-time statistics will be written
     * @throws NullPointerException if
     * {@code par == null || runStats == null}
     */
    public RecombHapPairSampler(Par par, RunStats runStats) {
        if (par==null) {
            throw new NullPointerException("par");
        }
        if (runStats==null) {
            throw new NullPointerException("runStats");
        }
        this.par = par;
        this.runStats = runStats;
        this.edgePairsPerMarker = 0;
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
     * @param hapPairs the target haplotype pairs used to build the haplotype
     * frequency model
     * @param useRevDag {@code true} if the order of markers should
     * be reversed when building the haplotype frequency model, and
     * {@code false} otherwise
     * @param gv the current scaled genotype probabilities for the target
     * samples or {@code null} if genotype probabilities are not to be estimated
     * @return the sampled haplotype pairs
     *
     * @throws IllegalArgumentException if {@code haps.isEmpty() == true}
     * @throws NullPointerException if  {@code cd == null || hapPairs == null}
     */
    public List<HapPair> sample(CurrentData cd, List<HapPair> hapPairs,
            boolean useRevDag, GenotypeValues gv) {
        SamplerData samplerData = new SamplerData(par, cd, hapPairs, useRevDag,
                runStats);
        int nSampledHaps = nCopies*cd.nTargetSamples();
        List<HapPair> sampledHaps = synchronizedEmptyList(nSampledHaps);
        if (gv!=null) {
            if (useRevDag) {
                gv = new RevGenotypeValues(gv);
            }
            sample(samplerData, sampledHaps, gv);
        }
        else {
            sample(samplerData, sampledHaps);
        }
        return new ArrayList<>(sampledHaps);
    }

    private static List<HapPair> synchronizedEmptyList(int capacity) {
        List<HapPair> sampledHaps = new ArrayList<>(capacity);
        sampledHaps = Collections.synchronizedList(sampledHaps);
        return sampledHaps;
    }

    @SuppressWarnings({"BroadCatchBlock", "TooBroadCatch"})
    private void sample(SamplerData samplerData,  List<HapPair> sampledHaps,
            GenotypeValues gv) {
        long t0 = System.nanoTime();
        int nThreads = samplerData.par().nthreads();
        boolean markersAreReversed = samplerData.markersAreReversed();
        Random rand = new Random(par.seed());
        final BlockingQueue<Integer> qIn = new ArrayBlockingQueue<>(3*nThreads);
        ExecutorService es = Executors.newFixedThreadPool(nThreads);
        for (int j=0; j<nThreads; ++j) {
            SingleBaumInterface sb = new RecombSingleBaum(samplerData,
                    rand.nextLong(), nCopies, par.lowmem());
            es.submit(new ConsumeSingleSamples(markersAreReversed, sb, qIn,
                        sampledHaps, gv));
        }
        try {
            for (int j=0, n=samplerData.nSamples(); j<n; ++j) {
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
    private void sample(SamplerData samplerData, List<HapPair> sampledHaps) {
        long t0 = System.nanoTime();
        int nThreads = samplerData.par().nthreads();
        boolean markersAreReversed = samplerData.markersAreReversed();
        Random rand = new Random(par.seed());
        final BlockingQueue<Integer> qIn = new ArrayBlockingQueue<>(3*nThreads);
        ExecutorService es = Executors.newFixedThreadPool(nThreads);
        for (int j=0; j<nThreads; ++j) {
            SingleBaumInterface sb = new RecombSingleBaum(samplerData,
                    rand.nextLong(), nCopies, par.lowmem());
            es.submit(new ConsumeSingleSamples(markersAreReversed, sb, qIn,
                        sampledHaps));
        }
        try {
            for (int j=0, n=samplerData.nSamples(); j<n; ++j) {
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
