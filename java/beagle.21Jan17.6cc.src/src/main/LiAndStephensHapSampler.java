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
import haplotype.SampleHapPairs;
import java.util.Queue;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import sample.LSHapBaum;
import sample.ImputationData;

/**
 * <p>Class {@code LiAndStephensHapSampler} estimates posterior allele probabilities.
 * </p>
 * <p>Instances of class {@code LiAndStephensHapSampler} are thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class LiAndStephensHapSampler {

    private final Par par;
    private final GeneticMap genMap;

    /**
     * Constructs a {@code LiAndStephensHapSampler} instance from the specified
     * data.
     * @param par the analysis parameters
     * @param genMap the genetic map or {@code null} if no genetic map is
     * specified.
     * @throws NullPointerException if {@code par == null}
     */
    public LiAndStephensHapSampler(Par par, GeneticMap genMap) {
        if (par==null) {
            throw new IllegalArgumentException("par==null");
        }
        this.par = par;
        this.genMap = genMap;
    }

    /**
     * Returns estimated allele probabilities for each target sample.
     * The contract for this method is undefined if the data in the
     * specified {@code shp} and {@code cd} parameters are inconsistent.
     * @param cd the input data for the current marker window
     * @param shp estimated target haplotypes at the genotyped markers
     * @return estimated allele probabilities for each target sample
     * @throws NullPointerException if {@code cd == null || shp == null}
     */
    public BasicAlleleProbs sample(CurrentData cd, SampleHapPairs shp)  {
        Queue<HapAlleleProbs> qOut = new ConcurrentLinkedQueue<>();
        multiThreadedHapSample(cd, shp, qOut, par.lowmem(), par.nthreads());
        HapAlleleProbs[] hapAlleleProbs = qOut.toArray(new HapAlleleProbs[0]);
        return new BasicAlleleProbs(hapAlleleProbs);
    }

    @SuppressWarnings({"BroadCatchBlock", "TooBroadCatch"})
    private void multiThreadedHapSample(CurrentData cd,
            SampleHapPairs targetHapPairs, Queue<HapAlleleProbs> qOut,
            boolean lowMem, int nThreads) {
        int qInSize = targetHapPairs.nSamples() + nThreads;
        final BlockingQueue<Integer> qIn = new ArrayBlockingQueue<>(qInSize);
        ExecutorService es = Executors.newFixedThreadPool(nThreads);
        ImputationData impData = new ImputationData(par, cd, targetHapPairs,
                genMap);
        try {
            for (int j=0, n=targetHapPairs.nSamples(); j<n; ++j) {
                qIn.put(j);
            }
            for (int j=0; j<nThreads; ++j) {
                qIn.put(LSHapSampler.POISON);
            }
            for (int j=0; j<nThreads; ++j) {
                LSHapBaum hb = new LSHapBaum(impData, lowMem);
                es.submit(new LSHapSampler(hb, qIn, qOut));
            }
            es.shutdown();
            es.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
        }
        catch (Throwable e) {
            Utilities.exit("RecombHapSampler: ERROR", e);
        }
    }

    private class LSHapSampler implements Runnable {

        public static final int POISON = -1;

        private final LSHapBaum baum;
        private final BlockingQueue<Integer> qIn;
        private final Queue<HapAlleleProbs> qOut;

        /*
         * Constructs a {@code ProduceHapSample} instance.
         *
         * @param markersAreReversed {@code true} if sampled haplotypes
         * will have their marker order reversed and {@code false} otherwise.
         * @param baum a thread-confined instance of class
         * {@code sample.HaplotypeBaum}.
         * @param qIn a thread-safe input work queue.
         * @param hapList a thread-safe list for storing sampled haplotype pairs.
         * @param gv a thread-safe object which stores genotype posterior probabilities.
         *
         * @throws NullPointerException if any parameter is {@code null}.
         */
        public LSHapSampler(LSHapBaum baum, BlockingQueue<Integer> qIn,
                Queue<HapAlleleProbs> qOut) {
            if (baum==null) {
                throw new NullPointerException("baum=null");
            }
            if (qIn==null) {
                throw new IllegalArgumentException("qIn==null");
            }
            this.baum = baum;
            this.qIn = qIn;
            this.qOut = qOut;
        }

        @Override
        @SuppressWarnings({"BroadCatchBlock", "TooBroadCatch"})
        public void run() {
            try {
                int sample = qIn.take();
                while (sample != POISON) {
                    int hap1 = 2*sample;
                    int hap2 = 2*sample+1;
                    qOut.add(baum.randomHapSample(hap1));
                    qOut.add(baum.randomHapSample(hap2));
                    sample = qIn.take();
                }
            }
            catch (Throwable e) {
                Utilities.exit("RecombHapSampler.LSHapSampler: ERROR", e);
            }
        }
    }
}
