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

import blbutil.IntPair;
import blbutil.Utilities;
import dag.Dag;
import haplotype.HapPairs;
import haplotype.SampleHapPairs;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import main.GeneticMap;
import vcf.GL;

/**
 * <p>Class {@code HaploidIbd} implements the Refined IBD algorithm.
 * The Refined IBD algorithm detects candidate haplotype IBD segments with the
 * Germline Algorithm and then evaluates candidate IBD segments using a
 * likelihood ratio test.
 * </p>
 * <p>Instances of class {@code HaploidIbd} are immutable.
 *</p>
 * Reference: Gusev A, Lowe JK, Stoffel M, Daly MJ, Altshuler D, Breslow JL,
 *      Friedman JM, Pe'er I.  Whole population, genomewide mapping
 *      of hidden relatedness.  Genome Research 2009;19(2):318-26.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class HaploidIbd {

    private final GeneticMap genMap;
    private final int ibdTrim;
    private final float minIbdLod;
    private final float minFreqLod;
    private final float minCm;

    /**
     * Constructs a new {@code HaploidIbd} instance from the specified data.
     * @param genMap the genetic map
     * @param ibdTrim the number of markers to trim from an IBS segment
     * when computing the IBD versus non-IBD likelihood ratio
     * @param minIbdLod the minimum IBD LOD score of reported IBD segments
     * @param minCm the minimum cM length of reported IBD segments
     *
     * @throws IllegalArgumentException if {@code ibdTrim < 0 }
     * @throws IllegalArgumentException if
     * {@code ibdLod <= 0.0f || Float.isFinite(ibdLod) == false}
     * @throws IllegalArgumentException if
     * {@code minCm <= 0.0f || Float.isFinite(minCm) == false}
     * @throws NullPointerException if {@code genMap == null}
     */
    public HaploidIbd(GeneticMap genMap, int ibdTrim, float minIbdLod,
            float minCm) {
        if (genMap==null) {
            throw new IllegalArgumentException(GeneticMap.class.toString());
        }
        if (ibdTrim < 0) {
            throw new IllegalArgumentException(String.valueOf(ibdTrim));
        }
        if (minIbdLod <= 0.0 || Float.isFinite(minIbdLod) == false) {
            throw new IllegalArgumentException(String.valueOf(minIbdLod));
        }
        if (minCm <= 0.0 || Float.isFinite(minCm) == false) {
            throw new IllegalArgumentException(String.valueOf(minCm));
        }
        this.genMap = genMap;
        this.ibdTrim = ibdTrim;
        this.minIbdLod = minIbdLod;
        this.minFreqLod = minIbdLod;
        this.minCm = minCm;
    }

    /**
     * Runs the Refined IBD algorithm, and returns a map whose keys are
     * ordered pairs of haplotype indices and whose values are thread-safe
     * lists of IBD segments for each haplotype pair. The minimum haplotype
     * index is listed first in each ordered pair of haplotype indices.
     *
     * @param gl the HMM emission probabilities
     * @param dag the HMM transition probabilities
     * @param haps the sample haplotype pairs
     * @param nThreads the number of threads of execution that may be used
     * @return the detected IBD segments
     *
     * @throws IllegalArgumentException if {@code nThreads < 1}
     * @throws IllegalArgumentException if
     * {@code gl.samples().equals(haps.samples()) == false}
     * @throws IllegalArgumentException if
     * {@code gl.markers().equals(dag.markers()) == false
                || gl.markers().equals(haps.markers()) == false}
     * @throws NullPointerException if
     * {@code gl == null || dag == null || haps == null}
     */
    @SuppressWarnings({"BroadCatchBlock", "TooBroadCatch"})
    public Map<IntPair, List<IbdSegment>> run(GL gl, Dag dag,
            SampleHapPairs haps, final int nThreads) {
        checkParameters(gl, dag, haps);
        double[] pos = genMap.genPos(dag.markers());
        IbsHapSegments ibsSegments = new IbsHapSegments(haps, pos, minCm);
        ConcurrentMap<IntPair, List<IbdSegment>> ibdMap
                = new ConcurrentHashMap<>();

        final BlockingQueue<Integer> qIn = new ArrayBlockingQueue<>(5*nThreads);
        ExecutorService es = Executors.newFixedThreadPool(nThreads);
        for (int j=0; j<nThreads; ++j) {
            IbdBaum baum = new IbdBaum(dag, gl);
            es.submit(new ProduceIbd(haps, baum, ibsSegments, qIn, ibdMap,
                    ibdTrim, minIbdLod));
        }
        try {
            for (int hap=0, n=haps.nHaps(); hap<n; ++hap) {
                qIn.put(hap);
            }
            for (int j=0; j<nThreads; ++j) {
               qIn.put(ProduceIbd.POISON);
            }
            es.shutdown();
            es.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
        }
        catch (Throwable e) {
            Utilities.exit("ERROR", e);
        }
        return ibdMap;
    }

    private void checkParameters(GL gl, Dag dag, SampleHapPairs haps) {
        if (gl.samples().equals(haps.samples())==false) {
            throw new IllegalArgumentException("inconstent samples");
        }
        if (gl.markers().equals(dag.markers())==false
                || gl.markers().equals(haps.markers())==false) {
            throw new IllegalArgumentException("inconsistent markers");
        }
    }

    private static double freqLod(int hap, int start, int end, int ibdTrim,
            Dag dag, HapPairs haps) {
        int trimmedStart = start + ibdTrim;
        int trimmedEnd = end - ibdTrim;
        if (trimmedStart >= trimmedEnd) {
            return 0.0f;
        }
        else {
            return IbdBaum.freqLod(hap, trimmedStart, trimmedEnd, haps, dag);
        }
    }

    private static double ibdLod(IbdBaum ibdBaum, int hap1, int hap2, int start,
            int end, int ibdTrim) {
        int trimmedStart = start + ibdTrim;
        int trimmedEnd = end - ibdTrim;
        if (trimmedStart >= trimmedEnd) {
            return 0.0f;
        }
        else {
            int sample1 = hap1/2;
            int sample2 = hap2/2;
            return ibdBaum.ibdLod(sample1, sample2, trimmedStart, trimmedEnd);
        }
    }

    private class ProduceIbd implements Runnable {

        public static final int POISON = -37;

        private final SampleHapPairs haps;
        private final IbdBaum baum;
        private final IbsHapSegments ibsHapSegments;
        private final BlockingQueue<Integer> qIn;
        private final ConcurrentMap<IntPair, List<IbdSegment>> ibdMap;
        private final int ibdTrim;
        private final float minIbdLod;

        public ProduceIbd(SampleHapPairs haps, IbdBaum baum,
                IbsHapSegments ibsHapSegments, BlockingQueue<Integer> qIn,
                ConcurrentMap<IntPair, List<IbdSegment>> ibdMap, int ibdTrim,
                float minIbdLod) {
            if (ibdTrim < 0) {
                throw new IllegalArgumentException("trim < 0: " + ibdTrim);
            }
            if (minIbdLod <= 0.0 || Float.isNaN(minIbdLod)) {
                throw new IllegalArgumentException("ibdlod: " + minIbdLod);
            }
            this.haps = haps;
            this.baum = baum;
            this.ibsHapSegments = ibsHapSegments;
            this.qIn = qIn;
            this.ibdMap = ibdMap;
            this.ibdTrim = ibdTrim;
            this.minIbdLod = minIbdLod;
        }

        /*
         * Takes haplotype indices from a thread-safe work-queue and stores
         * detected IBD segments that between the haplotype and
         * haplotypes with larger index in {@code this.ibdMap}.  The method
         * exits when {@code ProduceSingleSamples.POISON} is taken from the
         * work queue.
         *
         * @throws IndexOutOfBounds exception if a negative integer
         * other than {@code ProduceSingleSamples.POISON} is taken from the
         * work queue
         */
        @Override
        @SuppressWarnings({"BroadCatchBlock", "TooBroadCatch"})
        public void run() {
            try {
                int hap = qIn.take();
                while (hap!=POISON) {
                    List<HapSegment> ibsSegs = ibsHapSegments.find(hap);
                    for (int j=0, n=ibsSegs.size(); j<n; ++j) {
                        HapSegment hs = ibsSegs.get(j);
                        if (hap < hs.hap()) {
                            int start = hs.start();
                            int end = hs.end();
                            double freqLod = HaploidIbd.freqLod(hap, start,
                                    (end+1), ibdTrim, baum.dag(), haps);
                            if (freqLod >= minFreqLod) {
                                float ibdLod;
                                if ( (hap/2) == (hs.hap()/2) ) {
                                    int sample = hap/2;
                                    ibdLod = (float) baum.hbdLod(sample, start, (end+1));
                                }
                                else {
                                    ibdLod = (float) HaploidIbd.ibdLod(baum, hap,
                                            hs.hap(), start, (end+1), ibdTrim);
                                }
                                if (ibdLod >= minIbdLod) {
                                    IntPair hapPair = new IntPair(hap, hs.hap());
                                    List<IbdSegment> list = ibdMap.get(hapPair);
                                    if (list==null) {
                                        list = Collections.synchronizedList(
                                                new ArrayList<IbdSegment>(2));
                                        ibdMap.putIfAbsent(hapPair, list);
                                        list = ibdMap.get(hapPair);
                                    }
                                    IbdSegment segment = new IbdSegment(hapPair,
                                            baum.gl().marker(start),
                                            baum.gl().marker(end),
                                            ibdLod, start, end );
                                    list.add(segment);
                                }
                            }
                        }
                    }
                    hap = qIn.take();
                }
            }
            catch (Throwable e) {
                Utilities.exit("ProduceSingleSamples: ERROR", e);
            }
        }
    }
}
