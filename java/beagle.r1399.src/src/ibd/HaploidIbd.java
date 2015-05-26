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
import vcf.GL;

/**
 * <p>Class {@code HaploidIbd} implements the Refined IBD algorithm.
 * The Refined IBD algorithm detects candidate haplotype IBD segments with the
 * Germline Algorithm and then evaluates candidate IBD segments using a
 * likelihood ratio test.
 *</p>
 * Reference: Gusev A, Lowe JK, Stoffel M, Daly MJ, Altshuler D, Breslow JL,
 *      Friedman JM, Pe'er I.  Whole population, genomewide mapping
 *      of hidden relatedness.  Genome Research 2009;19(2):318-26.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class HaploidIbd {

    private final int ibdTrim;
    private final float minIbdLod;
    private final float minIbsLength;   // positions from Dag.posArray()
    private final float minFreqLod;     // for shared haplotype

    /**
     * Constructs a new {@code HaploidIbd} instance.
     * @param ibdTrim the number of markers to trim from an IBS segment
     * when computing a IBD vs. non-IBD likelihood ratio.
     * @param minIbdLod the minimum IBD LOD score of reported IBD segments
     *
     * @throws IllegalArgumentException if {@code ibdTrim<0}
     * @throws IllegalArgumentException if
     * {@code ibdLod<=0.0f || Float.isNaN(ibdLod)==true}
     */
    public HaploidIbd(int ibdTrim, float minIbdLod) {
        if (ibdTrim < 0) {
            throw new IllegalArgumentException("trim: " + ibdTrim);
        }
        if (minIbdLod <= 0.0 || Float.isNaN(minIbdLod)) {
            throw new IllegalArgumentException("minIbdlod: " + minIbdLod);
        }
        this.ibdTrim = ibdTrim;
        this.minIbdLod = minIbdLod;
        this.minIbsLength = 0.8f*minIbdLod;
        this.minFreqLod = minIbdLod;
    }

    /**
     * Runs the Refined IBD algorithm, and returns the detected IBD segments.
     *
     * @param gl the HMM emission probabilities.
     * @param dag the HMM transition probabilities.
     * @param haps the sample haplotype pairs.
     * @param nThreads the number of threads of execution that may be used.
     * @return a map whose keys are pairs of haplotype indices and whose
     * values are thread-safe lists of IBD segments for the haplotype pairs.
     *
     * @throws IllegalArgumentException if {@code nThreads<1}
     * @throws IllegalArgumentException if
     * {@code gl.samples().equals(haps.samples())==false}
     * @throws IllegalArgumentException if
     * {@code gl.markers().equals(dag.markers())==false
                || gl.markers().equals(haps.markers())==false}
     * @throws NullPointerException if
     * {@code gl==null || dag==null || haps==null}
     */
    @SuppressWarnings({"BroadCatchBlock", "TooBroadCatch"})
    public Map<IntPair, List<IbdSegment>> run(GL gl, Dag dag,
            SampleHapPairs haps, final int nThreads) {
        checkParameters(gl, dag, haps);
        double[] pos = dag.posArray();
        IbsHapSegments ibsSegments = new IbsHapSegments(haps, pos, minIbsLength);
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
            Utilities.exit("\"HapSampler: ERROR\"", e);
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
            return IbdBaum.freqLod(hap, trimmedStart, trimmedEnd, dag, haps);
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

        /**
         * Constructs a {@code ProduceIbd} instance.
         *
         * @param haps the haplotypes.
         * @param baum a thread-confined instance of class {@code ibd.IbdBaum}.
         * @param qIn a thread-safe input work queue.
         * @param ibdMap a thread-safe map whose keys are pairs of haplotype
         * indices, and whose values are thread-safe lists of IBD segments
         * for the haplotype pairs.
         * @param ibdTrim the number of markers to trim from an IBS segment
         * when computing an IBD vs. non-IBD likelihood ratio.
         * @param minIbdLod the minimum IBD LOD score of reported IBD segments.
         *
         * @throws IllegalArgumentException if {@code ibdTrim<0}
         * @throws IllegalArgumentException if
         * {@code minIbdLod<=0.0 || Float.isNan(isIbdLod)}
         * @throws NullPointerException if any parameter is {@code null}
         */
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

        /**
         * Takes haplotype indices from a thread-safe work-queue and stores
         * detected IBD segments that between the haplotype and
         * haplotypes with large index in {@code this.ibdMap}.  The method
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
