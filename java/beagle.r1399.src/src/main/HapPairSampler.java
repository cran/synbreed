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
import haplotype.RevHapPair;
import haplotype.Weights;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import sample.DuoBaum;
import sample.HapBaum;
import sample.ProduceHapSamples;
import sample.ProduceSingleSamples;
import sample.SingleBaum;
import sample.TrioBaum;
import vcf.AL;
import vcf.GL;
import vcf.Markers;
import vcf.RevAL;
import vcf.RevGL;

/**
 * Class {@code HapPairSampler} samples haplotype pairs and estimates posterior
 * genotype probabilities.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class HapPairSampler {

    private final Parameters par;
    private final RunStats runStats;

    /**
     * Constructs a new {@code HapPairSampler} instance.
     * @param par the analysis parameters.
     * @param runStats the object to which run-time statistics will be written.
     * @throws NullPointerException if {@code par==null || runStats==null}
     */
    public HapPairSampler(Parameters par, RunStats runStats) {
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
     * sampled under a model with markers in linkage equilibrium.
     *
     * @param fam the parent-offspring relationships.
     * @param freqGL the genotype likelihoods that will be used to estimate
     * allele frequencies.
     * @param emitGL the HMM emission probabilities.
     * @return a list of sampled haplotype pairs.
     *
     * @throws IllegalArgumentException if
     * {@code freqGL.markers().equals(emitGL.markers())==false}
     * @throws IllegalArgumentException if
     * {@code fam.samples().equals(emitGL.samples())==false}
     * @throws NullPointerException if any parameter is {@code null}
     */
    public List<HapPair> initialHaps(NuclearFamilies fam, GL freqGL, GL emitGL) {
        if (freqGL.markers().equals(emitGL.markers())==false) {
            throw new IllegalArgumentException("inconsistent markers");
        }
        if (fam.samples().equals(emitGL.samples())==false) {
            throw new IllegalArgumentException("inconsistent samples");
        }
        boolean useRevDag = false;
        float minAlleleFreq = 0.0001f;
        Dag dag = new LinkageEquilibriumDag(freqGL, minAlleleFreq);
        List<HapPair> sampledHaps = new ArrayList<>();
        sampledHaps = Collections.synchronizedList(sampledHaps);
        singleSample(fam, dag, emitGL, useRevDag, par.nsamples(),
                sampledHaps, par.nthreads());
        duoSample(fam, dag, emitGL, useRevDag, par.nsamples(), sampledHaps);
        trioSample(fam, dag, emitGL, useRevDag, par.nsamples(), sampledHaps);
        return new ArrayList<>(sampledHaps);
    }

    /**
     * Returns a list of sampled haplotype pairs.
     * @param useRevDag {@code true} if the order of markers should
     * be reversed when building the haplotype frequency model, and
     * {@code false} otherwise.
     * @param gl the HMM emission probabilities.
     * @param haps the haplotypes that will be used to build a haplotype
     * frequency model.
     * @param fam the parent-offspring relationships.
     * @param weights the per-haplotype weights.
     * @return a list of sampled haplotype pairs.
     *
     * @throws IllegalArgumentException if {@code haps.isEmpty()==true}
     * @throws IllegalArgumentException if
     * {@code gl.markers().equals(haps.markers(j))==false} for any
     * {@code j} satisfying {@code 0<=j && j<haps.size()}
     * @throws IllegalArgumentException if
     * {@code fam.samples().equals(gl.samples())==false}
     * @throws NullPointerException if any parameter is null
     */
    public List<HapPair> sample(boolean useRevDag, GL gl, List<HapPair> haps,
            NuclearFamilies fam, Weights weights) {
        if (fam.samples().equals(gl.samples())==false) {
            throw new IllegalArgumentException("inconsistent samples");
        }
        if (gl.markers().equals(hapsMarkers(haps))==false) {
            throw new IllegalArgumentException("inconsistent markers");
        }
        List<HapPair> sampledHaps = new ArrayList<>(haps.size());
        HapPairs dagHaps = new BasicHapPairs(haps, useRevDag);
        float[] wts = weights.get(dagHaps);
        Dag singleDag = singleDag(fam, dagHaps, wts);
        Dag duoDag = duoDag(fam, singleDag, dagHaps, wts);
        Dag trioDag = trioDag(fam, singleDag, duoDag, dagHaps, wts);
        setDagStats(singleDag, duoDag, trioDag);
        if (useRevDag) {
            gl = new RevGL(gl);
        }
        sampledHaps = Collections.synchronizedList(sampledHaps);
        singleSample(fam, singleDag, gl, useRevDag, par.nsamples(), sampledHaps,
                par.nthreads());
        duoSample(fam, duoDag, gl, useRevDag, par.nsamples(), sampledHaps);
        trioSample(fam, trioDag, gl, useRevDag, par.nsamples(), sampledHaps);
        return new ArrayList<>(sampledHaps);
    }

    /**
     * Returns a list of sampled haplotype pairs.
     * @param useRevDag {@code true} if the order of markers should
     * be reversed when building the haplotype frequency model, and
     * {@code false} otherwise.
     * @param gl the HMM emission probabilities.
     * @param haps the haplotypes that will be used to build a haplotype
     * frequency model.
     * @param gv the object in which posterior genotype probabilities will
     * be stored.
     * @param fam the parent-offspring relationships.
     * @param weights the per-haplotype weights.
     * @return a list of sampled haplotype pairs.
     *
     * @throws IllegalArgumentException if {@code haps.isEmpty()==true}
     * @throws IllegalArgumentException if
     * {@code fam.samples().equals(gl.samples())==false
                || fam.samples().equals(gv.samples())==false}
     * @throws IllegalArgumentException if
     * {@code gl.markers().equals(haps.markers(j))==false} for any
     * {@code j} satisfying {@code 0<=j && j<haps.size()}
     * @throws IllegalArgumentException if
     * {@code gl.markers().equals(gv.markers())==false}
     * @throws NullPointerException if any parameter is null
     */
     public List<HapPair> sample(boolean useRevDag, GL gl, List<HapPair> haps,
             GenotypeValues gv, NuclearFamilies fam, Weights weights) {
        if (fam.samples().equals(gl.samples())==false
                || fam.samples().equals(gv.samples())==false) {
            throw new IllegalArgumentException("inconsistent samples");
        }
        if (gl.markers().equals(hapsMarkers(haps))==false
                || gl.markers().equals(gv.markers())==false) {
            throw new IllegalArgumentException("inconsistent markers");
        }
        HapPairs dagHaps = new BasicHapPairs(haps, useRevDag);
        float[] wts = weights.get(dagHaps);
        Dag singleDag = singleDag(fam, dagHaps, wts);
        Dag duoDag = duoDag(fam, singleDag, dagHaps, wts);
        Dag trioDag = trioDag(fam, singleDag, duoDag, dagHaps, wts);
        setDagStats(singleDag, duoDag, trioDag);
        if (useRevDag) {
            gl = new RevGL(gl);
            gv = new RevGenotypeValues(gv);
        }

        List<HapPair> sampledHaps = new ArrayList<>(haps.size());
        sampledHaps = Collections.synchronizedList(sampledHaps);
        singleSample(fam, singleDag, gl, useRevDag,par.nsamples(), sampledHaps,
                gv, par.nthreads());
        duoSample(fam, duoDag, gl, useRevDag, par.nsamples(), sampledHaps, gv);
        trioSample(fam, trioDag, gl, useRevDag, par.nsamples(), sampledHaps, gv);
        return new ArrayList<>(sampledHaps);
    }

    /**
     * Returns a list of sampled haplotype pairs.
     * @param useRevDag {@code true} if the order of markers should
     * be reversed when building the haplotype frequency model, and
     * {@code false} otherwise.
     * @param al the HMM emission probabilities.
     * @param haps the haplotypes that will be used to build a haplotype
     * frequency model.
     * @param gv the object in which posterior genotype probabilities will
     * be stored.
     * @param weights the per-haplotype weights.
     * @return a list of sampled haplotype pairs.
     *
     * @throws IllegalArgumentException if {@code haps.isEmpty()==true}
     * @throws IllegalArgumentException if
     * {@code gl.markers().equals(haps.markers(j))==false} for any
     * {@code j} satisfying {@code 0<=j && j<haps.size()}
     * @throws IllegalArgumentException if
     * {@code gl.markers().equals(gv.markers())==false}
     * @throws NullPointerException if any parameter is null
     */
     public List<HapPair> sample(boolean useRevDag, AL al,
             List<HapPair> haps, GenotypeValues gv, Weights weights) {
        if (al.markers().equals(hapsMarkers(haps))==false
                || al.markers().equals(gv.markers())==false) {
            throw new IllegalArgumentException("inconsistent markers");
        }
        HapPairs dagHaps = new BasicHapPairs(haps, useRevDag);
        float[] wts = weights.get(dagHaps);
        Dag singleDag =  makeDag(dagHaps, wts, par.buildwindow(),
                par.singlescale());
        setDagStats(singleDag, null, null);
        if (useRevDag) {
            al = new RevAL(al);
            gv = new RevGenotypeValues(gv);
        }

        List<HapPair> sampledHaps = new ArrayList<>(haps.size());
        sampledHaps = Collections.synchronizedList(sampledHaps);
        hapSample(singleDag, al, useRevDag, par.nsamples(), sampledHaps, gv,
                par.nthreads());
        return new ArrayList<>(sampledHaps);
    }

    private Markers hapsMarkers(List<HapPair> haps) {
        if (haps.isEmpty()) {
            throw new IllegalArgumentException("haps.isEmpty()");
        }
        else {
            return haps.get(0).markers();
        }
    }

    private Dag singleDag(NuclearFamilies fam, HapPairs haps,
            float[] weights) {
        if (fam.nSingles()>0) {
            return makeDag(haps, weights, par.buildwindow(), par.singlescale());
        }
        else {
            return null;
        }
    }

    private Dag duoDag(NuclearFamilies fam,
            Dag singleDag, HapPairs haps, float[] weights) {
        if (fam.nDuos() > 0) {
            if (fam.nSingles()>0 && par.singlescale()==par.duoscale()) {
                return singleDag;
            }
            else {
                return makeDag(haps, weights, par.buildwindow(), par.duoscale());
            }
        }
        else {
            return null;
        }
    }

    private Dag trioDag(NuclearFamilies fam, Dag singleDag,
            Dag duoDag, HapPairs haps, float[] weights) {
        if (fam.nTrios() > 0) {
            if (fam.nSingles()>0 && par.singlescale()==par.trioscale()) {
                return singleDag;
            }
            else if (fam.nDuos()>0 && par.duoscale()==par.trioscale()) {
                return duoDag;
            }
            else {
                return makeDag(haps, weights, par.buildwindow(), par.trioscale());
            }
        }
        else {
            return null;
        }
    }

    private void setDagStats(Dag singleDag, Dag duoDag,
            Dag trioDag) {
        runStats.setSingleDagStats(singleDag);
        runStats.setDuoDagStats(duoDag);
        runStats.setTrioDagStats(trioDag);
    }

    @SuppressWarnings({"BroadCatchBlock", "TooBroadCatch"})
    private void singleSample(NuclearFamilies fam, Dag dag, GL gl,
            boolean markersAreReversed, int nSamples, List<HapPair> sampledHaps,
            int nThreads) {
        if (fam.nSingles() > 0) {
            long t0 = System.currentTimeMillis();
            Random rand = new Random(par.seed());
            final BlockingQueue<Integer> qIn = new ArrayBlockingQueue<>(3*nThreads);
            ExecutorService es = Executors.newFixedThreadPool(nThreads);
            for (int j=0; j<nThreads; ++j) {
                SingleBaum sb = new SingleBaum(dag, gl, rand.nextLong(), nSamples);
                es.submit(new ProduceSingleSamples(markersAreReversed, sb, qIn,
                        sampledHaps));
            }
            try {
                for (int j=0, n=fam.nSingles(); j<n; ++j) {
                    qIn.put(fam.single(j));
                }
                for (int j=0; j<nThreads; ++j) {
                    qIn.put(ProduceSingleSamples.POISON);
                }
                es.shutdown();
                es.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
            }
            catch (Throwable e) {
                Utilities.exit("\"HapSampler: ERROR\"", e);
            }
            runStats.singleSampleMillis(System.currentTimeMillis() - t0);
        }
    }

    private void duoSample(NuclearFamilies fam, Dag dag,
            GL gl, boolean markersAreReversed,
            int nSamples, List<HapPair> sampledHaps) {
        if (fam.nDuos() > 0) {
            long t0 = System.currentTimeMillis();
            DuoBaum duoBaum = new DuoBaum(dag, gl, par.seed(), nSamples);
            for (int j=0, n=fam.nDuos(); j<n; ++j) {
                List<HapPair> newHaps = duoBaum.sample(
                        fam.duoParent(j), fam.duoOffspring(j));
                storeHaps(sampledHaps, markersAreReversed, newHaps);
            }
            runStats.duoSampleMillis(System.currentTimeMillis() - t0);
        }
    }

    private void trioSample(NuclearFamilies fam, Dag dag,
            GL gl, boolean markersAreReversed,
            int nSamples, List<HapPair> sampledHaps) {
        if (fam.nTrios() > 0) {
            long t0 = System.currentTimeMillis();
            TrioBaum trioBaum = new TrioBaum(dag, gl, par.seed(), nSamples);
            for (int j=0, n=fam.nTrios(); j<n; ++j) {
                List<HapPair> newHaps = trioBaum.sample(
                        fam.trioFather(j), fam.trioMother(j), fam.trioOffspring(j));
                storeHaps(sampledHaps, markersAreReversed, newHaps);
            }
            runStats.trioSampleMillis(System.currentTimeMillis() - t0);
        }
    }

    @SuppressWarnings({"BroadCatchBlock", "TooBroadCatch"})
    private void hapSample(Dag dag, AL al, boolean markersAreReversed,
            int nCopies, List<HapPair> sampledHaps, GenotypeValues gv,
            int nThreads) {
        long t0 = System.currentTimeMillis();
        Random rand = new Random(par.seed());
        final BlockingQueue<Integer> qIn = new ArrayBlockingQueue<>(3*nThreads);
        ExecutorService es = Executors.newFixedThreadPool(nThreads);
        for (int j=0; j<nThreads; ++j) {
            HapBaum hb = new HapBaum(dag, al, rand.nextLong(),
                    nCopies);
            es.submit(new ProduceHapSamples(markersAreReversed, hb, qIn,
                    sampledHaps, gv));
        }
        try {
            for (int j=0, n=gv.nSamples(); j<n; ++j) {
                qIn.put(j);
            }
            for (int j=0; j<nThreads; ++j) {
                qIn.put(ProduceSingleSamples.POISON);
            }
            es.shutdown();
            es.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
        }
        catch (Throwable e) {
            Utilities.exit("\"HapSampler: ERROR\"", e);
        }
        runStats.singleSampleMillis(System.currentTimeMillis() - t0);
    }

    @SuppressWarnings({"BroadCatchBlock", "TooBroadCatch"})
    private void singleSample(NuclearFamilies fam, Dag dag,
            GL gl, boolean markersAreReversed, int nCopies,
            List<HapPair> sampledHaps, GenotypeValues gv,
            int nThreads) {
        if (fam.nSingles() > 0) {
            long t0 = System.currentTimeMillis();
            Random rand = new Random(par.seed());
            final BlockingQueue<Integer> qIn = new ArrayBlockingQueue<>(3*nThreads);
            ExecutorService es = Executors.newFixedThreadPool(nThreads);
            for (int j=0; j<nThreads; ++j) {
                SingleBaum sb = new SingleBaum(dag, gl, rand.nextLong(), nCopies);
                es.submit(new ProduceSingleSamples(markersAreReversed, sb, qIn,
                        sampledHaps, gv));
            }
            try {
                for (int j=0, n=fam.nSingles(); j<n; ++j) {
                    qIn.put(fam.single(j));
                }
                for (int j=0; j<nThreads; ++j) {
                    qIn.put(ProduceSingleSamples.POISON);
                }
                es.shutdown();
                es.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
            }
            catch (Throwable e) {
                Utilities.exit("\"HapSampler: ERROR\"", e);
            }
            runStats.singleSampleMillis(System.currentTimeMillis() - t0);
        }
    }

    private void duoSample(NuclearFamilies fam, Dag dag,
            GL curEm, boolean markersAreReversed, int nCopies,
            List<HapPair> sampledHaps, GenotypeValues gv) {
        if (fam.nDuos() > 0) {
            long t0 = System.currentTimeMillis();
            int gprobsLength = curEm.markers().sumGenotypes();
            double[] gprobsA = new double[gprobsLength];
            double[] gprobsB = new double[gprobsLength];
            DuoBaum duoBaum = new DuoBaum(dag, curEm, par.seed(), nCopies);
            for (int j=0, n=fam.nDuos(); j<n; ++j) {
                List<HapPair> newHaps = duoBaum.sample(
                        fam.duoParent(j), fam.duoOffspring(j), gprobsA, gprobsB);
                storeHaps(sampledHaps, markersAreReversed, newHaps);
                gv.add(fam.duoParent(j), gprobsA);
                gv.add(fam.duoOffspring(j), gprobsB);
            }
            runStats.duoSampleMillis(System.currentTimeMillis() - t0);
        }
    }

    private void trioSample(NuclearFamilies fam, Dag dag,
            GL curEm, boolean markersAreReversed, int nCopies,
            List<HapPair> sampledHaps, GenotypeValues gv) {
        if (fam.nTrios() > 0) {
            long t0 = System.currentTimeMillis();
            int gprobsLength = curEm.markers().sumGenotypes();
            double[] gprobsA = new double[gprobsLength];
            double[] gprobsB = new double[gprobsLength];
            double[] gprobsC = new double[gprobsLength];
            TrioBaum trioBaum = new TrioBaum(dag, curEm, par.seed(), nCopies);
            for (int j=0, n=fam.nTrios(); j<n; ++j) {
                List<HapPair> newHaps = trioBaum.sample(
                        fam.trioFather(j), fam.trioMother(j), fam.trioOffspring(j),
                        gprobsA, gprobsB, gprobsC);
                storeHaps(sampledHaps, markersAreReversed, newHaps);
                gv.add(fam.trioFather(j), gprobsA);
                gv.add(fam.trioMother(j), gprobsB);
                gv.add(fam.trioOffspring(j), gprobsC);
            }
            runStats.trioSampleMillis(System.currentTimeMillis() - t0);
        }
    }

    private void storeHaps(List<HapPair> sampledHaps,
            boolean reverseDag, List<HapPair> newHaps) {
        if (reverseDag) {
            for (HapPair hp : newHaps) {
                sampledHaps.add(new RevHapPair(hp));
            }
        }
        else {
            sampledHaps.addAll(newHaps);
        }
    }

    private Dag makeDag(HapPairs haps, float[] weights, int window, float scale) {
        long t0 = System.currentTimeMillis();
        Dag dag = MergeableDag.dag(haps, weights, window, scale);
        runStats.buildMillis(System.currentTimeMillis()-t0);
        return dag;
    }
}
