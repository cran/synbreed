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

import blbutil.Const;
import blbutil.IntPair;
import dag.Dag;
import dag.MergeableDag;
import haplotype.BasicHapPairs;
import haplotype.BasicSampleHapPairs;
import haplotype.ConsensusPhaser;
import haplotype.GLSampleHapPairs;
import haplotype.GenotypeCorrection;
import haplotype.HapPair;
import haplotype.HapPairs;
import haplotype.SampleHapPairs;
import haplotype.WrappedHapPair;
import ibd.HaploidIbd;
import ibd.IbdSegment;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import vcf.GL;
import vcf.MaskedEndsGL;
import vcf.NoPhaseGL;

/**
 * Class {@code MainHelper} is an auxiliary class with methods called by
 * the {@code main.Main} class.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class MainHelper {

    private final Par par;
    private final HapPairSampler hapSampler;
    private final RecombHapPairSampler recombSampler;
    private final GeneticMap genMap;
    private final RunStats runStats;

    /**
     * Constructs a new {@code MainHelper} instance.
     * @param par the command line parameters
     * @param genMap the genetic map or {@code null} if no genetic map
     * was specified
     * @param runStats the class for collecting and printing run-time statistics
     * @throws NullPointerException
     * if {@code (par == null || runStarts == null)}
     */
    MainHelper(Par par, GeneticMap genMap, RunStats runStats) {
        if (runStats==null) {
            throw new NullPointerException("runStats==null");
        }
        if (genMap==null) {
            double scaleFactor = 1e-6;
            genMap = new PositionMap(scaleFactor);
        }
        this.par = par;
        this.hapSampler = new HapPairSampler(par, runStats);
        this.recombSampler = new RecombHapPairSampler(par, runStats);
        this.genMap = genMap;
        this.runStats = runStats;
    }

    /**
     * Phases the current window of genotype data.
     * @param cd the current window of data
     * @param gv the current scaled genotype probabilities for the target
     * samples or {@code null} if genotype probabilities are not to be estimated
     * @return the phased genotype data
     * @throws IllegalArgumentException if
     * {@code gv != null && gv.markers().equals(cd.targetMarkers() == false)}
     * @throws IllegalArgumentException if
     * {@code gv != null && gv.samples().equals(cd.targetSamples() == false)}
     * @throws NullPointerException if {@code cd == null}
     */
    SampleHapPairs phase(CurrentData cd, GenotypeValues gv) {
        checkParameters(cd, gv);
        if (cd.targetGL().isRefData()) {
            return new GLSampleHapPairs(cd.targetGL());
        }

        List<HapPair> hapPairs = hapSampler.initialHaps(cd);
        if (par.burnin_its()>0) {
            runStats.println(Const.nl + "Starting burn-in iterations");
            hapPairs = runBurnin1(cd, hapPairs);
        }
        if (par.phase_its()>0) {
            boolean estGprobs = (par.gt()==null && par.niterations()==0);
            hapPairs = runBurnin2(cd, hapPairs, (estGprobs ? gv : null));
        }
        if (par.niterations()>0) {
            runStats.println(Const.nl + "Starting phasing iterations");
            hapPairs = runRecomb(cd, hapPairs, gv);
        }
        else {
            hapPairs = ConsensusPhaser.run(hapPairs);
        }
        return new BasicSampleHapPairs(cd.targetSamples(), hapPairs);
    }

    private void checkParameters(CurrentData cd, GenotypeValues gv) {
        if (gv!=null && gv.markers().equals(cd.targetMarkers())==false) {
            throw new IllegalArgumentException(String.valueOf(gv));
        }
        if (gv!=null && gv.samples().equals(cd.targetSamples())==false) {
            throw new IllegalArgumentException(String.valueOf(gv));
        }
    }

    private List<HapPair> runBurnin1(CurrentData cd, List<HapPair> hapPairs) {
        GenotypeValues gv = null;
        for (int j=0; j<par.burnin_its(); ++j) {
            boolean useRevDag = (j & 1)==1;
            hapPairs = hapSampler.sample(cd, hapPairs, useRevDag, gv);
            runStats.printIterationUpdate(cd.window(), j+1);
        }
        return hapPairs;
    }

    private List<HapPair> runBurnin2(CurrentData cd, List<HapPair> hapPairs,
            GenotypeValues gv) {
        List<HapPair> cumHapPairs = new ArrayList<>();
        int start = par.burnin_its();
        int end = start + par.phase_its();
        for (int j=start; j<end; ++j) {
            boolean useRevDag = (j & 1)==1;
            hapPairs = hapSampler.sample(cd, hapPairs, useRevDag, gv);
            runStats.printIterationUpdate(cd.window(), j+1);
            cumHapPairs.addAll(hapPairs);
        }
        return cumHapPairs;
    }

    private List<HapPair> runRecomb(CurrentData cd, List<HapPair> hapPairs,
            GenotypeValues gv) {
        hapPairs = ConsensusPhaser.run(hapPairs);
        List<HapPair> cumHapPairs = new ArrayList<>();
        int start = par.burnin_its() + par.phase_its();
        int end = start + par.niterations();
        for (int j=start; j<end; ++j) {
            boolean useRevDag = (j & 1)==1;
            hapPairs = recombSampler.sample(cd, hapPairs, useRevDag, gv);
            runStats.printIterationUpdate(cd.window(), j+1);
            cumHapPairs.addAll(hapPairs);
        }
        hapPairs = ConsensusPhaser.run(cumHapPairs);
        hapPairs = correctGenotypes(cd, hapPairs);
        return hapPairs;
    }

    private List<HapPair> correctGenotypes(CurrentData cd, List<HapPair> hapPairs) {
        int start = cd.prevTargetSpliceStart();
        int end = cd.nextTargetSpliceStart();
        GL modGL = new MaskedEndsGL(cd.targetGL(), start, end);
//        File outFile = new File(par.out() + ".gterr");
//        boolean append = cd.window() > 1;
//        GenotypeCorrection.run(hapPairs, modGL, par.seed(), outFile, append);
        GenotypeCorrection.run(hapPairs, modGL, par.seed());
        return hapPairs;
    }

    /**
     * Applies the refined IBD algorithm to the specified data.
     * @param cd the current window of data
     * @param targetHapPairs the estimated haplotype pairs
     * @return the detected IBD segments
     * @throws NullPointerException if
     * {@code cd == null || targetHapPairs ==  null}
     */
    Map<IntPair, List<IbdSegment>> refinedIbd(CurrentData cd,
            SampleHapPairs targetHapPairs) {
        if (par.ibd()) {
            long t0 = System.nanoTime();
            int nSamples = cd.nRefSamples() + cd.nTargetSamples();
            float scale = par.adjustedIbdScale(nSamples);

            Dag dag = ibdDag(cd, targetHapPairs, scale);
            HaploidIbd hapIbd = new HaploidIbd(genMap, par.ibdtrim(),
                    par.ibdlod(), par.ibdcm());
            GL ibdGL = new NoPhaseGL(cd.targetGL());

            Map<IntPair, List<IbdSegment>> ibdMap =
                    hapIbd.run(ibdGL, dag, targetHapPairs, par.nthreads());
            long nanos = (System.nanoTime() - t0);
            runStats.ibdNanos(nanos);
            runStats.printRefinedIbdUpdate(scale, dag, nanos);
            return ibdMap;
        }
        else {
            return null;
        }
    }

    private Dag ibdDag(CurrentData cd, SampleHapPairs targetHaps,
            float scale) {
        float[] weights = cd.weights().get(targetHaps);
        float[] combWeights;
        HapPairs dagHaps;
        if (cd.nRefSamples() == 0) {
            dagHaps = targetHaps;
            combWeights = weights;
        }
        else {
            List<HapPair> hapsList = new ArrayList<>(cd.nAllSamples());
            cd.addRestrictedRefHapPairs(hapsList);
            for (int j=0, n=targetHaps.nSamples(); j<n; ++j) {
                hapsList.add(new WrappedHapPair(targetHaps, j));
            }
            dagHaps = new BasicHapPairs(hapsList);
            int nRefHaps = 2*cd.nRefSamples();
            combWeights = new float[dagHaps.nHaps()];
            Arrays.fill(combWeights, 0, nRefHaps, 1.0f);
            System.arraycopy(weights, 0, combWeights, nRefHaps, weights.length);
        }
        long t0 = System.nanoTime();
        int nInitLevels = 500;
        Dag ibdDag = MergeableDag.dag(dagHaps, combWeights, scale, nInitLevels);
        runStats.buildNanos(System.nanoTime() - t0);
        runStats.setDagStats(ibdDag);
        return ibdDag;
    }

    /**
     * Performs genotype imputation
     * @param cd the current window of data
     * @param shp the estimated target haplotype pairs.
     * @return imputed haplotypes
     * @throws NullPointerException if {@code cd == null || shp == null}
     */
    AlleleProbs LSImpute(CurrentData cd, SampleHapPairs shp) {
        if (cd.nMarkers()==cd.nTargetMarkers() || par.impute() == false) {
            return new SampleHapPairAlleleProbs(shp);
        }
        long t0 = System.nanoTime();
        LiAndStephensHapSampler recombHapSampler =
                new LiAndStephensHapSampler(par, genMap);

        BasicAlleleProbs alProbs = recombHapSampler.sample(cd, shp);
        runStats.imputationNanos(System.nanoTime() - t0);
        runStats.printImputationUpdate();
        return alProbs;
    }
}
