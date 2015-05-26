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

import beagleutil.Samples;
import blbutil.Const;
import blbutil.IntPair;
import blbutil.Utilities;
import dag.Dag;
import dag.MergeableDag;
import haplotype.BasicHapPairs;
import haplotype.BasicSampleHapPairs;
import haplotype.ConsensusPhasing;
import haplotype.HapPair;
import haplotype.HapPairs;
import haplotype.SampleHapPairs;
import haplotype.Weights;
import haplotype.WrappedHapPair;
import ibd.HaploidIbd;
import ibd.IbdSegment;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Random;
import vcf.AL;
import vcf.Data;
import vcf.GL;
import vcf.HapAL;
import vcf.ImputationGL;
import vcf.Markers;
import vcf.NoPhaseGL;

/**
 * Class {@code MainHelper} is an auxiliary class with methods called by
 * the {@code main.Main} class.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
class MainHelper {

    private final Parameters par;
    private final HapPairSampler hapSampler;
    private final NuclearFamilies fam;
    private final Weights weights;
    private final RunStats runStats;

    MainHelper(Parameters par,  GeneticMap genMap, NuclearFamilies fam,
            Weights weights, RunStats runStats, Random random) {
        if (runStats==null) {
            throw new NullPointerException("runStats==null");
        }
        this.par = par;
        this.hapSampler = new HapPairSampler(par, runStats);
        this.fam = fam;
        this.weights = weights;
        this.runStats = runStats;
    }

    SampleHapPairs sample(Data data, GenotypeValues restrictedGV) {
        if (par.burnin_its()>0) {
            runStats.println(Const.nl + "Starting burn-in iterations");
        }
        List<HapPair> sampledHaps = runBurnin1(data);

        if (par.phase_its()>0) {
            runStats.println(Const.nl + "Starting phasing iterations");
        }
        sampledHaps = runBurnin2(data, sampledHaps, restrictedGV);

        return merge(data.nonRefSamples(), sampledHaps);
    }

    private List<HapPair> runBurnin1(Data data) {
        checkPrephasedTarget(data.nonRefEmissions().isRefData());
        int startIt = 0;
        int endIt = par.burnin_its();
        boolean useRevDag = (startIt % 2)==0;
        List<HapPair> restrictedRefHaps = data.restrictedRefHaps();
        GL gl = data.nonRefEmissions();
        List<HapPair> sampledHaps = hapSampler.initialHaps(fam, gl, gl);
        for (int iter=startIt; iter<endIt; ++iter) {
            useRevDag = !useRevDag;
            sampledHaps.addAll(restrictedRefHaps);
            sampledHaps = hapSampler.sample(useRevDag, gl, sampledHaps, fam,
                    weights);
            runStats.printIterationUpdate(data.window(), iter+1);
        }
        return sampledHaps;
    }

    private List<HapPair> runBurnin2(Data data, List<HapPair> sampledHaps,
            GenotypeValues gv) {
        if (par.phase_its()==0) {
            return sampledHaps;
        }
        int startIt = par.burnin_its();
        int endIt = startIt + par.phase_its();
        boolean useRevDag = (startIt % 2)==0;
        List<HapPair> allSamples = new ArrayList<>(par.nsamples()*par.phase_its());
        List<HapPair> restrictedRefHaps = data.restrictedRefHaps();
        GL gl = data.nonRefEmissions();
        for (int iter=startIt; iter<endIt; ++iter) {
            useRevDag = !useRevDag;
            sampledHaps.addAll(restrictedRefHaps);
            if (gv==null) {
                sampledHaps = hapSampler.sample(useRevDag, gl, sampledHaps, fam,
                        weights);
            }
            else {
                sampledHaps = hapSampler.sample(useRevDag, gl, sampledHaps, gv,
                        fam, weights);
            }
            allSamples.addAll(sampledHaps);
            runStats.printIterationUpdate(data.window(), iter+1);
        }
        return allSamples;
    }

    Map<IntPair, List<IbdSegment>> refinedIbd(List<HapPair> refHaps,
            GL gl, SampleHapPairs nextHaps, Weights weights) {
        if (par.ibd()) {
            long time = System.nanoTime();
            int nSamples = refHaps.size() + nextHaps.nSamples();
            float scale = par.adjustedIbdScale(nSamples);
            float[] wts = weights.get(nextHaps);
            Dag dag = ibdDag(refHaps, nextHaps, wts, scale, par.buildwindow());
            HaploidIbd hapIbd = new HaploidIbd(par.ibdtrim(), par.ibdlod());
            GL ibdGL = new NoPhaseGL(gl);
            Map<IntPair, List<IbdSegment>> ibdMap;

            ibdMap = hapIbd.run(ibdGL, dag, nextHaps, par.nthreads());
            long millis = (System.nanoTime() - time)/Const.mega;
            runStats.ibdMillis(millis);
            runStats.printRefinedIbdUpdate(scale, dag, millis);
            return ibdMap;
        }
        else {
            return null;
        }
    }

    private Dag ibdDag(List<HapPair> refHaps, SampleHapPairs targetHaps,
            float[] weights, float ibdScale, int buildWindow) {
        float[] combWeights;
        HapPairs dagHaps;
        if (refHaps.isEmpty()) {
            dagHaps = targetHaps;
            combWeights = weights;
        }
        else {
            List<HapPair> hapsList = new ArrayList<>(refHaps);
            for (int j=0, n=targetHaps.nSamples(); j<n; ++j) {
                hapsList.add(new WrappedHapPair(targetHaps, j));
            }
            dagHaps = new BasicHapPairs(hapsList);
            int nRefHaps = 2*refHaps.size();
            combWeights = new float[dagHaps.nHaps()];
            Arrays.fill(combWeights, 0, nRefHaps, 1.0f);
            System.arraycopy(weights, 0, combWeights, nRefHaps, weights.length);
        }
        long t0 = System.currentTimeMillis();
        Dag ibdDag = MergeableDag.dag(dagHaps, combWeights, buildWindow, ibdScale);
        runStats.buildMillis(System.currentTimeMillis()-t0);
        runStats.setSingleDagStats(ibdDag);
        return ibdDag;
    }

    SampleHapPairs impute(Data data,  SampleHapPairs mergedHaps,
            GenotypeValues gv) {
        Markers markers = data.markers();
        GL refEmissions = data.refEmissions();
        if (markers.nMarkers() > mergedHaps.nMarkers() ) {
            int startIt = par.burnin_its() + par.phase_its();
            int endIt = startIt + par.impute_its();
            int size = par.nsamples()*par.impute_its();
            List<HapPair> allSamples = new ArrayList<>(size);
            float nonRefWt = Math.min(1.0f, 10.0f/gv.nSamples());
            NuclearFamilies noFams = new NuclearFamilies(gv.samples(), null); // remove Mendelian constraints
            GenotypeValues fixedGv = new FixedGenotypeValues(gv, mergedHaps.markers());
            float err = 0.0f;
            AL al = new HapAL(gv.markers(), mergedHaps, err);
            GL gl = new ImputationGL(markers,  mergedHaps);
            List<HapPair> sampledHaps = hapSampler.initialHaps(noFams, refEmissions, gl);

            runStats.println(Const.nl + "Starting imputation iterations");
            Weights imputeWeights = new Weights(noFams, nonRefWt);
            runImpIts(data, sampledHaps, al, startIt, endIt, imputeWeights,
                    allSamples, fixedGv);
            mergedHaps = merge(gv.samples(), allSamples);
        }
        return mergedHaps;
    }

    private List<HapPair> runImpIts(Data data, List<HapPair> modelHaps, AL al,
            int startIt, int endIt, Weights imputeWeights,
            List<HapPair> allSamples, GenotypeValues gv) {
        List<HapPair> refHaps = data.refHaps();
        boolean useRevDag = (startIt % 2)==0;
        for (int iter=startIt; iter<endIt; ++iter) {
            useRevDag = !useRevDag;
            modelHaps.addAll(refHaps);
            modelHaps = hapSampler.sample(useRevDag, al, modelHaps, gv,
                    imputeWeights);
            if (allSamples!=null) {
                allSamples.addAll(modelHaps);
            }
            runStats.printIterationUpdate(data.window(), iter+1);
        }
        return modelHaps;
    }

    static SampleHapPairs merge(Samples samples,
            List<HapPair> hapPairList) {
        List<HapPair> mergedList
                = ConsensusPhasing.consensusHaps(hapPairList);
        Collections.sort(mergedList, BasicSampleHapPairs.hapsComparator(samples));
        return new BasicSampleHapPairs(samples, mergedList);
    }

    private void checkPrephasedTarget(boolean targetIsRef) {
        if (par.burnin_its()==0 && par.phase_its()==0 && targetIsRef==false) {
            String s = "ERROR: If \"burnin-its=0 phase-its=0\", all target genotypes  "
                    + "are required to be" + Const.nl
                    + "       phased, to have the '|'allele separator,"
                    + " and to have no missing alleles.";
            Utilities.exit(Parameters.usage() + s);
        }
    }
}
