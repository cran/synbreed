/*
 * Copyright 2015 Brian L. Browning
 *
 * This file is part of Beagle
 *
 * IBD is licensed under the Apache License, Version 2.0 (the "License");
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
import haplotype.BasicSampleHapPairs;
import haplotype.ConsensusPhaser;
import haplotype.HapPair;
import haplotype.RevSampleHapPairs;
import haplotype.SampleHapPairs;
import haplotype.Weights;
import java.util.ArrayList;
import java.util.List;
import main.CurrentData;
import main.Par;
import main.RunStats;
import vcf.FuzzyGL;
import vcf.GL;
import vcf.Markers;
import vcf.RevGL;

/**
 * <p>Class {@code SamplerData} contains immutable input data for the
 * current marker window.
 * </p>
 * <p>Instances of class {@code SamplerData} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class SamplerData {

    private static final double MIN_CM_DIST = 1e-7;
    private static final int nInitLevels = 500;

    private final Par par;
    private final boolean revMarkers;
    private final RestrictedDag rdag;
    private final GL gl;
    private final float[] recombRate;

    /**
     * Constructs a new {@code SamplerData} instance from the specified data.
     * The contract for this method is undefined if the specified
     * {@code hapPairs} is inconsistent with the input data
     * contained in the {@code cd} parameter.
     *
     * @param par the analysis parameters
     * @param cd the input data for the current marker window
     * @param hapPairs the target haplotype pairs used to build the haplotype
     * frequency model
     * @param revMarkers {@code true} if the order of markers should
     * be reversed when building the haplotype frequency model, and
     * {@code false} otherwise
     * @param runStats the object to which run-time statistics will be written
     *
     * @throws IllegalArgumentException if {@code haps.isEmpty() == true}
     * @throws NullPointerException if any parameter is {@code null}
     */
    public SamplerData(Par par, CurrentData cd, List<HapPair> hapPairs,
            boolean revMarkers, RunStats runStats) {
        if (hapPairs.isEmpty()) {
            throw new IllegalArgumentException("hapPairs.isEmpty()");
        }
        this.par = par;
        this.revMarkers = revMarkers;
        this.rdag = restrictedDag(cd, hapPairs, revMarkers, runStats);
        this.gl = gl(cd, par.err(), revMarkers);
        this.recombRate = recombRate(par, cd, rdag.dag(), revMarkers);
    }

    private static float[] recombRate(Par par, CurrentData cd, Dag dag,
            boolean revMarkers) {
        float[] recombRate = cd.recombRate();
        if (recombRate != null && revMarkers) {
            for (int j=1, n=(recombRate.length + 1)/2; j<n; ++j) {
                int k = recombRate.length - j;
                float tmp = recombRate[j];
                recombRate[j] = recombRate[k];
                recombRate[k] = tmp;
            }
            recombRate[0] = 0;
        }
        else {
            recombRate = dagRecombRate(dag, par.mapscale());
        }
        return recombRate;
    }

    private RestrictedDag restrictedDag(CurrentData cd, List<HapPair> hapPairs,
            boolean revMarkers, RunStats runStats) {
        hapPairs = new ArrayList<>(hapPairs);   // xx defensive copy
        long t0 = System.nanoTime();
        Weights weights = cd.weights();
        List<HapPair> haps = ConsensusPhaser.run(hapPairs);
        cd.addRestrictedRefHapPairs(haps);

        SampleHapPairs dagHaps = new BasicSampleHapPairs(cd.allSamples(), haps);
        if (revMarkers) {
            dagHaps = new RevSampleHapPairs(dagHaps);
        }
        float[] wts = weights.get(dagHaps);
        RestrictedDag rdag = new RestrictedDag(dagHaps, wts, nInitLevels,
                par.modelscale(), par.ibdlength(), par.ibdextend());
        runStats.buildNanos(System.nanoTime() - t0);
        runStats.setDagStats(rdag.dag());
        return rdag;
    }

    private static float[] dagRecombRate(Dag dag, float xdist) {
        double[] bglDist = dag.posArray();
        for (int j=0; j<bglDist.length; ++j) {
            bglDist[j] *= 0.2;
        }
        double c = -2.0*xdist;
        float[] rr = new float[dag.nLevels()];
        rr[0] = 0.0f;
        double lastGenPos = bglDist[0];
        for (int j=1; j<rr.length; ++j) {
            double genPos = bglDist[j];
            double genDist = Math.max(Math.abs(genPos - lastGenPos), MIN_CM_DIST);
            rr[j] = (float) -Math.expm1(c*genDist);
            lastGenPos = genPos;
        }
        return rr;
    }

    private static GL gl(CurrentData cd, float err, boolean markersAreReversed) {
        GL gl = new FuzzyGL(cd.targetGL(), err);
        if (markersAreReversed) {
            gl = new RevGL(gl);
        }
        return gl;
    }

    /**
     * Returns {@code true} if the order of markers is reversed, and
     * {@code false} otherwise
     * @return {@code true} if the order of markers is reversed, and
     * {@code false} otherwise
     */
    public boolean markersAreReversed() {
        return revMarkers;
    }

    /**
     * Returns the number of markers.
     * @return the number of markers
     */
    public int nMarkers() {
        return gl.nMarkers();
    }

    /**
     * Returns the number of samples.
     * @return the number of samples
     */
    public int nSamples() {
        return gl.nSamples();
    }

    /**
     * Returns the number of haplotypes.
     * @return the number of haplotypes
     */
    public int nHaps() {
        return 2*gl.nSamples();
    }

    /**
     * returns the list of markers.
     * @return the list of markers
     */
    public Markers markers() {
        return gl.markers();
    }

    /**
     * Returns the analysis parameters.
     * @return the analysis parameters
     */
    public Par par() {
        return par;
    }

    /**
     * Returns the DAG model.
     * @return the DAG model
     */
    public RestrictedDag rdag() {
        return rdag;
    }

    /**
     * Returns the genotype likelihoods for the
     * target samples at the target data markers.
     * @return the genotype likelihoods for the
     * target samples at the target data markers.
     */
    public GL gl() {
        return gl;
    }

    /**
     * Returns the allele error rate
     * @return the allele error rate
     */
    public float err() {
        return par.err();
    }

    /**
     * Returns the probability of recombination between {@code (marker - 1)}
     * and {@code marker}.
     * @param marker a marker index
     * @return the probability of recombination between {@code (marker - 1)}
     * and {@code marker}
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}
     */
    public float pRecomb(int marker) {
        return recombRate[marker];
    }
}

