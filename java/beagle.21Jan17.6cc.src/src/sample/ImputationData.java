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
package sample;

import beagleutil.Samples;
import blbutil.IntArray;
import haplotype.SampleHapPairs;
import java.util.Arrays;
import main.CurrentData;
import main.GeneticMap;
import main.Par;
import vcf.Markers;

/**
 * <p>Class {@code ImputationData} contains the input data that is
 * required for imputation of ungenotyped markers in the imputation target.
 * </p>
 * <p>Instances of class {@code ImputationData} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class ImputationData {

    private static final double MIN_CM_DIST = 1e-7;

    private final SampleHapPairs refHapPairs;
    private final SampleHapPairs targHapPairs;
    private final RefHapSegs refHapSegs;
    private final IntArray[] refAlleles;
    private final IntArray[] targAlleles;
    private final float[] errProb;
    private final float[] pRecomb;
    private final float[] weight;
    private final int nClusters;

    /**
     * Constructs a new {@code ImputationData} instance from the specified data.
     * @param par the analysis parameters
     * @param cd the reference haplotype data for the current marker window
     * @param targetHapPairs the target haplotype pairs
     * @param map the genetic map
     *
     * @throws IllegalArgumentException if
     * {@code cd.targetMarkers().equals(targetHapPairs.markers() == false}
     * @throws IllegalArgumentException if
     * {@code cd.targetSamples().equals(targetHapPairs.samples()) == false}
     * @throws NullPointerException if any parameter is {@code null}
     */
    public ImputationData(Par par, CurrentData cd,
            SampleHapPairs targetHapPairs, GeneticMap map) {
        if (cd.targetMarkers().equals(targetHapPairs.markers())==false) {
            throw new IllegalArgumentException("inconsistent markers");
        }
        if (cd.targetSamples().equals(targetHapPairs.samples())==false) {
            throw new IllegalArgumentException("inconsistent samples");
        }
        int[] targClustEnd = targClustEnd(targetHapPairs.markers(), map, par.cluster());
        this.refAlleles = new IntArray[targClustEnd.length];
        this.targAlleles = new IntArray[targClustEnd.length];
        setCodedAlleles(cd.restrictedRefSampleHapPairs(), targetHapPairs,
                targClustEnd, refAlleles, targAlleles);
        this.nClusters = targClustEnd.length;
        this.refHapPairs = cd.refSampleHapPairs();
        this.refHapSegs = refHapSegs(refHapPairs, targClustEnd, cd.markerIndices());
        this.targHapPairs = targetHapPairs;
        this.errProb = err(par.err(), targClustEnd);
        this.pRecomb = ImputationData.pRecomb(refHapSegs, map, par.ne());
        this.weight = wts(refHapSegs, map);
    }

    private static int[] targClustEnd(Markers targetMarkers, GeneticMap genMap,
            float clusterDist) {
        int nMarkers = targetMarkers.nMarkers();
        int[] ends = new int[nMarkers];
        double startPos = genMap.genPos(targetMarkers.marker(0));
        int index = 0;
        for (int m=1; m<nMarkers; ++m) {
            double pos = genMap.genPos(targetMarkers.marker(m));
            if ((pos - startPos) > clusterDist)  {
                ends[index++] = m;
                startPos = pos;
            }
        }
        ends[index++] = nMarkers;
        return Arrays.copyOf(ends, index);
    }

    private static void setCodedAlleles(SampleHapPairs refHapPairs,
            SampleHapPairs targetHapPairs, int[] targEnd,
            IntArray[] refAlleles, IntArray[] targAlleles) {
        HaplotypeCoder coder = new HaplotypeCoder(refHapPairs, targetHapPairs);
        int start = 0;
        for (int j=0; j<targEnd.length; ++j) {
            IntArray[] ia = coder.run(start, targEnd[j]);
            refAlleles[j] = ia[0];
            targAlleles[j] = ia[1];
            start = targEnd[j];
        }
    }

    private static float[] err(float errRate, int[] targClustEnd) {
        float maxErrProb = 0.5f;
        float[] err = new float[targClustEnd.length];
        int start = 0;
        for (int j=0; j<err.length; ++j) {
            err[j] = errRate * (targClustEnd[j] - start);
            if (err[j] > maxErrProb) {
                err[j] = maxErrProb;
            }
            start = targClustEnd[j];
        }
        return err;
    }

    private static float[] pRecomb(RefHapSegs refHapSegs, GeneticMap map,
            float ne) {
        SampleHapPairs refHaps = refHapSegs.refHapPairs();
        Markers refMarkers = refHaps.markers();
        int nHaps = refHaps.nHaps();
        int[] midPos = midPos(refMarkers, refHapSegs);
        int chrom = refMarkers.marker(0).chromIndex();
        return pRcomb(chrom, midPos, nHaps, map, ne);
    }

    private static int[] midPos(Markers refMarkers, RefHapSegs refHapSegs) {
        int[] midPos = new int[refHapSegs.nSegs() - 1];
        for (int j=0; j<midPos.length; ++j) {
            int startPos = refMarkers.marker(refHapSegs.segStart(j+1)).pos();
            int endPos = refMarkers.marker(refHapSegs.segEnd(j) - 1).pos();
            midPos[j] = (startPos + endPos) / 2;
        }
        return midPos;
    }

    private static float[] pRcomb(int chrom, int[] midPos, int nHaps,
            GeneticMap map, float ne) {
        float[] rr = new float[midPos.length];
        double c = -(0.04*ne/nHaps);    // 0.04 = 4/(100 cM/M)
        double lastGenPos = map.genPos(chrom, midPos[0]);
        rr[0] = 0f;
        for (int j=1; j<rr.length; ++j) {
            double genPos = map.genPos(chrom, midPos[j]);
            double genDist = Math.max(Math.abs(genPos - lastGenPos), MIN_CM_DIST);
            rr[j] = (float) -Math.expm1(c*genDist);
            lastGenPos = genPos;
        }
        return rr;
    }

    private static float[] wts(RefHapSegs refHapSegs, GeneticMap map) {
        Markers refMarkers = refHapSegs.refHapPairs().markers();
        double[] cumPos = cumPos(refMarkers, map);
        int nMarkers = refMarkers.nMarkers();
        int nClusters = refHapSegs.nSegs() - 1;
        float[] wts = new float[cumPos.length];
        if (nClusters > 0) {
            Arrays.fill(wts, 0, refHapSegs.segStart(1), Float.NaN);
        }
        for (int j=1; j<nClusters; ++j) {
            int start = refHapSegs.segStart(j);
            int end = refHapSegs.segEnd(j-1);
            int nextStart = refHapSegs.segStart(j+1);
            double nextStartPos = cumPos[nextStart];
            double totalLength = nextStartPos - cumPos[end - 1];
            Arrays.fill(wts, start, end, 1f);
            for (int m=end; m<nextStart; ++m) {
                wts[m] = (float) ( (nextStartPos - cumPos[m]) / totalLength );
            }
        }
        Arrays.fill(wts, refHapSegs.segStart(nClusters), nMarkers, Float.NaN);
        return wts;
    }

    private static double[] cumPos(Markers markers, GeneticMap map) {
        double[] cumPos = new double[markers.nMarkers()];
        double lastGenPos = map.genPos(markers.marker(0));
        cumPos[0] = 0.0;
        for (int j=1; j<cumPos.length; ++j) {
            double genPos = map.genPos(markers.marker(j));
            double genDist = Math.max(Math.abs(genPos - lastGenPos), MIN_CM_DIST);
            cumPos[j] = cumPos[j-1] + genDist;
            lastGenPos = genPos;
        }
        return cumPos;
    }

    private static RefHapSegs refHapSegs(SampleHapPairs refHapPairs,
            int[] targClustEnd, int[] targToRef) {
        int n = targClustEnd.length;
        int[] refClusterStart = new int[n+1];
        int[] refClusterEnd = new int[n+1];

        refClusterStart[0] = 0;
        refClusterEnd[0] = targToRef[targClustEnd[0] - 1] + 1;

        refClusterStart[1] = targToRef[0];
        refClusterEnd[1] = targToRef[targClustEnd[1] - 1] + 1;

        for (int j=2; j<n; ++j) {
            refClusterStart[j] = targToRef[targClustEnd[j-2]];
            refClusterEnd[j] = targToRef[targClustEnd[j] - 1] + 1;
        }

        refClusterStart[n] = targToRef[targClustEnd[n-2]];
        refClusterEnd[n] = refHapPairs.nMarkers();

        return new RefHapSegs(refHapPairs, refClusterStart, refClusterEnd);
    }

    /**
     * Return the reference haplotype pairs.
     * @return the reference haplotype pairs
     */
    public SampleHapPairs refHapPairs() {
        return refHapPairs;
    }

    /**
     * Return the target haplotype pairs.
     * @return the target haplotype pairs
     */
    public SampleHapPairs targHapPairs() {
        return targHapPairs;
    }

    /**
     * Return the reference haplotype segments.
     * @return the reference haplotype segments
     */
    public RefHapSegs refHapSegs() {
        return refHapSegs;
    }

    /**
     * Return the number of target marker clusters.
     * @return the number of target marker clusters
     */
    public int nClusters() {
        return nClusters;
    }

    /**
     * Returns the list of target samples.
     * @return the list of target samples
     */
    public Samples targetSamples() {
        return targHapPairs.samples();
    }

    /**
     * Returns the specified reference allele.
     * @param marker a marker index
     * @param haplotype a haplotype index
     * @return the specified reference allele
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nClusters()}
     * @throws IndexOutOfBoundsException if
     * {@code haplotype < 0 || haplotype >= this.refHapPairs().nHaps()}
     */
    public int refAllele(int marker, int haplotype) {
        return refAlleles[marker].get(haplotype);
    }

    /**
     * Returns the specified target allele.
     * @param marker a marker index
     * @param haplotype a haplotype index
     * @return the specified target allele
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nClusters()}
     * @throws IndexOutOfBoundsException if
     * {@code haplotype < 0 || haplotype >= targ.refHapPairs().nHaps()}
     */
    public int targetAllele(int marker, int haplotype) {
        return targAlleles[marker].get(haplotype);
    }

    /**
     * Returns the specified allele error probability.
     * @param marker the marker index
     * @return the specified allele error probability
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nClusters()}
     */
    public float errProb(int marker) {
        return errProb[marker];
    }

    /**
     * Returns {@code (1f - this.errProb(marker))}.
     * @param marker a marker index
     * @return {@code (1f - this.errProb(marker))}
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nClusters()}
     */
    public float noErrProb(int marker) {
        return 1f - errProb[marker];
    }

    /**
     * Return the probability of recombination between the specified
     * marker and the previous marker, or returns {@code 0}
     * if {@code (marker == 0)}.
     * @param marker a marker index
     * @return the specified recombination probability
     * @throws IllegalArgumentException if
     * {@code marker < 0 || marker >= this.refHapPairs().nMarkers()}
     */
    public float pRecomb(int marker) {
        return pRecomb[marker];
    }

    /**
     * Return the specified weight.
     * @param marker a marker index
     * @return the specified weight
     * @throws IllegalArgumentException if
     * {@code marker < 0 || marker >= this.refHapPairs().nMarkers()}
     */
    public double weight(int marker) {
        return weight[marker];
    }
}
