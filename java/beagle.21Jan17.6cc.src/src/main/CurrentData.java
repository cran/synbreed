/*
 * Copyright (C) 2015 browning
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

import beagleutil.Samples;
import haplotype.BasicSampleHapPairs;
import haplotype.HapPair;
import haplotype.SampleHapPairs;
import haplotype.Weights;
import java.util.List;
import vcf.Data;
import vcf.GL;
import vcf.Markers;
import vcf.SplicedGL;

/**
 * <p>Class {@code CurrentData} represents input data for the current marker
 * window.  All marker indices returned my methods of class {@code CurrentData}
 * are indexed with respect to the current marker window.
 * </p>
 * <p>Instances of class {@code CurrentData} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class CurrentData {

    private static final float MIN_GEN_DIST = 1e-7f;

    private final int window;
    private final SampleHapPairs initHaps;
    private final int prevSpliceStart;
    private final int nextOverlapStart;
    private final int nextSpliceStart;
    private final int nextTargetSpliceStart;
    private final int nextTargetOverlapStart;

    private final GL targetGL;
    private final NuclearFamilies families;
    private final Weights weights;

    private final Samples refSamples;
    private final Samples targetSamples;
    private final Samples allSamples;

    private final Markers markers;
    private final Markers targetMarkers;
    private final int[] targetMarkerIndex;
    private final int[] markerIndex;

    private final List<HapPair> restRefHapPairs;
    private final SampleHapPairs refSampleHapPairs;
    private final SampleHapPairs restrictedRefSampleHapPairs;

    private final float[] recombRate;

    /**
     * Constructs a new {@code CurrentData} instance from the specified
     * data.
     *
     * @param par the analysis parameters
     * @param genMap the genetic map or {@code null} if no
     * genetic map is specified
     * @param data input data for the current marker window
     * @param overlapHaps haplotype constraints in the overlap with previous
     * window or {@code null} if no such constraints exist
     * @param families the parent-offspring relationships
     *
     * @throws IllegalArgumentException if
     * {@code data.targetSamples().equals(families.samples()) == false}
     * @throws IllegalArgumentException if
     * {@code (overlapHaps != null
     * && data.targetSamples().equals(overlapHaps.samples()) == false)}
     * @throws IllegalArgumentException if
     * {@code (overlapHaps != null &&
     * overlapHaps.marker(j).equals(data.targetGL().marker(j) == false)}
     * for some {@code j} satisfying
     * {@code (0 <= j && j <= overlapHaps.nMarkers())}
     * @throws NullPointerException if any parameter is {@code null}
     */
    public CurrentData(Par par, GeneticMap genMap, Data data,
            SampleHapPairs overlapHaps, NuclearFamilies families) {
        if (families.samples().equals(data.targetSamples())==false) {
            throw new IllegalArgumentException("inconsistent samples");
        }
        if (overlapHaps != null
                && data.targetSamples().equals(overlapHaps.samples())==false) {
            throw new IllegalArgumentException("inconsistent samples");
        }
        this.window = data.window();
        this.initHaps = overlapHaps;
        this.prevSpliceStart = data.overlap()/2;
        this.nextOverlapStart = CurrentData.this.nextOverlapStart(data, par.overlap());
        this.nextSpliceStart = (data.nMarkers() + nextOverlapStart)/2;
        this.nextTargetOverlapStart = targetIndex(data, nextOverlapStart);
        this.nextTargetSpliceStart = targetIndex(data, nextSpliceStart);

        this.families = families;
        this.weights = new Weights(families);
        this.targetGL = (overlapHaps==null) ? data.targetGL() :
                new SplicedGL(overlapHaps, data.targetGL());

        this.refSamples = data.refSamples();
        this.targetSamples = data.targetSamples();
        this.allSamples = data.allSamples();
        this.markers = data.markers();
        this.targetMarkers = data.targetMarkers();
        this.targetMarkerIndex = refToTargetMarker(data);
        this.markerIndex = targetToRefMarker(data);

        this.restRefHapPairs = data.restrictedRefHapPairs();
        this.refSampleHapPairs = data.refSampleHapPairs();
        this.restrictedRefSampleHapPairs = refSamples != null ?
                new BasicSampleHapPairs(refSamples, restRefHapPairs) : null;
        this.recombRate = recombRate(targetMarkers, genMap, par.mapscale());
    }

    /* Returns the index of the first marker in the overlap */
    private static int nextOverlapStart(Data data, int targetOverlap) {
        if (targetOverlap < 0) {
            throw new IllegalArgumentException(String.valueOf(targetOverlap));
        }
        if (targetOverlap==0 || data.lastWindowOnChrom()) {
            return data.nMarkers();
        }
        Markers markers = data.markers();
        int nextOverlap = Math.max(0, data.nMarkers() - targetOverlap);
        while (nextOverlap>0
                && markers.marker(nextOverlap).pos()
                    == markers.marker(nextOverlap - 1).pos()) {
            --nextOverlap;
        }
        return nextOverlap;
    }

    /* Returns the index of the first marker after the next splice point */
    private static int nextSpliceStart(Data data, int overlap) {
        if (data.canAdvanceWindow() && data.lastWindowOnChrom()==false) {
            return data.nMarkers() - overlap + (overlap/2);
        }
        else {
            return data.nMarkers();
        }
    }

    /* first target index on or after specified ref index */
    private static int targetIndex(Data data, int refIndex) {
        int i=0;
        while (i<data.nTargetMarkers() && data.markerIndex(i)<refIndex) {
            ++i;
        }
        return i;
    }

    private static float[] recombRate(Markers markers, GeneticMap map,
            float mapScale) {
        if (map==null) {
            return null;
        }
        else {
            double c = -2.0*mapScale;
            float[] rr = new float[markers.nMarkers()];
            rr[0] = 0.0f;
            double lastGenPos = map.genPos(markers.marker(0));
            for (int j=1; j<rr.length; ++j) {
                double genPos = map.genPos(markers.marker(j));
                double genDist = Math.max(Math.abs(genPos - lastGenPos), MIN_GEN_DIST);
                rr[j] = (float) -Math.expm1(c*genDist);
                lastGenPos = genPos;
            }
            return rr;
        }
    }

    /**
     * Returns the marker window index.
     * @return the marker window index
     */
    public int window() {
        return window;
    }

    /**
     * Returns the first marker index in the overlap between this
     * marker window and the next marker window, or
     * returns {@code this.nMarkers()} there is no overlap.
     * @return the first marker index in the overlap between this
     * marker window and the next marker window
     */
    public int nextOverlapStart() {
        return nextOverlapStart;
    }

    /**
     * Returns the first target marker index in the overlap between this
     * marker window and the next marker window, or
     * returns {@code this.nMarkers()} if there is no overlap or if there are
     * no target markers in the overlap.
     * @return the first target marker index in the overlap between this
     * marker window and the next marker window
     */
    public int nextTargetOverlapStart() {
        return nextTargetOverlapStart;
    }

    /**
     * Returns the first marker index after the splice point with
     * the previous marker window. Returns 0 if the current marker window
     * is the first marker window.
     * @return the first marker index after the splice point with
     * the previous marker window
     */
    public int prevSpliceStart() {
        return prevSpliceStart;
    }

    /**
     * Returns the first marker index after the splice point between this
     * marker window and the next marker window, or returns
     * {@code this.nMarkers()} if there is no overlap or if there are
     * no markers after the splice point.
     * @return the first marker index after the next splice point
     */
    public int nextSpliceStart() {
        return nextSpliceStart;
    }

    /**
     * Returns the first target marker index after the splice point with
     * the previous marker window. Returns 0 if the current marker window
     * is the first marker window.
     * @return the first target marker index after the splice point with
     * the previous marker window
     */
    public int prevTargetSpliceStart() {
        return initHaps==null ? 0 : initHaps.nMarkers();
    }

    /**
     * Returns the first target marker index after the splice point between this
     * marker window and the next marker window, or returns
     * {@code this.nTargetMarkers()} if there is no overlap or if there are
     * no target markers after the splice point
     * @return the first target marker index after the next splice point
     */
    public int nextTargetSpliceStart() {
        return nextTargetSpliceStart;
    }

    /**
     * Returns the target data haplotype pairs in the segment of the current
     * marker window preceding the splice point with the previous marker window:
     * {@code this.targetMarkers().restrict(0, this.prevTargetSplice())}
     * @return the target data haplotype pairs in the segment of the current
     * marker window preceding the splice point with the previous marker window
     */
    public SampleHapPairs initHaps() {
        return initHaps;
    }

    private int[] refToTargetMarker(Data data) {
        int[] ia = new int[data.nMarkers()];
        for (int j=0; j<ia.length; ++j) {
            ia[j] = data.targetMarkerIndex(j);
        }
        return ia;
    }

    private static int[] targetToRefMarker(Data data) {
        int[] ia = new int[data.nTargetMarkers()];
        for (int j=0; j<ia.length; ++j) {
            ia[j] = data.markerIndex(j);
        }
        return ia;
    }

    /**
     * Returns the parent-offspring relationships.
     * @return the parent-offspring relationships
     */
    public NuclearFamilies families() {
        return families;
    }

    /**
     * Returns the per-haplotype weights.
     * @return the per-haplotype weights
     */
    public Weights weights() {
        return weights;
    }

    /**
     * Returns the number of reference samples.
     * @return the number of reference samples
     */
    public int nRefSamples() {
        return refSamples == null ? 0 : refSamples.nSamples();
    }

    /**
     * Returns the list of reference samples, or {@code null} if
     * there are no reference samples.
     * @return the list of reference samples, or {@code null} if
     * there are no reference samples
     */
    public Samples refSamples() {
        return refSamples;
    }

    /**
     * Returns the number of target samples.
     * @return the number of target samples
     */
    public int nTargetSamples() {
        return targetSamples.nSamples();
    }

    /**
     * Returns the list of target samples.
     * @return the list of target samples
     */
    public Samples targetSamples() {
        return targetSamples;
    }


    /**
     * Returns the number of reference and target samples.
     * @return the number of reference and target samples
     */
    public int nAllSamples() {
        return allSamples.nSamples();
    }

     /**
      * Returns a list of all target and reference samples.
      * Target samples are listed first in the same order as the list returned
      * by {@code this.targetSamples()}. Reference samples are listed last
      * in the same order as the list returned by {@code this.refSamples()}.
      * @return a list of all target and reference samples
      */
    public Samples allSamples() {
        return allSamples;
    }

    /**
     * Returns the number of target data markers.
     * @return the number of target data markers
     */
    public int nTargetMarkers() {
        return targetMarkers.nMarkers();
    }

    /**
     * Returns the list of target data markers.
     * @return the list of target data markers
     */
    public Markers targetMarkers() {
        return targetMarkers;
    }

    /**
     * Returns the number of reference data markers.
     * @return the number of reference data markers
     */
     public int nMarkers() {
         return markers.nMarkers();
     }

    /**
     * Returns the list of reference data markers.
     * @return the list of reference data markers
     */
     public Markers markers() {
         return markers;
     }

    /**
     * Returns the index of the specified marker in the reference data markers.
     * @param targetMarker index of a marker in the list of target data markers
     * @return the index of the specified marker in the reference data markers
     * @throws IndexOutOfBoundsException if
     * {@code targetMarker < 0 || targetMarker >= this.nTargetMarkers()}
     */
    public int markerIndex(int targetMarker) {
        return markerIndex[targetMarker];
    }

    /**
     * Returns an array of length {@code this.nTargetMarkers()} which maps
     * the {@code k}-th marker in the list of target data markers to the
     * index of the marker in the list of reference data markers.
     * @return an array of length {@code this.nTargetMarkers()} which maps
     * the {@code k}-th marker in the list of target data markers to the
     * index of the marker in the list of reference data markers
     */
    public int[] markerIndices() {
        return markerIndex.clone();
    }

    /**
     * Returns the index of the specified marker in the target data, or
     * returns -1 if the marker is not present in the target data.
     * @param marker index of a marker in the reference data
     * @return the index of the specified marker in the target data, or
     * returns -1 if the marker is not present in the target data
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}.
     */
     public int targetMarkerIndex(int marker) {
         return targetMarkerIndex[marker];
     }

    /**
     * Returns an array of length {@code this.nMarkers()} whose {@code k}-th
     * element is the index of the {@code k}-th marker in the list of target
     * markers or is -1 if the marker is not present in the target data.
     * @return an array of length {@code this.nMarkers()} whose {@code k}-th
     * element is the index of the {@code k}-th marker in the list of target
     * markers or is -1 if the marker is not present in the target data
     */
     public int[] targetMarkerIndices() {
         return targetMarkerIndex.clone();
     }

    /**
     * Add the reference haplotype pairs that are restricted
     * to the target data markers to the specified list.
     * @param list a list of haplotype pairs for target data markers
     * @throws NullPointerException if {@code list == null}
     */
    public void addRestrictedRefHapPairs(List<HapPair> list) {
        list.addAll(restRefHapPairs);
    }

    /**
     * Returns a list of reference haplotype pairs that are restricted
     * to the target data markers, or returns {@code null}
     * if there are no reference samples.
     * @return a list of reference haplotype pairs that are restricted
     * to the target data markers
     */
    public SampleHapPairs restrictedRefSampleHapPairs() {
        return restrictedRefSampleHapPairs;
    }

    /**
     * Returns a list of reference haplotype pairs, or returns {@code null}
     * if there are no reference samples.
     * @return a list of reference haplotype pairs
     */
    public SampleHapPairs refSampleHapPairs() {
        return refSampleHapPairs;
    }

    /**
     * Returns the genotype likelihoods for the
     * target samples at the target data markers.
     * @return the genotype likelihoods for the
     * target samples at the target data markers.
     */
    public GL targetGL() {
        return targetGL;
    }

    /**
     * Returns an array whose initial element is {@code 0} and whose
     * {@code j}-th element for {@code j > 0} is the recombination rate
     * between the target markers with indices {@code (j - 1)} and {@code j}.
     *
     * @return inter-marker recombination rates for the target markers
     */
    public float[] recombRate() {
        return recombRate==null ? null : recombRate.clone();
    }
}
