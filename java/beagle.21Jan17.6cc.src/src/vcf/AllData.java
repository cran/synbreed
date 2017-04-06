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
package vcf;

import beagleutil.SampleIds;
import beagleutil.Samples;
import blbutil.SampleFileIt;
import haplotype.HapPair;
import haplotype.RefHapPairs;
import haplotype.SampleHapPairs;
import haplotype.WrappedHapPair;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * <p>Class {@code AllData} represents a sliding window of
 * reference and target VCF records.
 * </p>
 * <p>Instances of class {@code AllData} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class AllData implements Data {

    private int window = 0;
    private VcfEmission[] refData;
    private SampleHapPairs refSampleHapPairs;
    private GL refEmissions;
    private VcfEmission[] targetData;  // missing markers as null entries
    private int[] refIndices;
    private int[] targetIndices;
    private GL targetEmissions;
    private final Samples allSamples;

    private final List<HapPair> refHapPairs;
    private final List<HapPair> targetRefHapPairs; // at target markers
    private final VcfWindow refWindow;
    private final RestrictedVcfWindow targetWindow;

    /**
     * Constructs and returns a new {@code AllData} instance from VCF records
     * returned by the specified {@code SampleFileIt} objects.
     *
     * @param refIt an iterator that returns reference VCF records
     * @param targetIt an iterator that returns target VCF records
     * @return a new {@code AllData} instance.
     *
     * @throws IllegalArgumentException if either the reference data or
     * target data contain no samples
     * @throws IllegalArgumentException if a format error is detected
     * in a string VCF record
     * @throws NullPointerException if {@code refIt == null || targetIt == null}
     */
    public static AllData allData(SampleFileIt<VcfEmission> refIt,
            SampleFileIt<? extends VcfEmission> targetIt) {
        if (refIt.samples().nSamples()==0 && targetIt.samples().nSamples()==0) {
            throw new IllegalArgumentException("nSamples==0");
        }
        VcfWindow refWindow = new VcfWindow(refIt);
        RestrictedVcfWindow targetWindow = new RestrictedVcfWindow(targetIt);
        return new AllData(refWindow, targetWindow);
    }

    private AllData(VcfWindow refWind, RestrictedVcfWindow targetWind) {
        checkSampleOverlap(refWind.samples(), targetWind.samples());
        this.refWindow = refWind;
        this.targetWindow = targetWind;

        this.refData = new VcfEmission[0];
        this.refSampleHapPairs = null;
        this.refEmissions = new RefGL(refWind.samples(), refData);
        this.targetData = new VcfEmission[0];
        this.refIndices = new int[0];
        this.targetIndices = new int[0];
        this.targetEmissions = new BasicGL(targetWind.samples(), targetData);
        this.allSamples = allSamples(refWind.samples(), targetWind.samples());

        this.refHapPairs = new ArrayList<>(0);
        this.targetRefHapPairs = new ArrayList<>(0);
    }

    private static Samples allSamples(Samples ref, Samples target) {
        /*
           Target samples are listed first so that sample indices agree
           with sample indices in target data genotype likelihoods.
        */
        int nRef = ref.nSamples();
        int nTarget = target.nSamples();
        int[] idIndices = new int[nRef + nTarget];
        for (int j=0; j<nTarget; ++j) {
            idIndices[j] = target.idIndex(j);
        }
        for (int j=0; j<nRef; ++j) {
            idIndices[nTarget + j] = ref.idIndex(j);
        }
        return new Samples(idIndices);
    }

    private static void checkSampleOverlap(Samples ref, Samples nonRef) {
        int nRef = ref.nSamples();
        int nNonRef = nonRef.nSamples();
        int n = nRef + nNonRef;
        int[] idIndices = new int[n];
        for (int j=0; j<nRef; ++j) {
            idIndices[j] = ref.idIndex(j);
        }
        for (int j=0; j<nNonRef; ++j) {
            idIndices[nRef + j] = nonRef.idIndex(j);
        }
        Arrays.sort(idIndices);
        for (int j=1; j<idIndices.length; ++j) {
            if (idIndices[j-1]==idIndices[j]) {
                String s = "Overlap between reference and non-reference samples: "
                        + SampleIds.instance().id(idIndices[j-1]);
                throw new IllegalArgumentException(s);
            }
        }
    }

    @Override
    public boolean lastWindowOnChrom() {
        return refWindow.lastWindowOnChrom();
    }

    @Override
    public boolean canAdvanceWindow() {
        return refWindow.canAdvanceWindow();
    }

    @Override
    public void advanceWindow(int overlap, int windowSize) {
        Samples refSamples = refWindow.samples();
        refData = refWindow.advanceWindow(overlap, windowSize);
        refEmissions = new RefGL(refSamples, refData);
        refSampleHapPairs = new RefHapPairs(refEmissions.markers(), refSamples, refData);
        targetData = targetWindow.advanceWindow(refEmissions.markers());
        refIndices = refIndices(targetData);
        targetIndices = targetIndices(targetData);
        targetEmissions = targetEmissions(targetWindow.samples(),
                targetData, refIndices);
        ++window;
        setRefHaplotypes(refEmissions.markers(), refData);
        setTargetRefHaplotypes(targetEmissions.markers(), refData, refIndices);
    }

    @Override
    public int window() {
        return window;
    }

    private static int[] refIndices(VcfEmission[] vma) {
        int nonNullCnt = 0;
        for (VcfEmission vm : vma) {
            if (vm!=null) {
                ++nonNullCnt;
            }
        }
        int[] inclusionMap = new int[nonNullCnt];
        int index = 0;
        for (int j=0; j<vma.length; ++j) {
            if (vma[j]!=null) {
                inclusionMap[index++] = j;
            }
        }
        if (index != inclusionMap.length) {
            throw new IllegalStateException("vma modification detected");
        }
        return inclusionMap;
    }

    private static int[] targetIndices(VcfEmission[] vma) {
        int[] inclusionMap = new int[vma.length];
        int index = 0;
        for (int j=0; j<inclusionMap.length; ++j) {
            if (vma[j]!=null) {
                inclusionMap[j] = index++;
            }
            else {
                inclusionMap[j] = -1;
            }
        }
        return inclusionMap;
    }

    private static GL targetEmissions(Samples samples,
            VcfEmission[] vma, int[] refMarkerIndex) {
        VcfEmission[] restricted = new VcfEmission[refMarkerIndex.length];
        for (int j=0; j<refMarkerIndex.length; ++j) {
            restricted[j] = vma[refMarkerIndex[j]];
        }
        return new BasicGL(samples, restricted);
    }

    private void setRefHaplotypes(Markers refMarkers, VcfEmission[] refData) {
        refHapPairs.clear();
        SampleHapPairs refHaplotypes =
                new RefHapPairs(refMarkers, refWindow.samples(), refData);
        for (int j=0, n=refHaplotypes.nSamples(); j<n; ++j) {
            refHapPairs.add(new WrappedHapPair(refHaplotypes, j));
        }
    }

    private void setTargetRefHaplotypes(Markers targetMarkers, VcfEmission[] refData,
            int[] refMarkerIndices) {
        assert targetMarkers.nMarkers()==refMarkerIndices.length;
        targetRefHapPairs.clear();
        VcfEmission[] vma = new VcfEmission[refMarkerIndices.length];
        for (int j=0; j<refMarkerIndices.length; ++j) {
                vma[j] = refData[refMarkerIndices[j]];
        }
        SampleHapPairs refHaplotypes
                = new RefHapPairs(targetMarkers, refWindow.samples(), vma);
        for (int j=0, n=refHaplotypes.nSamples(); j<n; ++j) {
            targetRefHapPairs.add(new WrappedHapPair(refHaplotypes, j));
        }
    }

    @Override
    public int targetOverlap() {
        return targetWindow.overlap();
    }

    @Override
    public int overlap() {
        return refWindow.overlap();
    }

    @Override
    public int nTargetMarkers() {
        return targetEmissions.markers().nMarkers();
    }

    @Override
    public int nTargetMarkersSoFar() {
        return targetWindow.cumMarkerCnt();
    }

    @Override
    public Markers targetMarkers() {
        return targetEmissions.markers();
    }


    @Override
    public int nMarkers() {
        return refEmissions.nMarkers();
    }

    @Override
    public int nMarkersSoFar() {
        return refWindow.cumMarkerCnt();
    }

    @Override
    public Markers markers() {
        return refEmissions.markers();
    }

    @Override
    public int targetMarkerIndex(int refIndex) {
        return targetIndices[refIndex];
    }

    @Override
    public int markerIndex(int nonRefIndex) {
        return refIndices[nonRefIndex];
    }

    @Override
    public int nTargetSamples() {
        return targetEmissions.nSamples();
    }

    @Override
    public Samples targetSamples() {
        return targetEmissions.samples();
    }

    @Override
    public int nRefSamples() {
        return refWindow.nSamples();
    }

    @Override
    public Samples refSamples() {
        return refWindow.samples();
    }

    @Override
    public int nAllSamples() {
        return allSamples.nSamples();
    }

    @Override
    public Samples allSamples() {
        return allSamples;
    }


    @Override
    public GL targetGL() {
       return targetEmissions;
    }

    @Override
    public List<HapPair> restrictedRefHapPairs() {
        return new ArrayList<>(targetRefHapPairs);
    }

    @Override
    public List<HapPair> refHapPairs() {
        return new ArrayList<>(refHapPairs);
    }

    @Override
    public SampleHapPairs refSampleHapPairs() {
        return refSampleHapPairs;
    }

    @Override
    public void close() {
        refWindow.close();
        targetWindow.close();
    }

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(20);
        sb.append(this.getClass().toString());
        return sb.toString();
    }
}
