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

import blbutil.SampleFileIterator;
import beagleutil.ChromInterval;
import beagleutil.Samples;
import blbutil.Filter;
import blbutil.FilterUtils;
import haplotype.HapPair;
import haplotype.RefHapPairs;
import haplotype.SampleHapPairs;
import haplotype.WrappedHapPair;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Class {@code Data} represents a sliding marker window for
 * reference and non-reference samples.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class AllData implements Data {

    private int window = 0;
    private VcfEmission[] refData;
    private GL refEmissions;
    private VcfEmission[] targetData;  // missing markers as null entries
    private int[] refIndices;
    private int[] targetIndices;
    private GL targetEmissions;

    private final List<HapPair> refHaps;
    private final List<HapPair> targetRefHaps; // at target markers
    private final VcfWindow refWindow;
    private final VcfWindow targetWindow;

    private AllData(VcfWindow refWindow, VcfWindow targetWindow) {
        checkSampleOverlap(refWindow.samples(), targetWindow.samples());
        this.refWindow = refWindow;
        this.targetWindow = targetWindow;

        this.refData = new VcfEmission[0];
        this.refEmissions = new RefGL(refWindow.samples(), refData);
        this.targetData = new VcfEmission[0];
        this.refIndices = new int[0];
        this.targetIndices = new int[0];
        this.targetEmissions = new BasicGL(targetWindow.samples(), targetData);

        this.refHaps = new ArrayList<>(0);
        this.targetRefHaps = new ArrayList<>(0);
    }

    private static Filter<Marker> restrictToVcfMarkers(File vcfFile,
            Filter<Marker> markerFilter, ChromInterval chromInterval) {
        Set<Marker> includedMarkers = new HashSet<>(50000);
        SampleFileIterator<VcfRecord> it = new VcfIterator(vcfFile);
        if (chromInterval != null) {
            it = new IntervalVcfIterator(it, chromInterval);
        }
        if (markerFilter != null) {
            it = new FilteredVcfIterator(it, markerFilter);
        }
        while (it.hasNext()) {
            includedMarkers.add(it.next().marker());
        }
        return FilterUtils.includeFilter(includedMarkers);
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
                        + ref.id(idIndices[j-1]);
                throw new IllegalArgumentException(s);
            }
        }
    }

    /**
     * Constructs and returns a new {@code Data} instance
     * from phased genotypes for reference samples and from called genotypes
     * for non-reference samples.
     *
     * @param ref a file in VCF format with phased, non-missing genotypes.
     * The VCF records for each chromosome must be contiguous and sorted
     * in order of increasing position.
     * @param nonRef a file in VCF format with GT format field data.
     * The VCF records for each chromosome must be contiguous and sorted
     * in order of increasing position.
     * @param sampleFilter a sample filter, or {@code null} if there
     * is no sample filter.
     * @param markerFilter a marker filter, or {@code null} if there
     * is no marker filter.
     * @param chromInterval the the chromosome interval to read, or
     * {@code null} if there is no interval restriction.
     * @param usePhase {@code true} if phase information in the specified
     * VCF file should be used, and {@code false} otherwise.
     * @param impute {@code true} if markers in the reference file
     * that are missing from the non-reference file should be imputed, and
     * {@code false} otherwise.
     * @return a new {@code Data} instance.
     * @throws IllegalArgumentException if any VCF file is incorrectly
     * formatted.
     * @param pedFile a linkage-format pedigree file, or {@code null}
     * if no pedigree relationships are known.  A pedigree file must have
     * at least 4 white-space delimited columns.  The first column of the
     * pedigree file (family ID) is ignored.  The second, third, and fourth
     * columns are the individual ID, father's ID, and mother's ID respectively.
     *
     * @throws IllegalArgumentException if either VCF file contains no samples.
     * @throws IllegalArgumentException if any VCF header line
     * does not conform to the VCF specification, or if the first
     * VCF record does not conform to the VCF specification.
     * @throws IllegalArgumentException if the specified
     * {@code pedFile} is not {@code null} and contains
     * a non-blank line having less than 4 white-space delimited fields
     * or contains duplicate individual identifiers in the second column.
     * @throws NullPointerException if {@code ref==null || nonRef==null}.
     */
    public static Data gt(File ref, File nonRef, Filter<String> sampleFilter,
            Filter<Marker> markerFilter, ChromInterval chromInterval,
            File pedFile, boolean usePhase, boolean impute) {
        if (impute==false) {
            markerFilter = restrictToVcfMarkers(nonRef, markerFilter,
                    chromInterval);
        }
        SampleFileIterator<VcfRecord> filtRefIt = VcfIterator.filteredIterator(
                ref, sampleFilter, markerFilter, chromInterval);
        SampleFileIterator<VcfEmission> refIt = new VcfRefIterator(filtRefIt);
        VcfWindow refWindow = new VcfWindow(refIt);

        SampleFileIterator<VcfRecord> filtNonRefIt =
                VcfIterator.filteredIterator(nonRef, sampleFilter, markerFilter,
                chromInterval);
        SampleFileIterator<VcfEmission> targetIt = VcfEmissionIterator.gt(
                filtNonRefIt, pedFile, usePhase);
        VcfWindow targetWindow = new VcfWindow(targetIt);

        return new AllData(refWindow, targetWindow);
    }

    /**
     * Constructs and returns a new {@code Data} instance from
     * phased genotypes for reference samples and from genotype likelihoods
     * for non-reference samples.
     *
     * @param ref a file in VCF format with phased, non-missing genotypes.
     * The VCF records for each chromosome must be contiguous and sorted
     * in order of increasing position.
     * @param nonRef a file in VCF format with GL or PL format field data.
     * The VCF records for each chromosome must be contiguous and sorted
     * in order of increasing position.  If a record has both GL and PL
     * format codes, the PL format field data will be ignored.
     * @param sampleFilter a sample filter, or {@code null} if there
     * is no sample filter.
     * @param markerFilter a marker filter, or {@code null} if there
     * is no marker filter.
     * @param chromInterval the the chromosome interval to read, or
     * {@code null} if there is no interval restriction.
     * @param pedFile a linkage-format pedigree file, or {@code null}
     * if no pedigree relationships are known.  A pedigree file must have
     * at least 4 white-space delimited columns.  The first column of the
     * pedigree file (family ID) is ignored.  The second, third, and fourth
     * columns are the individual ID, father's ID, and mother's ID respectively.
     * @param maxLR maximum likelihood ratio.  If the likelihood ratio between
     * two possible genotypes is larger than {@code maxLR}, then the
     * smaller likelihood is set to 0.0, unless the change would create a
     * Mendelian inconsistency in a parent-offspring trio or duo.  In such
     * a case the unmodified likelihoods are used for all members of the
     * inconsistent duo or trio.
     * @param impute {@code true} if markers in the reference file
     * that are missing from the non-reference file should be imputed, and
     * {@code false} otherwise.
     * @return a new {@code Data} instance.
     *
     * @throws IllegalArgumentException if either VCF file contains no samples.
     * @throws IllegalArgumentException if any VCF header line
     * does not conform to the VCF specification, or if the first
     * VCF record does not conform to the VCF specification.
     * @throws IllegalArgumentException if the specified
     * {@code pedFile} is not {@code null} and contains
     * a non-blank line having less than 4 white-space delimited fields
     * or contains duplicate individual identifiers in the second column.
     * @throws IllegalArgumentException if
     * {@code Float.isNaN(maxLR) || maxLR<=1.0f}.
     * @throws NullPointerException if {@code ref==null || nonRef==null}.
     */
    public static Data gl(File ref, File nonRef, Filter<String> sampleFilter,
            Filter<Marker> markerFilter, ChromInterval chromInterval,
            File pedFile, float maxLR, boolean impute) {
        if (impute==false) {
            markerFilter = restrictToVcfMarkers(nonRef, markerFilter,
                    chromInterval);
        }
        SampleFileIterator<VcfRecord> filtRefIt = VcfIterator.filteredIterator(
                ref, sampleFilter, markerFilter, chromInterval);
        SampleFileIterator<VcfEmission> refIt = new VcfRefIterator(filtRefIt);
        VcfWindow refWindow = new VcfWindow(refIt);

        SampleFileIterator<VcfRecord> filtNonRefIt =
                VcfIterator.filteredIterator(nonRef, sampleFilter, markerFilter,
                chromInterval);
        SampleFileIterator<VcfEmission> targetIt = VcfEmissionIterator.gl(
                filtNonRefIt, pedFile, maxLR);
        VcfWindow targetWindow = new VcfWindow(targetIt);

        return new AllData(refWindow, targetWindow);
    }

    /**
     * Constructs and returns a new {@code Data} instance from
     * phased genotypes for reference samples and from genotype likelihoods
     * or called genotypes for non-reference samples.
     * @param ref a file in VCF format with phased, non-missing genotypes.
     * The VCF records for each chromosome must be contiguous and sorted
     * in order of increasing position.
     * @param nonRef a file in VCF format with GL, PL, or GT format field data.
     * The VCF records for each chromosome must be contiguous and sorted
     * in order of increasing position.  If a record has a GL or PL format
     * code, GT format field data will be ignored.  If a record has both
     * GL and PL format codes, the PL format field data will be ignored.
     * @param sampleFilter a sample filter, or {@code null} if there
     * is no sample filter.
     * @param markerFilter a marker filter, or {@code null} if there
     * is no marker filter.
     * @param chromInterval the the chromosome interval to read, or
     * {@code null} if there is no interval restriction.
     * @param pedFile a linkage-format pedigree file, or {@code null}
     * if no pedigree relationships are known.  A pedigree file must have
     * at least 4 white-space delimited columns.  The first column of the
     * pedigree file (family ID) is ignored.  The second, third, and fourth
     * columns are the individual ID, father's ID, and mother's ID respectively.
     * @param usePhase {@code true} if phase information in the specified
     * VCF record will be used when a sample's genotype emission probabilities
     * are determined by a called genotype, and {@code false} if
     * phase information in the specified VCF file record will be ignored.
     * @param maxLR maximum likelihood ratio.  If the likelihood ratio between
     * two possible genotypes is larger than {@code maxLR}, then the
     * smaller likelihood is set to 0.0, unless the change would create a
     * Mendelian inconsistency in a parent-offspring trio or duo.  In such
     * a case the unmodified likelihoods are used for all members of the
     * inconsistent duo or trio.
     * @param impute {@code true} if markers in the reference file
     * that are missing from the non-reference file should be imputed, and
     * {@code false} otherwise.
     * @return a new {@code Data} instance from
     * phased genotypes for reference samples and from genotype likelihoods
     * or called genotypes for non-reference samples.
     *
     * @throws IllegalArgumentException if either VCF file contains no samples.
     * @throws IllegalArgumentException if any VCF header line
     * does not conform to the VCF specification, or if the first
     * VCF record does not conform to the VCF specification.
     * @throws IllegalArgumentException if the specified
     * {@code pedFile} is not {@code null} and contains
     * a non-blank line having less than 4 white-space delimited fields
     * or contains duplicate individual identifiers in the second column.
     * @throws IllegalArgumentException if
     * {@code Float.isNaN(maxLR) || maxLR<=1.0f}.
     * @throws NullPointerException if {@code ref==null || nonRef==null}.
     */
    public static Data gtgl(File ref, File nonRef, Filter<String> sampleFilter,
            Filter<Marker> markerFilter, ChromInterval chromInterval,
            File pedFile, boolean usePhase, float maxLR, boolean impute) {
        if (impute==false) {
            markerFilter = restrictToVcfMarkers(nonRef, markerFilter, chromInterval);
        }
        boolean preferGL = false;
        SampleFileIterator<VcfRecord> filtRefIt = VcfIterator.filteredIterator(
                ref, sampleFilter, markerFilter, chromInterval);
        SampleFileIterator<VcfEmission> refIt = new VcfRefIterator(filtRefIt);
        VcfWindow refWindow = new VcfWindow(refIt);

        SampleFileIterator<VcfRecord> filtNonRefIt =
                VcfIterator.filteredIterator(nonRef, sampleFilter, markerFilter,
                chromInterval);
        SampleFileIterator<VcfEmission> targetIt = VcfEmissionIterator.gtgl(
                filtNonRefIt, pedFile, usePhase, maxLR, preferGL);
        VcfWindow targetWindow = new VcfWindow(targetIt);

        return new AllData(refWindow, targetWindow);
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
    public void advanceWindow(int requestedOverlap, int windowSize) {
        refData = refWindow.advanceWindow(requestedOverlap, windowSize);
        refEmissions = new RefGL(refWindow.samples(), refData);
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
        refHaps.clear();
        SampleHapPairs refHaplotypes =
                new RefHapPairs(refMarkers, refWindow.samples(), refData);
        for (int j=0, n=refHaplotypes.nSamples(); j<n; ++j) {
            refHaps.add(new WrappedHapPair(refHaplotypes, j));
        }
    }

    private void setTargetRefHaplotypes(Markers targetMarkers, VcfEmission[] refData,
            int[] refMarkerIndices) {
        assert targetMarkers.nMarkers()==refMarkerIndices.length;
        targetRefHaps.clear();
        VcfEmission[] vma = new VcfEmission[refMarkerIndices.length];
        for (int j=0; j<refMarkerIndices.length; ++j) {
                vma[j] = refData[refMarkerIndices[j]];
        }
        SampleHapPairs refHaplotypes
                = new RefHapPairs(targetMarkers, refWindow.samples(), vma);
        for (int j=0, n=refHaplotypes.nSamples(); j<n; ++j) {
            targetRefHaps.add(new WrappedHapPair(refHaplotypes, j));
        }
    }

    @Override
    public Markers markers() {
        return refEmissions.markers();
    }

    @Override
    public Markers nonRefMarkers() {
        return targetEmissions.markers();
    }

    @Override
    public int nMarkers() {
        return refEmissions.nMarkers();
    }

    @Override
    public int nNonRefMarkers() {
        return targetEmissions.markers().nMarkers();
    }

    @Override
    public int markerIndex(int nonRefIndex) {
        return refIndices[nonRefIndex];
    }

    @Override
    public int nonRefMarkerIndex(int refIndex) {
        return targetIndices[refIndex];
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
    public int nNonRefSamples() {
        return targetEmissions.nSamples();
    }

    @Override
    public Samples nonRefSamples() {
        return targetEmissions.samples();
    }

    @Override
    public int overlap() {
        return refWindow.overlap();
    }

    @Override
    public int nonRefOverlap() {
        return targetWindow.overlap();
    }

    @Override
    public int cumMarkerCnt() {
        return refWindow.cumMarkerCnt();
    }

    @Override
    public GL refEmissions() {
        return refEmissions;
    }

    @Override
    public GL nonRefEmissions() {
       return targetEmissions;
    }

    @Override
    public List<HapPair> restrictedRefHaps() {
        return new ArrayList<>(targetRefHaps);
    }

    @Override
    public List<HapPair> refHaps() {
        return new ArrayList<>(refHaps);
    }

    @Override
    public void close() {
        refWindow.close();
        targetWindow.close();
    }

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}.
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("vcf.AllData");
        return sb.toString();
    }
}
