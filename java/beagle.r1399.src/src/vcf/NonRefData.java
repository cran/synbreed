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

import beagleutil.ChromInterval;
import beagleutil.Samples;
import blbutil.Filter;
import blbutil.SampleFileIterator;
import haplotype.HapPair;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * Class {@code NonRefWindow} represents a sliding marker window for
 * non-reference samples.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class NonRefData implements Data {

    private final VcfWindow vcfWindow;

    private int window = 0;
    private Markers markers;
    private VcfEmission[] markerData;
    private GL gl;

    private NonRefData(VcfWindow vcfWindow) {
        this.vcfWindow = vcfWindow;
        this.markers = new Markers(new Marker[0]);
        this.markerData = new VcfEmission[0];
        this.gl = new BasicGL(vcfWindow.samples(), markerData);
    }

    /**
     * Constructs and returns a new {@code Data} instance
     * from called genotypes for non-reference samples.
     *
     * @param vcfFile a file in VCF format with GT format field data.
     * The VCF records for each chromosome must be contiguous and sorted
     * in order of increasing position.
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
     * VCF file should be used, and {@code false} otherwise.
     * @return a new {@code Data} instance.
     *
     * @throws IllegalArgumentException if the VCF file contains no samples.
     * @throws IllegalArgumentException if any VCF header line
     * does not conform to the VCF specification, or if the first
     * VCF record does not conform to the VCF specification.
     * @throws IllegalArgumentException if the specified
     * {@code pedFile} is not {@code null} and contains
     * a non-blank line having less than 4 white-space delimited fields
     * or contains duplicate individual identifiers in the second column.
     * @throws NullPointerException if {@code vcfFile==null}.
     */
    public static Data gt(File vcfFile, Filter<String> sampleFilter,
            Filter<Marker> markerFilter, ChromInterval chromInterval,
            File pedFile, boolean usePhase) {

        SampleFileIterator<VcfRecord> filtIt = VcfIterator.filteredIterator(
                vcfFile, sampleFilter, markerFilter, chromInterval);
        SampleFileIterator<VcfEmission> targetIt = VcfEmissionIterator.gt(filtIt,
                pedFile, usePhase);
        return new NonRefData(new VcfWindow(targetIt));
    }

    /**
     * Constructs and returns a new {@code Data} instance
     * from genotype likelihoods for non-reference samples.
     *
     * @param vcfFile a file in VCF format with GL or PL format field data.
     * The VCF records for each chromosome must be contiguous and sorted
     * in order of increasing position.  If a record has GL and PL
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
     * @return a new {@code Data} instance.
     *
     * @throws IllegalArgumentException if the VCF file contains no samples.
     * @throws IllegalArgumentException if any VCF header line
     * does not conform to the VCF specification, or if the first
     * VCF record does not conform to the VCF specification.
     * @throws IllegalArgumentException if the specified
     * {@code pedFile} is not {@code null} and contains
     * a non-blank line having less than 4 white-space delimited fields
     * or contains duplicate individual identifiers in the second column.
     * @throws IllegalArgumentException if
     * {@code Float.isNaN(maxLR) || maxLR<1.0f}.
     * @throws NullPointerException if {@code vcfFile==null}.
     */
    public static Data gl(File vcfFile, Filter<String> sampleFilter,
            Filter<Marker> markerFilter, ChromInterval chromInterval,
            File pedFile, float maxLR) {
        SampleFileIterator<VcfRecord> filtIt = VcfIterator.filteredIterator(
                vcfFile, sampleFilter, markerFilter, chromInterval);
        SampleFileIterator<VcfEmission> targetIt = VcfEmissionIterator.gl(filtIt,
                pedFile, maxLR);
        return new NonRefData(new VcfWindow(targetIt));
    }

    /**
     * Constructs and returns a new {@code Data} instance
     * from genotype likelihoods or from called genotypes for non-reference
     * samples.
     *
     * @param vcfFile a file in VCF format.  The VCF records for
     * each chromosome must be contiguous and sorted in order of
     * increasing position.  If a record has a GL or PL format
     * code, GT format field data will be ignored.  If a record has
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
     * @return a new {@code Data} instance.
     *
     * @throws IllegalArgumentException if the VCF file contains no samples.
     * @throws IllegalArgumentException if any VCF header line
     * does not conform to the VCF specification, or if the first
     * VCF record does not conform to the VCF specification.
     * @throws IllegalArgumentException if the specified
     * {@code pedFile} is not {@code null} and contains
     * a non-blank line having less than 4 white-space delimited fields
     * or contains duplicate individual identifiers in the second column.
     * @throws IllegalArgumentException if
     * {@code Float.isNaN(maxLR) || maxLR<1.0f}.
     * @throws NullPointerException if {@code vcfFile==null}.
     */
    public static Data gtgl(File vcfFile, File pedFile, boolean usePhase,
            float maxLR, Filter<String> sampleFilter, Filter<Marker> markerFilter,
            ChromInterval chromInterval) {
        boolean preferGL = false;
        SampleFileIterator<VcfRecord> filtIt = VcfIterator.filteredIterator(
                vcfFile, sampleFilter, markerFilter, chromInterval);
        SampleFileIterator<VcfEmission> targetIt = VcfEmissionIterator.gtgl(
                filtIt, pedFile, usePhase, maxLR, preferGL);
        return new NonRefData(new VcfWindow(targetIt));
    }

    @Override
    public boolean lastWindowOnChrom() {
        return vcfWindow.lastWindowOnChrom();
    }

    @Override
    public boolean canAdvanceWindow() {
        return vcfWindow.canAdvanceWindow();
    }

    @Override
    public void advanceWindow(int overlap, int windowSize) {
        markerData = vcfWindow.advanceWindow(overlap, windowSize);
        markers = extractMarkers(markerData);
        gl = new BasicGL(vcfWindow.samples(), markerData);
        ++window;
    }

    private static Markers extractMarkers(VcfEmission[] markerData) {
        Marker[] ma = new Marker[markerData.length];
        for (int j=0; j<ma.length; ++j) {
            ma[j] = markerData[j].marker();
        }
        return new Markers(ma);
    }

    @Override
    public int window() {
        return window;
    }

    @Override
    public int overlap() {
        return vcfWindow.overlap();
    }

    @Override
    public int nonRefOverlap() {
        return vcfWindow.overlap();
    }

    @Override
    public int cumMarkerCnt() {
        return vcfWindow.cumMarkerCnt();
    }

    @Override
    public Markers markers() {
        return markers;
    }

    @Override
    public Markers nonRefMarkers() {
        return markers;
    }

    @Override
    public int nMarkers() {
        return markers.nMarkers();
    }

    @Override
    public int nNonRefMarkers() {
        return markers.nMarkers();
    }

    @Override
    public int markerIndex(int nonRefIndex) {
        if (nonRefIndex < 0 || nonRefIndex >= markers.nMarkers()) {
            throw new ArrayIndexOutOfBoundsException(nonRefIndex);
        }
        return nonRefIndex;
    }

    @Override
    public int nonRefMarkerIndex(int refIndex) {
        if (refIndex < 0 || refIndex >= markers.nMarkers()) {
            throw new ArrayIndexOutOfBoundsException(refIndex);
        }
        return refIndex;
    }

    @Override
    public int nRefSamples() {
        return 0;
    }

    @Override
    public Samples refSamples() {
        return null;
    }

    @Override
    public int nNonRefSamples() {
        return vcfWindow.nSamples();
    }

    @Override
    public Samples nonRefSamples() {
       return vcfWindow.samples();
    }

    @Override
    public GL refEmissions() {
        return null;
    }

    @Override
    public GL nonRefEmissions() {
        return gl;
    }

    @Override
    public List<HapPair> restrictedRefHaps() {
        // no reference haplotypes to add
        return new ArrayList<>();
    }

    @Override
    public List<HapPair> refHaps() {
        // no reference haplotypes to add
        return new ArrayList<>();
    }

    @Override
    public void close() {
       vcfWindow.close();
    }

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}.
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("vcf.NonRefData");
        return sb.toString();
    }
}
