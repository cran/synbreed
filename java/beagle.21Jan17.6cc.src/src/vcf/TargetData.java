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

import beagleutil.Samples;
import blbutil.SampleFileIt;
import haplotype.HapPair;
import haplotype.SampleHapPairs;
import java.util.ArrayList;
import java.util.List;

/**
 * <p>Class {@code TargetData} represents a sliding window of
 * target VCF records.
 * </p>
 * <p>Instances of class {@code TargetData} are not thread-safe.
 * </p>
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class TargetData implements Data {

    private final VcfWindow vcfWindow;

    private int window = 0;
    private Markers markers;
    private VcfEmission[] markerData;
    private GL gl;

    /**
     * Constructs and returns a new {@code TargetData} instance from
     * VcfRecords returned by the specified {@code SampleFileIt} objects.
     *
     * @param it an iterator that returns the target VCF records
     * @return a new {@code TargetData} instance
     *
     * @throws IllegalArgumentException if the data returned by
     * the specified iterator contains no samples
     * @throws IllegalArgumentException if a format error is detected
     * in a string VCF record
     * @throws NullPointerException if {@code it == null}
     */
    public static TargetData targetData(SampleFileIt<? extends VcfEmission> it) {
        if (it.samples().nSamples()==0) {
            throw new IllegalArgumentException("nSamples==0");
        }
        return new TargetData(new VcfWindow(it));
    }

    private TargetData(VcfWindow vcfWindow) {
        this.vcfWindow = vcfWindow;
        this.markers = Markers.create(new Marker[0]);
        this.markerData = new VcfEmission[0];
        this.gl = new BasicGL(vcfWindow.samples(), markerData);
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
        return Markers.create(ma);
    }

    @Override
    public int window() {
        return window;
    }


    @Override
    public int targetOverlap() {
        return vcfWindow.overlap();
    }

    @Override
    public int overlap() {
        return vcfWindow.overlap();
    }

    @Override
    public int nTargetMarkers() {
        return markers.nMarkers();
    }

    @Override
    public int nTargetMarkersSoFar() {
        return vcfWindow.cumMarkerCnt();
    }

    @Override
    public Markers targetMarkers() {
        return markers;
    }

    @Override
    public int nMarkers() {
        return markers.nMarkers();
    }

    @Override
    public int nMarkersSoFar() {
        return vcfWindow.cumMarkerCnt();
    }

    @Override
    public Markers markers() {
        return markers;
    }

    @Override
    public int targetMarkerIndex(int refIndex) {
        if (refIndex < 0 || refIndex >= markers.nMarkers()) {
            throw new ArrayIndexOutOfBoundsException(refIndex);
        }
        return refIndex;
    }

    @Override
    public int markerIndex(int nonRefIndex) {
        if (nonRefIndex < 0 || nonRefIndex >= markers.nMarkers()) {
            throw new ArrayIndexOutOfBoundsException(nonRefIndex);
        }
        return nonRefIndex;
    }

    @Override
    public int nTargetSamples() {
        return vcfWindow.nSamples();
    }

    @Override
    public Samples targetSamples() {
       return vcfWindow.samples();
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
    public int nAllSamples() {
        return nTargetSamples();
    }

    @Override
    public Samples allSamples() {
        return targetSamples();
    }

    @Override
    public GL targetGL() {
        return gl;
    }

    @Override
    public List<HapPair> restrictedRefHapPairs() {
        // no reference haplotypes to add
        return new ArrayList<>();
    }

    @Override
    public List<HapPair> refHapPairs() {
        // no reference haplotypes to return
        return new ArrayList<>();
    }


    @Override
    public SampleHapPairs refSampleHapPairs() {
        // no reference haplotypes to return
        return null;
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
