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
import beagleutil.Samples;
import blbutil.Const;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import main.Logger;

/**
 * Class {@code VcfWindow} represents a sliding window of VCF records.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class VcfWindow {

    private final SampleFileIterator<VcfEmission> it;
    private int overlap;
    private int cumMarkerCnt;
    private VcfEmission[] data;

    private VcfEmission next;

    /**
     * Constructs a new {@code VcfWindow} instance.
     * @param it an iterator that returns {@code VcfEmission} objects.

     * @throws IllegalArgumentException if {@code it.hasNext()==false}.
     * @throws IllegalArgumentException if a format error in the input data
     * is encountered.
     * @throws NullPointerException if {@code it==null}.
     */
    public VcfWindow(SampleFileIterator<VcfEmission> it) {
        if (it.hasNext()==false) {
            StringBuilder sb = new StringBuilder(100);
            sb.append("No VCF records found.");
            if (it.file()!=null) {
                sb.append(" (file: ");
                sb.append(it.file());
                sb.append(")");
            }
            sb.append(Const.nl);
            sb.append("Check that the chromosome identifiers are the same in each input VCF");
            sb.append(Const.nl);
            sb.append("file and in the \'chrom=\' command line argument (if \'chrom=\' is used).");
            throw new IllegalArgumentException(sb.toString());
        }
        this.overlap = 0;
        this.cumMarkerCnt = 0;
        this.it = it;
        this.data = new VcfEmission[0];
        this.next = it.next();
    }

    /**
     * Returns {@code true} if the sliding marker window is the last
     * marker window for the chromosome and returns {@code false}
     * otherwise.
     * @return {@code true} if the sliding marker window is the last
     * marker window for the chromosome and {@code false} otherwise.
     */
    public boolean lastWindowOnChrom() {
        return next==null || data.length==0
                || sameChrom(next, data[data.length-1])==false;
    }

    /**
     * Returns {@code true} if the sliding marker window can advance
     * and returns {@code false} otherwise.
     * @return {@code true} if the sliding marker window can advance
     * and returns {@code false} otherwise.
     */
    public boolean canAdvanceWindow() {
        return next!=null;
    }

    /**
     * Advances the sliding marker window, and returns the advanced window.
     *
     * @param overlap the requested number of markers of overlap between the
     * marker window immediately before the method invocation and the
     * marker window immediately after the method returns.
     * If the marker window immediately before method invocation contains
     * the final marker on a chromosome, there will be 0 markers
     * in the overlap.
     * @param windowSize the requested number of the markers in the
     * advanced marker window.
     * @return the advanced window.
     *
     * @throws IllegalArgumentException if a format error in the input data
     * is encountered.
     * @throws IllegalArgumentException if
     * {@code overlap<0 || overlap>=windowSize}
     * @throws IllegalStateException if
     * {@code this.canAdvanceWindow()==false}
     */
    public VcfEmission[] advanceWindow(int overlap, int windowSize) {
        checkParameters(overlap, windowSize);
        if (canAdvanceWindow()==false) {
            throw new IllegalStateException("canAdvanceWindow()==false");
        }
        if (overlap > data.length) {
            overlap = data.length;
        }
        if (data.length==0 || sameChrom(next, data[data.length-1])==false) {
            overlap = 0;
        }
        assert next!=null;
        VcfEmission first = next;
        List<VcfEmission> newWindow = new ArrayList<>(windowSize);
        for (int j=(data.length-overlap); j<data.length; ++j) {
            newWindow.add(data[j]);
        }
        for (int j=overlap;
                (j<windowSize && next!=null && sameChrom(next, first)); ++j) {
            newWindow.add(next);
            next = readRecord(it);
        }
        // add all markers at the same marker position
        VcfEmission last = newWindow.get(newWindow.size()-1);
        while (next!=null && sameChromAndPosition(last, next)) {
            newWindow.add(next);
            next = readRecord(it);
        }
        this.overlap = overlap;
        this.data = newWindow.toArray(new VcfEmission[0]);
        this.cumMarkerCnt += (data.length - overlap);
        return data.clone();
    }

    private void checkParameters(int overlap, int windowSize) {
        if (overlap < 0 || overlap >= windowSize) {
            String s = "overlap=" + overlap + "windowSize=" + windowSize;
            throw new IllegalArgumentException(s);
        }
    }

    /**
     * Advances the sliding marker window, and returns the advanced window.
     * The returned array will have length {@code markers.nMarkers()}.
     * Markers not found in the data source will have {@code null}
     * entries in the returned array.
     *
     * @param markers a set of markers in the advanced window.
     * @return the advanced marker window.
     *
     * @throws IllegalArgumentException if {@code markers.nMarkers()==0}
     * @throws IllegalArgumentException if any two of the specified markers
     * are on different chromosomes.
     * @throws IllegalArgumentException if the specified markers do not
     * cause the marker window to advance.
     * @throws IllegalArgumentException if a format error in the input data
     * is encountered.
     * @throws IllegalArgumentException if the input data does not contain
     * any of the specified markers.
     * @throws IllegalStateException if {@code canAdvanceWindow()==false}.
     */
    public VcfEmission[] advanceWindow(Markers markers) {
        if (canAdvanceWindow()==false) {
            throw new IllegalStateException("canAdvanceWindow()==false");
        }
        checkMarkers(markers);
        VcfEmission[] prevWindow = data;
        VcfEmission[] nextWindow = new VcfEmission[markers.nMarkers()];
        int firstNewMarkerIndex = firstNewMarkerIndex(markers, prevWindow);
        this.overlap = copyOverlap(markers, firstNewMarkerIndex, prevWindow,
                nextWindow);
        boolean addedNew = readDataForNewmarkers(markers, nextWindow,
                firstNewMarkerIndex);
        if (overlap==0 && addedNew==false) {
            String s = missingTargetMarkersErr(markers);
            throw new IllegalArgumentException(s);
        }
        this.data = nextWindow;
        return nextWindow.clone();
    }

    private void checkMarkers(Markers markers) {
        if (markers.nMarkers()==0) {
            throw new IllegalArgumentException("markers.nMarkers()=0");
        }
        Marker first = markers.marker(0);
        Marker last = markers.marker(markers.nMarkers()-1);
        if (first.chromIndex() != last.chromIndex()) {
            String s = "inconsistent chromosomes:" + Const.nl
                    + first + Const.nl + last;
            throw new IllegalArgumentException(s);
        }
    }

    private static int firstNewMarkerIndex(Markers markers,
            VcfEmission[] prevData) {
        assert markers.nMarkers()>0;
        int lastChrom = Integer.MIN_VALUE;
        int lastPos = Integer.MIN_VALUE;
        int lastNonNullIndex = lastNonNullIndex(prevData);
        if (lastNonNullIndex >= 0) {
            lastChrom = prevData[lastNonNullIndex].marker().chromIndex();
            lastPos = prevData[lastNonNullIndex].marker().pos();
        }
        if (markers.marker(0).chromIndex() != lastChrom) {
            return 0;
        }
        else {
            int n = markers.nMarkers();
            int index=0;
            Marker m = markers.marker(index);
            while (m!=null && m.chromIndex()==lastChrom && m.pos()<=lastPos ) {
                m = ++index<n ? markers.marker(index) : null;
            }
            if (index==markers.nMarkers()) {
                String s = "markers do not advance current marker window";
                throw new IllegalArgumentException(s);
            }
            return index;
        }
    }

    private static int lastNonNullIndex(VcfEmission[] data) {
        int index = data.length - 1;
        while (index>=0 && data[index]==null) {
            --index;
        }
        return index;
    }

    private static String missingTargetMarkersErr(Markers markers) {
        return "no target markers found in interval: "
                + interval(markers)
                + Const.nl
                + "Check that target VCF file markers have identical CHROM, POS, REF,"
                + Const.nl
                + "and ALT fields as the corresponding reference VCF file markers.";
    }

    private static String interval(Markers markers) {
        Marker a = markers.marker(0);
        Marker b = markers.marker(markers.nMarkers()-1);
        String start = a.chrom() + Const.colon + a.pos();
        String end = String.valueOf(b.pos());
        if (a.chromIndex()!=b.chromIndex()) {
            end = b.chromIndex() + Const.colon + end;

        }
        return (start + Const.hyphen + end);
    }

    private static int copyOverlap(Markers markers, int firstNewMarkerIndex,
            VcfEmission[] prev, VcfEmission[] next) {
        int overlap = 0;
        Map<Marker, VcfEmission> map = new HashMap<>();
        for (VcfEmission vm : prev) {
            if (vm!=null) {
                map.put(vm.marker(), vm);
            }
        }
        for (int j=0; j<firstNewMarkerIndex; ++j) {
            next[j] = map.get(markers.marker(j));
            if (next[j]!=null) {
                ++overlap;
            }
        }
        return overlap;
    }
    /*
     * Returns true if one or more new markers were added, and returns false
     * otherwise.
     */
    private boolean readDataForNewmarkers(Markers markers,
            VcfEmission[] nextWindow, int firstNewIndex) {
        int addedCnt = 0;
        if (firstNewIndex==0 && markers.nMarkers()>0) {
            int chromIndex = markers.marker(0).chromIndex();
            while (next!=null && next.marker().chromIndex()!=chromIndex) {
                next = readRecord(it);
            }
        }
        for (int j=firstNewIndex; j<nextWindow.length; ++j) {
            boolean success = readDataForMarker(markers, j);
            if (success) {
                assert next.marker().equals(markers.marker(j));
                ++addedCnt;
                nextWindow[j] = next;
                next = readRecord(it);
            }
        }
        cumMarkerCnt += addedCnt;
        return (addedCnt>0);
    }

    /* returns true if data for marker was read */
    private boolean readDataForMarker(Markers markers, int index) {
        Marker m = markers.marker(index);
        while (next!=null
                && next.marker().chromIndex()==m.chromIndex()
                && next.marker().pos()<m.pos() ) {
            next = readRecord(it);
        }
        while (next!=null
                && next.marker().chromIndex()==m.chromIndex()
                && next.marker().pos()==m.pos()
                && markers.contains(next.marker())==false) {
            next = readRecord(it);
        }
        return next!=null && m.equals(next.marker());
    }

    private boolean sameChrom(VcfEmission a, VcfEmission b) {
        return a.marker().chromIndex()==b.marker().chromIndex();
    }

    private boolean sameChromAndPosition(VcfEmission a, VcfEmission b) {
        return a.marker().chromIndex()==b.marker().chromIndex()
                && a.marker().pos()==b.marker().pos();
    }

    private VcfEmission readRecord(SampleFileIterator<VcfEmission> it) {
        VcfEmission e = it.hasNext() ? it.next() : null;
        if (e!=null && e.isMissingData()) {
            String s = "Skipping VCF record with no data: " + e.marker();
            Logger.getInstance().println(s);
            e = it.hasNext() ? it.next() : null;
        }
        return e;
    }

    /**
     * Returns the file from which data are read, or returns
     * {@code null} if the file source is unknown or if the
     * source is standard input.
     * @return the file from which data are read, or returns
     * {@code null} if the file source is unknown or if the
     * source is standard input.
     */
    public File file() {
        return it.file();
    }

    /**
     * Returns the list of samples.
     * @return the list of samples.
     */
    public Samples samples() {
        return it.samples();
    }

    /**
     * Returns the number of samples.
     * @return the number of samples.
     */
    public int nSamples() {
        return it.samples().nSamples();
    }

    /**
     * Returns the number of markers overlap between the current marker window
     * and the previous marker window.  Returns 0 if the current marker window
     * is the first marker window.
     *
     * @return the number of markers overlap between the current marker window
     * and the previous marker window.
     */
    public int overlap() {
        return overlap;
    }

    /**
     * Returns the number of markers in the union of the current marker window
     * and all previous marker windows.
     *
     * @return the number of markers in the union of the current marker window
     * and all previous marker windows.
     */
    public int cumMarkerCnt() {
        return cumMarkerCnt;
    }

    /**
     * Closes any I/O resources controlled by this
     * {@code VcfWindow} object.
     */
    public void close() {
        it.close();
    }

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}.
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("vcf.VcfWindow");
        return sb.toString();
    }
}
