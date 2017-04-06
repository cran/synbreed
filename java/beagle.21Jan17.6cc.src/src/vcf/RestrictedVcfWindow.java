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
import blbutil.Const;
import blbutil.SampleFileIt;
import java.io.Closeable;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * <p>Class {@code RestrictedVcfWindow} represents a sliding window of VCF
 * records.
 * </p>
 * <p>Instances of class {@code RestrictedVcfWindow} are not thread.safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class RestrictedVcfWindow  implements Closeable {

    private final SampleFileIt<? extends VcfEmission> it;
    private final List<VcfEmission> window;
    private Markers markers;
    private int overlap = 0;
    private int cumMarkerCnt;
    private VcfEmission next;

    /**
     * Construct a new {@code RestrictedVcfWindow} instance.
     * @param it an iterator that returns VCF records.
     * @throws IllegalArgumentException if {@code it.hasNext() == false}
     * @throws IllegalArgumentException if a format error is detected in
     * a VCF record
     * @throws NullPointerException if {@code it == null}
     */
    public RestrictedVcfWindow(SampleFileIt<? extends VcfEmission> it) {
        if (it.hasNext()==false) {
            throw new IllegalArgumentException("it.hasNext()==false");
        }
        this.it = it;
        this.window = new ArrayList<>();
        this.markers = null;
        this.overlap = 0;
        this.cumMarkerCnt = 0;
        this.next = it.next();
    }

    /**
     * Advances the sliding marker window, and returns the advanced window
     * as a {@code VcfEmission[]} object.
     * The returned array will have length {@code markers.nMarkers()}.
     * Markers not found in the data source will have {@code null}
     * entries in the returned array.
     *
     * @param nextMarkers the set of markers in the advanced window
     * @return the advanced marker window
     *
     * @throws IllegalArgumentException if {@code markers.nMarkers() == 0}
     * @throws IllegalArgumentException if any two of the specified markers
     * are on different chromosomes
     * @throws IllegalArgumentException if specified markers are
     * inconsistent with a sliding marker window
     * @throws IllegalArgumentException if the specified markers do not
     * advance the current marker window
     * @throws IllegalArgumentException if a format error is detected in a
     * VCF record
     * @throws IllegalArgumentException if the input data does not contain
     * any of the specified markers
     * @throws NullPointerException if {@code nextMarkers == null}
     */
    public VcfEmission[] advanceWindow(Markers nextMarkers) {
        checkMarkers(nextMarkers);
        advanceToCurrentChrom(nextMarkers);
        int fullOverlap = overlap(markers, nextMarkers);

        List<VcfEmission> newWindow = new ArrayList<>(nextMarkers.nMarkers());
        newWindow.addAll(window.subList(window.size() - fullOverlap, window.size()));
        this.overlap = countNonNull(newWindow);
        for (int j = fullOverlap, n=nextMarkers.nMarkers(); j<n; ++j) {
            Marker m = nextMarkers.marker(j);
            if (next!=null && next.marker().chromIndex()==m.chromIndex()) {
                while (next != null && next.marker().pos() < m.pos()) {
                    next = it.hasNext() ? it.next() : null;
                }
                while (next != null && next.marker().pos() == m.pos()
                        && next.marker().equals(m)==false) {
                    next = it.hasNext() ? it.next() : null;
                }
            }
            if (next != null && next.marker().equals(m)) {
                ++cumMarkerCnt;
                newWindow.add(next);
                next = it.hasNext() ? it.next() : null;
            }
            else {
                newWindow.add(null);
            }
        }
        this.markers = nextMarkers;
        this.window.clear();
        this.window.addAll(newWindow);
        if (countNonNull(newWindow) == 0) {
            missingMarkersErr(nextMarkers);
        }
        return window.toArray(new VcfEmission[0]);
    }

    private void checkMarkers(Markers mkrs) {
        if (mkrs.nMarkers()==0) {
            throw new IllegalArgumentException("markers do not advance window");
        }
        Marker start = mkrs.marker(0);
        Marker end = mkrs.marker(mkrs.nMarkers()-1);
        if (this.markers != null) {
            Marker m = this.markers.marker(this.markers.nMarkers()-1);
            if (m.chromIndex()==end.chromIndex() && m.pos()>=end.pos()) {
                String s = "markers do not advance window";
                throw new IllegalArgumentException(s);
            }
        }
        if (start.chromIndex() != end.chromIndex()) {
            String s = "inconsistent chromosomes:" + Const.nl
                    + start + Const.nl + end;
            throw new IllegalArgumentException(s);
        }
    }

    private void advanceToCurrentChrom(Markers markers) {
        int chromIndex = markers.marker(0).chromIndex();
        while (next!=null && next.marker().chromIndex()!=chromIndex) {
            next = it.hasNext() ? it.next() : null;
        }
    }

    private static int overlap(Markers prev, Markers next) {
        if (prev==null
                || prev.marker(0).chromIndex() != next.marker(0).chromIndex()) {
            return 0;
        }
        Marker startMarker = next.marker(0);
        int startPos = startMarker.pos();
        int index = prev.nMarkers() - 1; // index of first overlap marker
        while (index >= 0 && prev.marker(index).pos() > startPos) {
           --index;
        }
        while (index >= 0 && prev.marker(index).equals(startMarker)==false) {
            --index;
        }
        if (index < 0) {
            throw new IllegalArgumentException("inconsistent markers");
        }
        int overlap = prev.nMarkers() - index;
        for (int j = 1; j < overlap; ++j) {
            if (prev.marker(index + j).equals(next.marker(j))==false) {
                throw new IllegalArgumentException("inconsistent markers");
            }
        }
        return overlap;
    }

    private static <E> int countNonNull(List<E> list) {
        int cnt = 0;
        for (int j=0, n=list.size(); j<n; ++j) {
            if (list.get(j)!=null) {
                ++cnt;
            }
        }
        return cnt;
    }

    private static void missingMarkersErr(Markers markers) {
        StringBuilder sb = new StringBuilder(500);
        sb.append(Const.nl);
        sb.append("ERROR: Reference and target files have no markers in common"
                +  " in interval: ");
        sb.append(Const.nl);
        sb.append("       ");
        sb.append(interval(markers));
        sb.append(Const.nl);
        sb.append(Const.nl);
        sb.append("Common markers must have identical CHROM, POS, REF, and ALT"
                + " fields.");
        sb.append(Const.nl);
        sb.append("Exiting program.");
        sb.append(Const.nl);
        blbutil.Utilities.exit(sb.toString());
    }

    private static String interval(Markers markers) {
        Marker a = markers.marker(0);
        Marker b = markers.marker(markers.nMarkers()-1);
        assert a.chromIndex() == b.chromIndex();
        return a.chrom() + Const.colon + a.pos() + Const.hyphen + b.pos();
    }

    /**
     * Returns the file from which VCF records are read, or returns
     * {@code null} if the source is standard input.
     * @return the file from which VCF records are read, or
     * {@code null} if the source is standard input
     */
    public File file() {
        return it.file();
    }

    /**
     * Returns the list of samples.
     * @return the list of samples
     */
    public Samples samples() {
        return it.samples();
    }

    /**
     * Returns the number of samples.
     * @return the number of samples
     */
    public int nSamples() {
        return it.samples().nSamples();
    }

    /**
     * Returns the number of VCF records in the overlap between the current
     * window and the previous window.  Returns 0 if the current window
     * is the first window.
     *
     * @return the number of VCF records in the overlap between the current
     * window and the previous window
     */
    public int overlap() {
        return overlap;
    }

    /**
     * Returns the number of VCF records in the union of the current window
     * and all previous windows.
     *
     * @return the number of VCF records in the union of the current window
     * and all previous windows
     */
    public int cumMarkerCnt() {
        return cumMarkerCnt;
    }

    /**
     * Releases any I/O resources controlled by this object.
     */
    @Override
    public void close() {
        it.close();
    }

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(100);
        sb.append(this.getClass().toString());
        sb.append(" - next:");
        sb.append(next);
        return sb.toString();
    }
}
