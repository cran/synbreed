/*
 * Copyright (C) 2014 Brian L. Browning
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
package vcf;

import beagleutil.Samples;

/**
 * <p>Class {@code MaskedEndsGL} is a wrapper for a {@code GL}
 * instance that masks the genotype emission probabilities for a
 * user-specified number of starting and ending markers. The {@code gl()},
 * {@code allele1()}, and {@code allele2()} methods return
 * {@code 1.0f}, {@code -1}, and {@code -1} respectively if a genotype
 * emission probability for a marker is masked.
 * </p>
 * <p>Instances of class {@code MaskedEndsGL} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class MaskedEndsGL implements GL {

    private final GL gl;
    private final int start;
    private final int end;

    /**
     * Constructs a new {@code MaskedEndsGL} instance.
     * @param gl genotype emission probabilities for all markers
     * @param start the starting marker index (inclusive) of the markers
     * whose genotype emission probabilities are not masked
     * @param end the ending marker index (exclusive) of the markers
     * whose genotype emission probabilities are not masked
     * @throws IllegalArgumentException if
     * {@code start < 0 || start > end || end > gl.nMarkers()}
     * @throws NullPointerException if {@code gl == null}
     */
    public MaskedEndsGL(GL gl, int start, int end) {
        if (start<0 || start>end || end>gl.nMarkers()) {
            String s = "start=" + start + " end=" + end;
            throw new IllegalArgumentException(s);
        }
        this.gl = gl;
        this.start = start;
        this.end = end;
    }

    @Override
    public boolean isRefData() {
        if ((start>0 || end<gl.nMarkers()) && start<end) {
            return false;
        }
        else {
            return gl.isRefData();
        }
    }

    private void checkMarkerAndSample(int marker, int sample) {
        if (marker<0 || marker>=gl.nMarkers()) {
            throw new IndexOutOfBoundsException("marker: " + marker);
        }
        if (sample<0 || sample>=gl.nSamples()) {
            throw new IndexOutOfBoundsException("sample: " + sample);
        }
    }

    private void checkAllele(int marker, int allele) {
        if (allele<0 || allele>=gl.marker(marker).nAlleles()) {
            String s = "marker=" + marker + " allele=" + allele;
            throw new IndexOutOfBoundsException(s);
        }
    }

    @Override
    public float gl(int marker, int sample, int allele1, int allele2) {
        if (marker<start || marker>=end) {
            checkMarkerAndSample(marker, sample);
            checkAllele(marker, allele1);
            checkAllele(marker, allele2);
            return 1.0f;
        }
        else {
            return gl.gl(marker, sample, allele1, allele2);
        }
    }

    @Override
    public boolean isPhased(int marker, int sample) {
        if (marker<start || marker>=end) {
            checkMarkerAndSample(marker, sample);
            return false;
        }
        else {
            return gl.isPhased(marker, sample);
        }
    }

    @Override
    public int allele1(int marker, int sample) {
        if (marker<start || marker>=end) {
            checkMarkerAndSample(marker, sample);
            return -1;
        }
        else {
            return gl.allele1(marker, sample);
        }
    }

    @Override
    public int allele2(int marker, int sample) {
        if (marker<start || marker>=end) {
            checkMarkerAndSample(marker, sample);
            return -1;
        }
        else {
            return gl.allele2(marker, sample);
        }
    }

    @Override
    public int allele(int marker, int hap) {
        if (marker<start || marker>=end) {
            checkMarkerAndSample(marker, hap/2);
            return -1;
        }
        else {
            return gl.allele(marker, hap);
        }
    }

    @Override
    public Marker marker(int marker) {
        return gl.marker(marker);
    }

    @Override
    public Markers markers() {
        return gl.markers();
    }

    @Override
    public int nMarkers() {
        return gl.nMarkers();
    }

    @Override
    public int nHaps() {
        return gl.nHaps();
    }

    @Override
    public int nSamples() {
        return gl.nSamples();
    }

    @Override
    public Samples samples() {
        return gl.samples();
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(10000);
        sb.append(this.getClass().toString());
        sb.append(": nSamples=");
        sb.append(this.nSamples());
        return sb.toString();
    }
}
