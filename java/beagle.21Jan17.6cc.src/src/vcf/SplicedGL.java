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
import haplotype.SampleHapPairs;

/**
 * <p>Class {@code SplicedGL} represents genotype emission probabilities
 * for a set of samples. The genotype emission probabilities are determined
 * by a {@code SampleHapPairs} instance for the initial markers, and are
 * determined by a {@code GL} instance for the remaining markers.
 * The {@code isRefData()} method of the {@code SplicedGL} class
 * returns the same value as the {@code isRefData()} method of
 * the {@code GL} instance.
 * </p>
 * <p>Instances of class {@code SplicedGL} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class SplicedGL implements GL {

    private final int overlap;
    private final SampleHapPairs haps;
    private final GL gl;

    /**
     * Constructs a new {@code SplicedGL} instance.
     * @param haps sample haplotype pairs for the initial
     * markers
     * @param gl genotype emission probabilities for all markers
     * @throws IllegalArgumentException if
     * {@code haps.nMarkers() >= gl.nMarkers()}
     * @throws IllegalArgumentException if
     * {@code haps.marker(j).equals(gl.marker(j)) == false} for any {@code j}
     * satisfying {@code 0 <= j && j < haps.nMarkers()}
     * @throws IllegalArgumentException if
     * {@code haps.samples().equals(gl.samples()) == false}
     * @throws NullPointerException if {@code haps == null || gl == null}
     */
    public SplicedGL(SampleHapPairs haps, GL gl) {
        if (haps.nMarkers()>=gl.nMarkers()) {
            throw new IllegalArgumentException("inconsistent markers");
        }
        for (int j=0, n=haps.nMarkers(); j<n; ++j) {
            if (haps.marker(j).equals(gl.marker(j))==false) {
                throw new IllegalArgumentException("inconsistent markers");
            }
        }
        if (haps.samples().equals(gl.samples())==false) {
            throw new IllegalArgumentException("inconsistent samples");
        }
        this.overlap = haps.nMarkers();
        this.haps = haps;
        this.gl = gl;
    }

    @Override
    public boolean isRefData() {
        return gl.isRefData();
    }

    @Override
    public float gl(int marker, int sample, int allele1, int allele2) {
        if (marker<overlap) {
            int a1 = haps.allele1(marker, sample);
            int a2 = haps.allele2(marker, sample);
            return (allele1==a1 && allele2==a2) ? 1.0f : 0.0f;
        }
        else {
            return gl.gl(marker, sample, allele1, allele2);
        }
    }

    @Override
    public boolean isPhased(int marker, int sample) {
        if (marker<overlap) {
            return true;
        }
        else {
            return gl.isPhased(marker, sample);
        }
    }

    @Override
    public int allele1(int marker, int sample) {
        if (marker<overlap) {
            return haps.allele1(marker, sample);
        }
        else {
            return gl.allele1(marker, sample);
        }
    }

    @Override
    public int allele2(int marker, int sample) {
        if (marker<overlap) {
            return haps.allele2(marker, sample);
        }
        else {
            return gl.allele2(marker, sample);
        }
    }

    @Override
    public int allele(int marker, int hap) {
        if (marker<overlap) {
            return haps.allele(marker, hap);
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
        sb.append("SplicedGL: nSamples=");
        sb.append(this.nSamples());
        return sb.toString();
    }
}
