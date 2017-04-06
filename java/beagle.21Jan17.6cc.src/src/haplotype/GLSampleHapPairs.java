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
package haplotype;

import beagleutil.Samples;
import vcf.GL;
import vcf.Marker;
import vcf.Markers;

/**
 * <p>Class {@code GLSampleHapPairs} wraps a {@code GL} instance that stores
 * phased, non-missing genotypes.
 * </p>
 * Instances of class {@code GLSampleHapPairs} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class GLSampleHapPairs implements SampleHapPairs {

    private final GL gl;

    /**
     * Constructs a new {@code GLSampleHapPairs} instance from the
     * specified data.
     * @param gl phased, non-missing genotypes
     *
     * @throws IllegalArgumentException if {@code gl.isRefData() == false}
     * @throws NullPointerException if {@code gl == null}
     */
    public GLSampleHapPairs(GL gl) {
        if (gl.isRefData() == false) {
            throw new IllegalArgumentException("gl.isRefData()==false");
        }
        this.gl = gl;
    }

    @Override
    public Samples samples() {
        return gl.samples();
    }

    @Override
    public int nSamples() {
        return gl.nSamples();
    }

    @Override
    public int allele(int marker, int haplotype) {
        int sample = haplotype/2;
        if ( (haplotype & 1) ==  0) {
            return gl.allele1(marker, sample);
        }
        else {
            return gl.allele2(marker, sample);
        }
    }

    @Override
    public int allele1(int marker, int hapPair) {
            return gl.allele1(marker, hapPair);
    }

    @Override
    public int allele2(int marker, int hapPair) {
                    return gl.allele2(marker, hapPair);
    }

    @Override
    public int nMarkers() {
        return gl.nMarkers();
    }

    @Override
    public Markers markers() {
        return gl.markers();
    }

    @Override
    public Marker marker(int marker) {
        return gl.marker(marker);
    }

    @Override
    public int nHaps() {
        return 2*gl.nSamples();
    }

    @Override
    public int nHapPairs() {
        return gl.nSamples();
    }

    @Override
    public Samples samples(int hapPair) {
        if (hapPair < 0 || hapPair >= gl.nSamples()) {
            throw new IndexOutOfBoundsException(String.valueOf(hapPair));
        }
        return gl.samples();
    }

    @Override
    public int sampleIndex(int hapPair) {
        if (hapPair < 0 || hapPair >= gl.nSamples()) {
            throw new IndexOutOfBoundsException(String.valueOf(hapPair));
        }
        return hapPair;
    }

    @Override
    public int nAlleles(int marker) {
        return gl.marker(marker).nAlleles();
    }

    @Override
    public boolean storesNonMajorIndices(int marker) {
        return false;
    }

    @Override
    public int majorAllele(int marker) {
        String s = "this.storesNonMajorIndices(marker)==false";
        throw new UnsupportedOperationException(s);
    }

    @Override
    public int alleleCount(int marker, int allele) {
        String s = "this.storesNonMajorIndices(marker)==false";
        throw new UnsupportedOperationException(s);
    }

    @Override
    public int hapIndex(int marker, int allele, int copy) {
        String s = "this.storesNonMajorIndices(marker)==false";
        throw new UnsupportedOperationException(s);
    }
}
