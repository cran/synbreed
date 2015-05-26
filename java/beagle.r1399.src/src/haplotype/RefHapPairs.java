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
import vcf.Marker;
import vcf.Markers;
import vcf.VcfEmission;

/**
 * <p>Class {@code RefHapPairs} stores a list of reference haplotype pairs.
 * </p>
 * Instances of class {@code RefHapPairs} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class RefHapPairs implements SampleHapPairs {

    private final Markers markers;
    private final Samples samples;
    private final VcfEmission[] phasedMarkers;

    /**
     * Constructs a {@code RefHapPairs} instance for the specified
     * reference genotype data.
     * @param markers the sequence of markers.
     * @param samples the sequence of samples.
     * @param phasedMarkers the sequence of per-marker genotype data.
     *
     * @throws IllegalArgumentException if
     * {@code markers.nMarkers()!=phasedMarkers.length}.
     * @throws IllegalArgumentException if
     * {@code phasedMarkers[k].samples().equals(samples)==false} for some
     * {@code 0<=k  && k<phasedMarkers.length}.
     * @throws IllegalArgumentException if
     * {@code phasedMarkers[k].marker().equals(markers.marker(k))==false}
     * for some {@code 0<=k && k<phasedMarkers.length}.
     * @throws IllegalArgumentException if
     * {@code phasedMarkers[k].isRefData()==false}
     * for some {@code 0<=k && k<phasedMarkers.length}.
     * @throws NullPointerException if
     * {@code markers==null || samples==null || phasedMarkers==null}.
     */
    public RefHapPairs(Markers markers, Samples samples,
            VcfEmission[] phasedMarkers) {
        checkPhasedMarkers(markers, samples, phasedMarkers);
        this.markers = markers;
        this.samples = samples;
        this.phasedMarkers = phasedMarkers.clone();
    }

    private static void checkPhasedMarkers(Markers markers, Samples samples,
            VcfEmission[] phasedMarkers) {
        if (markers.nMarkers()!=phasedMarkers.length) {
            String s = "markers.nMarkers()=" + markers.nMarkers()
                    + " phasedMarkers.length=" + phasedMarkers.length;
            throw new IllegalArgumentException(s);
        }
        for (int j=0; j<phasedMarkers.length; ++j) {
            if (phasedMarkers[j].samples().equals(samples)==false) {
                String s = "sample inconsistency at index " + j;
                throw new IllegalArgumentException(s);
            }
            if (phasedMarkers[j].marker().equals(markers.marker(j))==false) {
                String s = "marker inconsistency at index " + j;
                throw new IllegalArgumentException(s);
            }
            if (phasedMarkers[j].isRefData()==false) {
                String s = "non-reference data at marker index " + j;
                throw new IllegalArgumentException(s);
            }
        }
    }

    @Override
    public byte allele1(int marker, int hapPair) {
        return phasedMarkers[marker].allele1(hapPair);
    }

    @Override
    public byte allele2(int marker, int hapPair) {
        return phasedMarkers[marker].allele2(hapPair);
    }

    @Override
    public byte allele(int marker, int haplotype) {
        int hapPair = haplotype/2;
        if ((haplotype & 1)==0) {
            return phasedMarkers[marker].allele1(hapPair);
        }
        else {
            return phasedMarkers[marker].allele2(hapPair);
        }
    }

    @Override
    public int nMarkers() {
       return markers.nMarkers();
    }

    @Override
    public Markers markers() {
        return markers;
    }

    @Override
    public Marker marker(int marker) {
        return markers.marker(marker);
    }

    @Override
    public int nHaps() {
        return 2*samples.nSamples();
    }

    @Override
    public int nHapPairs() {
        return samples.nSamples();
    }

    @Override
    public int nSamples() {
        return samples.nSamples();
    }

    @Override
    public Samples samples() {
        return samples;
    }

    @Override
    public int idIndex(int hapPair) {
        return samples.idIndex(hapPair);
    }
}
