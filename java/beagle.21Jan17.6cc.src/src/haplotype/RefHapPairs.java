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
 * <p>Class {@code RefHapPairs} stores a list of samples and a
 * haplotype pair for each sample.
 * </p>
 * <p>Instances of class {@code RefHapPairs} are immutable.<p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class RefHapPairs implements SampleHapPairs {

    private final Markers markers;
    private final Samples samples;
    private final VcfEmission[] refVcfRecs;

    /**
     * Constructs a new {@code RefHapPairs} instance.
     * @param markers the sequence of markers
     * @param samples the sequence of samples
     * @param refVcfRecs the sequence of per-marker genotype data
     *
     * @throws IllegalArgumentException if
     * {@code markers.nMarkers() != refVcfRecs.length}
     * @throws IllegalArgumentException if
     * {@code refVcfRecs[k].samples().equals(samples) == false} for any
     * {@code k} satisfying {@code 0 <= k  && k < refVcfRecs.length}
     * @throws IllegalArgumentException if
     * {@code refVcfRecs[k].marker().equals(markers.marker(k)) == false}
     * for any {@code k} satisfying {@code 0 <= k && k < refVcfRecs.length}
     * @throws IllegalArgumentException if
     * {@code refVcfRecs[k].isRefData() == false} for any {@code k}
     * satisfying {@code 0 <= k && k < refVcfRecs.length}
     * @throws NullPointerException if
     * {@code markers == null || samples == null || refVcfRecs == null
     * || refVcfRecs[k] == null} for any {@code k} satisfying
     * {@code 0 <= k && k <= refVcfRecs.length}
     */
    public RefHapPairs(Markers markers, Samples samples,
            VcfEmission[] refVcfRecs) {
        checkPhasedMarkers(markers, samples, refVcfRecs);
        this.markers = markers;
        this.samples = samples;
        this.refVcfRecs = refVcfRecs.clone();
    }

    private static void checkPhasedMarkers(Markers markers, Samples samples,
            VcfEmission[] refVcfRecs) {
        if (markers.nMarkers()!=refVcfRecs.length) {
            String s = "markers.nMarkers()=" + markers.nMarkers()
                    + " refVcfRecs.length=" + refVcfRecs.length;
            throw new IllegalArgumentException(s);
        }
        for (int j=0; j<refVcfRecs.length; ++j) {
            if (refVcfRecs[j].samples().equals(samples)==false) {
                String s = "sample inconsistency at index " + j;
                throw new IllegalArgumentException(s);
            }
            if (refVcfRecs[j].marker().equals(markers.marker(j))==false) {
                String s = "marker inconsistency at index " + j;
                throw new IllegalArgumentException(s);
            }
            if (refVcfRecs[j].isRefData()==false) {
                String s = "non-reference data at marker index " + j;
                throw new IllegalArgumentException(s);
            }
        }
    }

    @Override
    public int allele1(int marker, int hapPair) {
        return refVcfRecs[marker].allele1(hapPair);
    }

    @Override
    public int allele2(int marker, int hapPair) {
        return refVcfRecs[marker].allele2(hapPair);
    }

    @Override
    public int allele(int marker, int haplotype) {
        int hapPair = haplotype/2;
        if ((haplotype & 1)==0) {
            return refVcfRecs[marker].allele1(hapPair);
        }
        else {
            return refVcfRecs[marker].allele2(hapPair);
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
    public Samples samples(int hapPair) {
        if (hapPair < 0 || hapPair >= samples.nSamples()) {
            throw new IndexOutOfBoundsException(String.valueOf(hapPair));
        }
        return samples;
    }

    @Override
    public int sampleIndex(int hapPair) {
        if (hapPair < 0 || hapPair >= samples.nSamples()) {
            throw new IndexOutOfBoundsException(String.valueOf(hapPair));
        }
        return hapPair;
    }

    @Override
    public int nAlleles(int marker) {
        return refVcfRecs[marker].nAlleles();
    }

    @Override
    public boolean storesNonMajorIndices(int marker) {
        return refVcfRecs[marker].storesNonMajorIndices();
    }

    @Override
    public int majorAllele(int marker) {
        return refVcfRecs[marker].majorAllele();
    }

    @Override
    public int alleleCount(int marker, int allele) {
        return refVcfRecs[marker].alleleCount(allele);
    }

    @Override
    public int hapIndex(int marker, int allele, int copy) {
        return refVcfRecs[marker].hapIndex(allele, copy);
    }
}
