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
import blbutil.Const;
import java.util.BitSet;
import vcf.Marker;
import vcf.Markers;

/**
 * <p>Class {@code BitHapPair} represents a pair of haplotypes for a sample.
 * The class stores alleles using {@code java.util.BitSet} objects.
 * </p>
 * Instances of class {@code BitHapPair} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class BitHapPair implements HapPair {

    private final Markers markers;
    private final Samples samples;
    private final int sampleIndex;
    private final BitSet alleles1;
    private final BitSet alleles2;

    /**
     * Constructs a new {@code BitHapPair} instance.
     * @param markers the sequence of markers
     * @param samples the list of samples
     * @param sampleIndex the sample index
     * @param alleles1 the sequence of allele indices for the first haplotype
     * @param alleles2 the sequence of alleles indices for the second haplotype
     *
     * @throws IllegalArgumentException if
     * {@code alleles1.length != markers.nMarkers()
     * || alleles2.length != markers.nMarkers()}
     * @throws IllegalArgumentException if {@code alleles1[k] < 0 ||
     * allele1[k] >= markers.marker(k).nAlleles()} for some {@code k} satisfying
     * {@code 0 <= k && k < markers.nMarkers()}
     * @throws IllegalArgumentException if {@code alleles2[k] < 0 ||
     * allele2[k] >= markers.marker(k).nAlleles()} for some {@code k} satisfying
     * {@code 0 <= k && k < markers.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code sampleIndex < 0 || sampleIndex >= samples.nSamples()}
     * @throws NullPointerException if
     * {@code marker == null || samples == null || alleles1 == null
     * || allele2 == null}
     */
    public BitHapPair(Markers markers, Samples samples, int sampleIndex,
            int[] alleles1, int[] alleles2) {
        if (alleles1.length != markers.nMarkers()
                || alleles2.length != markers.nMarkers()) {
            throw new IllegalArgumentException("inconsistent markers");
        }
        if (sampleIndex < 0 || sampleIndex >= samples.nSamples()) {
            throw new IndexOutOfBoundsException(String.valueOf(sampleIndex));
        }
        this.markers = markers;
        this.samples = samples;
        this.sampleIndex = sampleIndex;
        this.alleles1 = toBitSet(markers, alleles1);
        this.alleles2 = toBitSet(markers, alleles2);
    }

    private static BitSet toBitSet(Markers markers, int[] alleles) {
        int index = 0;
        BitSet bs = new BitSet(markers.sumHaplotypeBits());
        for (int k=0; k<alleles.length; ++k) {
            int allele = alleles[k];
            if (allele < 0 || allele >= markers.marker(k).nAlleles()) {
                String s = "allele \"" + allele + "\" out of bounds for marker: "
                        + markers.marker(k);
                throw new IllegalArgumentException(s);
            }
            int mask = 1;
            int nBits = markers.sumHaplotypeBits(k+1) - markers.sumHaplotypeBits(k);
            for (int l=0; l<nBits; ++l) {
                boolean b = (allele & mask)==mask;
                bs.set(index++, b);
                mask <<= 1;
            }
        }
        return bs;
    }

    @Override
    public int allele1(int marker) {
        return allele(alleles1, marker);
    }

    @Override
    public int allele2(int marker) {
        return allele(alleles2, marker);
    }

    private int allele(BitSet bitset, int marker) {
        int start = markers.sumHaplotypeBits(marker);
        int end = markers.sumHaplotypeBits(marker+1);
        if (end==(start+1)) {
            return bitset.get(start) ? 1 : 0;
        }
        int allele = 0;
        int mask = 1;
        for (int j=start; j<end; ++j) {
            if (bitset.get(j)) {
                allele += mask;
            }
            mask <<= 1;
        }
        return allele;
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
    public int nMarkers() {
        return markers.nMarkers();
    }

    @Override
    public Samples samples() {
        return samples;
    }

    @Override
    public int sampleIndex() {
        return sampleIndex;
    }

    /**
     * Returns a string representation of {@code this}.  The
     * exact details of the representation are unspecified and subject
     * to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("sampleIndex=");
        sb.append(sampleIndex);
        sb.append(Const.nl);
        sb.append(alleles1);
        sb.append(Const.nl);
        sb.append(alleles2);
        return sb.toString();
    }
}
