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
import java.util.List;
import vcf.Marker;
import vcf.Markers;

/**
 * <p>Class {@code BasicHapPairs} represents a list of haplotype pairs.
 * Each haplotype pair is guaranteed to have two non-missing
 * alleles at each marker.
 * </p>
 * Instances of class {@code BasicHapPairs} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class BasicHapPairs implements HapPairs {

    private final Markers markers;
    private final HapPair[] hapPairs;

    /**
     * Constructs a new {@code BasicHapPairs} instance corresponding to
     * the specified list of haplotype pairs.
     * @param hapPairList a list of haplotype pairs
     *
     * @throws IllegalArgumentException if
     * {@code hapPairList.isEmpty() == true}
     * @throws NullPointerException if
     * {@code (hapPairList == null || hapPairList[j] == null)} for any {@code j}
     * satisfying {@code (0 <= j && j < hapPairsList.size())}
     */
    public BasicHapPairs(List<HapPair> hapPairList) {
        if (hapPairList.isEmpty()) {
            throw new IllegalArgumentException("haps.isEmpy()==true");
        }
        this.markers = BasicSampleHapPairs.checkAndExtractMarkers(hapPairList);
        this.hapPairs = hapPairList.toArray(new HapPair[0]);
    }

    @Override
    public int allele1(int marker, int hapPair) {
        return hapPairs[hapPair].allele1(marker);
    }

    @Override
    public int allele2(int marker, int hapPair) {
        return hapPairs[hapPair].allele2(marker);
    }

    @Override
    public int allele(int marker, int haplotype) {
        int hapPair = haplotype / 2;
        if ((haplotype & 1) == 0) {
            return hapPairs[hapPair].allele1(marker);
        } else {
            return hapPairs[hapPair].allele2(marker);
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
        return 2*hapPairs.length;
    }

    @Override
    public int nHapPairs() {
        return hapPairs.length;
    }

    @Override
    public Samples samples(int hapPair) {
        return hapPairs[hapPair].samples();
    }

    @Override
    public int sampleIndex(int hapPair) {
        return hapPairs[hapPair].sampleIndex();
    }

    /**
     * Returns a string representation of {@code this}.  The exact details
     * of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(10000);
        sb.append('[');
        sb.append(this.getClass().toString());
        sb.append(": nHapPairs=");
        sb.append(this.nHapPairs());
        sb.append(']');
        return sb.toString();
    }
}
