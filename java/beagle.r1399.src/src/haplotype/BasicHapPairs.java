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

import java.util.List;
import vcf.Marker;
import vcf.Markers;

/**
 * <p>Class {@code BasicHapPairs} stores a list of haplotype pairs.
 * </p>
 * Instances of class {@code BasicHapPairs} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class BasicHapPairs implements HapPairs {

    private final Markers markers;
    private final HapPair[] hapPairs;
    private final boolean reverseMarkers;
    private final int lastMarker;

    /**
     * Constructs a {@code SimpleHapPairs} instance corresponding to
     * the specified list of haplotype pairs.
     * @param hapPairList a list of haplotype pairs.
     *
     * @throws IllegalArgumentException if
     * {@code hapPairList.isEmpty()==true}.
     * @throws NullPointerException if {@code hapPairList==null},
     * or if any element of {@code hapPairList} is null.
     */
    public BasicHapPairs(List<HapPair> hapPairList) {
        this(hapPairList, false);
    }

    /**
     * Constructs a {@code SimpleHapPairs} instance corresponding to
     * the specified list of haplotype pairs.
     * @param hapPairList a list of haplotype pairs.
     * @param reverseMarkers {@code true} if the marker order of the
     * specified haplotype pairs will be reversed, and {@code false}
     * otherwise.
     *
     * @throws IllegalArgumentException if
     * {@code hapPairList.isEmpty()==true}.
     * @throws NullPointerException if {@code hapPairList==null},
     * or if any element of {@code hapPairList} is null.
     */
    public BasicHapPairs(List<HapPair> hapPairList, boolean reverseMarkers) {
        if (hapPairList.isEmpty()) {
            throw new IllegalArgumentException("haps.isEmpy()==true");
        }
        Markers mrkrs = BasicSampleHapPairs.checkAndExtractMarkers(hapPairList);
        this.markers = (reverseMarkers ? mrkrs.reverse() : mrkrs);
        this.hapPairs = hapPairList.toArray(new HapPair[0]);
        this.reverseMarkers = reverseMarkers;
        this.lastMarker = markers.nMarkers() - 1;
    }


    @Override
    public byte allele1(int marker, int hapPair) {
        if (reverseMarkers) {
            marker = lastMarker - marker;
        }
        return hapPairs[hapPair].allele1(marker);
    }

    @Override
    public byte allele2(int marker, int hapPair) {
        if (reverseMarkers) {
            marker = lastMarker - marker;
        }
        return hapPairs[hapPair].allele2(marker);
    }

    @Override
    public byte allele(int marker, int haplotype) {
        if (reverseMarkers) {
            marker = lastMarker - marker;
        }
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
    public int idIndex(int hapPair) {
        return hapPairs[hapPair].idIndex();
    }

    /**
     * Returns a string representation of {@code this}.  The exact details
     * of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}.
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(10000);
        sb.append("BasicHapPairs: nHapPairs=");
        sb.append(this.nHapPairs());
        return sb.toString();
    }
}
