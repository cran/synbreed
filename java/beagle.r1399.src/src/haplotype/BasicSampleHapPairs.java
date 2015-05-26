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
import java.util.Comparator;
import java.util.List;
import vcf.Marker;
import vcf.Markers;

/**
 * <p>Class {@code SampleHaps} stores a list of samples and a
 * haplotype pair for each sample.
 * </p>
 * Class {@code SampleHaps} is immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class BasicSampleHapPairs implements SampleHapPairs {

    private final Samples samples;
    private final Markers markers;
    private final HapPair[] hapPairs;
    private final boolean reverseMarkers;
    private final int lastMarker;

    /**
     * Constructs a {@code SimpleSampleHapPairs} instance corresponding to
     * the specified data.
     * @param samples a list of samples.
     * @param hapPairList a list of haplotype pairs corresponding to the
     * specified list of samples.
     *
     * @throws IllegalArgumentException if
     * {@code hapPairList.isEmpty()==true}.
     * @throws IllegalArgumentException if
     * {@code hapPairList.get(j).markers()!=hapPairList.get(k).markers()}
     * for any {@code <= j && j<k && k<hapPairList.size}.
     * @throws IllegalArgumentException if
     * {@code samples.idIndex(j)!=hapPairList.get(j).idIndex()}
     * for some {@code 0<=j && j<hapPairList.size()}.
     * @throws NullPointerException if
     * {@code sample==null || hapPairList==null},
     * or if any element of {@code hapPairList} is null.
     * @see beagleutil.SampleIds
     */
    public BasicSampleHapPairs(Samples samples, List<HapPair> hapPairList) {
        this(samples, hapPairList, false);
    }

    /**
     * Constructs a {@code SimpleSampleHapPairs} instance corresponding to
     * the specified data.
     * @param samples a list of samples.
     * @param hapPairList a list of haplotype pairs corresponding to the
     * specified list of samples.
     * @param reverseMarkers {@code true} if the marker order of the
     * specified haplotype pairs will be reversed, and {@code false}
     * otherwise.
     *
     * @throws IllegalArgumentException if
     * {@code hapPairList.isEmpty()==true}.
     * @throws IllegalArgumentException if
     * {@code hapPairList.get(j).markers().equals(hapPairList.get(k).markers())==false}
     * for any {@code 0<=j && j<k && k<hapPairList.size}.
     * @throws IllegalArgumentException if
     * {@code samples.idIndex(j)!=hapPairList.get(j).idIndex()}
     * for some {@code 0<=j && j<hapPairList.size()}.
     * @throws NullPointerException if
     * {@code samples==null || hapPairList==null},
     * or if any element of {@code hapPairList} is null.
     * @see beagleutil.SampleIds
     */
    public BasicSampleHapPairs(Samples samples, List<HapPair> hapPairList,
            boolean reverseMarkers) {
        if (hapPairList.isEmpty()) {
            throw new IllegalArgumentException("haps.isEmpy()==true");
        }
        checkSamples(samples, hapPairList);
        Markers mkrs = BasicSampleHapPairs.checkAndExtractMarkers(hapPairList);
        this.samples = samples;
        this.markers = (reverseMarkers ? mkrs.reverse() : mkrs);
        this.hapPairs = hapPairList.toArray(new HapPair[0]);
        this.reverseMarkers = reverseMarkers;
        this.lastMarker = markers.nMarkers() - 1;
    }

    private void checkSamples(Samples samples,
            List<HapPair> hapPairList) {
        if (samples.nSamples()!= hapPairList.size()) {
            String s = "samples and hapPairList are inconsistent";
            throw new IllegalArgumentException(s);
        }
        for (int j=0, n=samples.nSamples(); j<n; ++j) {
            if (samples.idIndex(j)!=hapPairList.get(j).idIndex()) {
                String s = "samples and hapPairList are inconsistent";
                throw new IllegalArgumentException(s);
            }
        }
    }

    /**
     * Checks that all haplotype pairs have the same markers reference
     * and returns the shared markers.
     * @param hapPairList a list of haplotype pairs.
     * @return the markers shared by all the specified haplotype pairs.
     * @throws IllegalArgumentException if
     * {@code haps.get(j).markers().equals(haps.get(k).markers())==false}
     * for any {@code 0<=j && j<k && k<haps.size()}.
     * @throws NullPointerException if {@code hapPairList==null},
     * or if any element of {@code hapPairList} is null.
     */
    public static Markers checkAndExtractMarkers(
            List<HapPair> hapPairList) {
        if (hapPairList.isEmpty()) {
            return new Markers(new Marker[0]);
        }
        else {
            Markers m = hapPairList.get(0).markers();
            for (int j=1, n=hapPairList.size(); j<n; ++j) {
                if (hapPairList.get(j).markers().equals(m)==false) {
                    String s = "Markers are not consistent";
                    throw new IllegalArgumentException(s);
                }
            }
            return m;
        }
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
        int pairIndex = haplotype/2;
        if ((haplotype & 1)==0) {
            return hapPairs[pairIndex].allele1(marker);
        }
        else {
            return hapPairs[pairIndex].allele2(marker);
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
    public int nSamples() {
        return samples.nSamples();
    }

    @Override
    public Samples samples() {
        return samples;
    }

    @Override
    public int idIndex(int hapPair) {
        return hapPairs[hapPair].idIndex();
    }

    /**
     * Returns a {@code Comparator<HapPairInterface>}
     * whose {@code compare(hp1, hp2)} method returns -1, 0, or 1
     * depending on whether {@code hp1.idIndex()} is less than, equal,
     * or greater than {@code hp2.idIndex()}.
     * @return a {@code Comparator<HapPairInterface>}
     * whose {@code compare(hp1, hp2)} method compares two
     * haplotype pairs for order.
     */
    public static Comparator<HapPair> hapsComparator() {
        return new Comparator<HapPair>() {
            @Override
            public int compare(HapPair hp1, HapPair hp2) {
                int i1 = hp1.idIndex();
                int i2 = hp2.idIndex();
                if (i1==i2) {
                    return 0;
                }
                else {
                    return (i1 < i2) ? -1 : 1;
                }
            }
        } ;
    }

    /**
     * Returns a {@code Comparator<HapPairInterface>}
     * whose {@code compare(hp1, hp2)} method returns -1, 0, or 1
     * depending on whether {@code samples.index(hp1.idIndex())} is
     * less than, equal, or greater than
     * {@code samples.index(hp2.idIndex())}.
     * @param samples the list of samples used to compare {@code HapsPair}
     * objects.
     * @return a {@code Comparator<HapPairInterface>}
     * whose {@code compare(hp1, hp2)} method compares two
     * haplotype pairs for order.
     * @throws NullPointerException if {@code samples==null}.
     * @see beagleutil.SampleIds
     */
    public static Comparator<HapPair> hapsComparator(
            final Samples samples) {
        if (samples==null) {
            throw new NullPointerException("samples==null");
        }
        return new Comparator<HapPair>() {
            @Override
            public int compare(HapPair hp1, HapPair hp2) {
                int i1 = samples.index(hp1.idIndex());
                int i2 = samples.index(hp2.idIndex());
                if (i1==i2) {
                    return 0;
                }
                else {
                    return (i1 < i2) ? -1 : 1;
                }
            }
        } ;
    }
}
