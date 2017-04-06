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
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import vcf.Marker;
import vcf.Markers;

/**
 * <p>Class {@code BasicSampleHapPairs} stores a list of samples and a
 * haplotype pair for each sample.
 * </p>
 * <p>Instance of class {@code BasicSampleHapPairs} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class BasicSampleHapPairs implements SampleHapPairs {

    private final Markers markers;
    private final Samples samples;
    private final HapPair[] hapPairs;

    /**
     * Constructs a new {@code BasicSampleHapPairs} instance.
     * @param samples a list of samples
     * @param hapPairList a list of haplotype pairs corresponding to the
     * specified list of samples
     *
     * @throws IllegalArgumentException if
     * {@code hapPairList.isEmpty() == true}
     * @throws IllegalArgumentException if
     * {@code hapPairList.get(j).markers().equals(hapPairList.get(k).markers())
     * == false}
     * for any indices {@code j, k} satisfying
     * {@code 0 <= j && j < k && k < hapPairList.size()}
     * @throws IllegalArgumentException if the list of samples does not
     * match the list of samples determined by {@code hapPairList}
     * @throws NullPointerException if {@code samples == null}
     * @throws NullPointerException if
     * {@code (hapPairList == null || hapPairList(j) == null)}
     * for any {@code j} satisfying {@code (0 <= j && j < hapPairList.size())}
     */
    public BasicSampleHapPairs(Samples samples, List<HapPair> hapPairList) {
        if (hapPairList.isEmpty()) {
            throw new IllegalArgumentException("haps.isEmpy()==true");
        }
        Collections.sort(hapPairList, hapsComparator(samples));
        checkSamples(samples, hapPairList);
        this.markers = checkAndExtractMarkers(hapPairList);
        this.samples = samples;
        this.hapPairs = hapPairList.toArray(new HapPair[0]);
    }

    private void checkSamples(Samples samples, List<HapPair> hapPairs) {
        if (samples.nSamples()!= hapPairs.size()) {
            throw new IllegalArgumentException("inconsistent samples");
        }
        for (int j=0, n=hapPairs.size(); j<n; ++j) {
            if (samples.equals(hapPairs.get(j).samples())==false) {
                HapPair hp = hapPairs.get(j);
                int i1 = samples.idIndex(j);
                int i2 = hp.samples().idIndex(hp.sampleIndex());
                if (i1 != i2) {
                    throw new IllegalArgumentException("inconsistent samples");
                }
            }
        }
    }

    /**
     * Checks that all haplotype pairs have alleles for the same list of
     * markers, and returns the list of markers.
     * @param hapPairList a list of haplotype pairs
     * @return the list of markers shared by the specified haplotype pairs
     * @throws IllegalArgumentException if
     * {@code hapPiarList.get(j).markers().equals(hapPairList.get(k).markers())
     * == false}
     * for any indices {@code j, k} satisfying
     * {@code 0 <= j && j < k && k < hapPairList.size()}
     * @throws NullPointerException if
     * {@code hapPairList == null || hapPairList(j) == null}
     * for any {@code j} satisfying {@code 0 <= j && j < hapPairList.size()}
     */
    static Markers checkAndExtractMarkers(List<HapPair> hapPairList) {
        if (hapPairList.isEmpty()) {
            return Markers.create(new Marker[0]);
        }
        else {
            Markers m = hapPairList.get(0).markers();
            for (int j=1, n=hapPairList.size(); j<n; ++j) {
                if (hapPairList.get(j).markers().equals(m)==false) {
                    throw new IllegalArgumentException("inconsistent markers");
                }
            }
            return m;
        }
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
    public Samples samples(int hapPair) {
        if (hapPair < 0 || hapPair >= hapPairs.length) {
            throw new IndexOutOfBoundsException(String.valueOf(hapPair));
        }
        return samples;
    }

    @Override
    public int sampleIndex(int hapPair) {
        if (hapPair < 0 || hapPair >= hapPairs.length) {
            throw new IndexOutOfBoundsException(String.valueOf(hapPair));
        }
        return hapPair;
    }

    @Override
    public int nAlleles(int marker) {
        return markers.marker(marker).nAlleles();
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

    /**
     * Returns a {@code Comparator<HapPairInterface>}
     * whose {@code compare(hp1, hp2)} method returns -1, 0, or 1
     * depending on whether {@code samples.index(hp1.idIndex())} is
     * less than, equal, or greater than
     * {@code samples.index(hp2.idIndex())}.
     * @param samples the list of samples used to compare {@code HapsPair}
     * objects
     * @return a {@code Comparator<HapPairInterface>}
     * whose {@code compare(hp1, hp2)} method compares two
     * haplotype pairs for order
     * @throws NullPointerException if {@code samples == null}
     */
    private static Comparator<HapPair> hapsComparator(final Samples samples) {
        if (samples==null) {
            throw new NullPointerException("samples==null");
        }
        return (HapPair hp1, HapPair hp2) -> {
            int id1 = hp1.samples().idIndex(hp1.sampleIndex());
            int id2 = hp2.samples().idIndex(hp2.sampleIndex());
            int i1 = samples.index(id1);
            int i2 = samples.index(id2);
            if (i1 == -1 || i2 == -1) {
                String id;
                if (i1 == -1) {
                    id = hp1.samples().id(hp1.sampleIndex());
                }
                else {
                    id = hp2.samples().id(hp2.sampleIndex());
                }
                String s = "samples do not contain: " + id;
                throw new IllegalArgumentException(s);
            }
            if (i1==i2) {
                return 0;
            }
            else {
                return (i1 < i2) ? -1 : 1;
            }
        } ;
    }
}
