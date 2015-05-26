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
import java.util.ArrayList;
import java.util.List;
import vcf.Marker;
import vcf.Markers;

/**
 * Class {@code SampleHapPairsAligner} revises a {@code SampleHapPairs}
 * instance to be consistent with a {@code SampleHapPairs} instance
 * for the previous marker window.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class SampleHapPairsSplicer {

    private SampleHapPairsSplicer() {
        // private constructor prevents instantiation
    }

    /**
     * Returns a {@code SampleHapPairs} instance that has the same haplotype
     * pairs as {@code prev} for the first
     * {@code nextSplice} markers of {@code next.markers()}, and
     * the same haplotype pairs as {@code next} for the remaining
     * markers of {@code next.markers()}.  The order of the haplotypes
     * in the haplotype pairs copied from {@code next} is reversed when
     * necessary for {@code prev} and {@code next} to have an identical
     * ordered heterozygote genotype near the splice point.
     *
     * @param prev sample haplotype pairs for the previous marker window.
     * @param next sample haplotype pairs for the current marker window.
     * @param nextStart the first marker index in the data copied from
     * {@code next}.
     * @param overlap the number of overlapping markers for the
     * {@code prev.markers()} and {@code next.markers()}
     * @return a {@code SampleHapPairs} instance that has the same markers and
     * samples as {@code next}.
     *
     * @throws IllegalArgumentException if
     * {@code prev.samples().equals(next.samples())==false}
     * @throws IllegalArgumentException if
     * {@code overlap<0  || overlap>prev.nMarkers() || overlap>next.nMarkers()}
     * @throws IllegalArgumentException if {@code nextStart>overlap}
     * @throws IllegalArgumentException if
     * {@code prev.marker(prev.nMarker()-overlap+j).equals(next.marker(j))==false}
     * for any {@code j} satisfying {@code 0<j && j<overlap}
     */
    public static SampleHapPairs spliceNext(SampleHapPairs prev,
            SampleHapPairs next, int nextStart, int overlap) {
        checkSamples(prev, next);
        checkOverlap(prev, next, nextStart, overlap);
        if (overlap==0) {
            return next;
        }
        return new SplicedHapPairs(prev, next, nextStart, overlap);
    }

    private static void checkSamples(SampleHapPairs prev,
            SampleHapPairs next) {
        if (next.samples().equals(prev.samples()) == false) {
            throw new IllegalArgumentException("inconsistent samples");
        }
    }

    private static void checkOverlap(SampleHapPairs prev, SampleHapPairs next,
            int nextStart, int overlap) {
        if (overlap<0  || overlap>prev.nMarkers() || overlap>next.nMarkers()) {
            throw new IllegalArgumentException("overlap: " + overlap);
        }
        if (nextStart > overlap) {
            String s = "nextStart=" + nextStart + " overlap=" + overlap;
            throw new IllegalArgumentException(s);
        }
        int n = prev.nMarkers() - overlap;
        for (int j=0; j<overlap; ++j) {
            Marker prevMarker = prev.marker(n+j);
            Marker nextMarker = next.marker(j);
            if (prevMarker.equals(nextMarker)==false) {
                throw new IllegalArgumentException("inconsistent markers");
            }
        }
    }

    private static class SplicedHapPairs implements SampleHapPairs {

        private final SampleHapPairs splicedPrev;
        private final SampleHapPairs next;
        private final boolean[] hapPairIsSwitched;

        private SplicedHapPairs(SampleHapPairs prev, SampleHapPairs next,
                int nextStart, int overlap) {
            int n = prev.nMarkers();
            assert overlap<=n;
            this.splicedPrev = restrict(prev, n-overlap, n-overlap + nextStart);
            this.next = next;
            this.hapPairIsSwitched = switchHapOrder(prev, next, nextStart,
                    overlap);
        }

        @Override
        public byte allele1(int marker, int hapPair) {
            if (marker<splicedPrev.nMarkers()) {
                return splicedPrev.allele1(marker, hapPair);
            }
            else {
                if (hapPairIsSwitched[hapPair]) {
                    return next.allele2(marker, hapPair);
                } else {
                    return next.allele1(marker, hapPair);
                }
            }
        }

        @Override
        public byte allele2(int marker, int hapPair) {
            if (marker<splicedPrev.nMarkers()) {
                return splicedPrev.allele2(marker, hapPair);
            }
            else {
                if (hapPairIsSwitched[hapPair]) {
                    return next.allele1(marker, hapPair);
                } else {
                    return next.allele2(marker, hapPair);
                }
            }
        }

        @Override
        public byte allele(int marker, int haplotype) {
            int hapPair = haplotype / 2;
            if ((haplotype & 1) == 0) {
                return this.allele1(marker, hapPair);
            } else {
                return this.allele2(marker, hapPair);
            }
        }

        @Override
        public Markers markers() {
            return next.markers();
        }

        @Override
        public Marker marker(int marker) {
            return next.marker(marker);
        }

        @Override
        public int nHaps() {
            return 2*nSamples();
        }

        @Override
        public int nHapPairs() {
            return next.nSamples();
        }

        @Override
        public Samples samples() {
            return next.samples();
        }

        @Override
        public int idIndex(int hapPair) {
            return next.idIndex(hapPair);
        }

        @Override
        public int nSamples() {
            return next.nSamples();
        }

        @Override
        public int nMarkers() {
            return next.nMarkers();
        }
    }

    private static SampleHapPairs restrict(SampleHapPairs shp, int start,
            int end) {
        int size = end - start;
        List<HapPair> l = new ArrayList<>(size);
        byte[] a1 = new byte[size];
        byte[] a2 = new byte[size];
        Markers markers = shp.markers().restrict(start, end);
        for (int sample=0, nSamples=shp.nSamples(); sample<nSamples; ++sample) {
            for (int j=0; j<size; ++j) {
                a1[j] = shp.allele1(start+j, sample);
                a2[j] = shp.allele2(start+j, sample);
            }
            int idIndex = shp.idIndex(sample);
            l.add(new BitHapPair(markers, idIndex, a1, a2));
        }
        return new BasicSampleHapPairs(shp.samples(), l);
    }

    private static boolean[] switchHapOrder(SampleHapPairs prev,
            SampleHapPairs next, int nextStart, int overlap) {
        boolean[] switchHapOrder = new boolean[next.nHapPairs()];
        for (int hapPair=0, n=next.nHapPairs(); hapPair<n; ++hapPair) {
            boolean phaseUnknown = true;
            for (int k=1, m=2*overlap; k<=m && phaseUnknown; ++k) {
                int nextIndex = nextStart + ((k&1)==0 ? -k : k) / 2;
                if (nextIndex>=0 && nextIndex<overlap) {
                    int prevIndex = prev.nMarkers() - overlap + nextIndex;
                    byte a1 = next.allele1(nextIndex, hapPair);
                    byte a2 = next.allele2(nextIndex, hapPair);
                    if (a1!=a2) {
                        assert a1>=0 && a2>=0;
                        byte b1 = prev.allele1(prevIndex, hapPair);
                        byte b2 = prev.allele2(prevIndex, hapPair);
                        if ((a1==b1 && a2==b2)) {
                            switchHapOrder[hapPair] = false;
                            phaseUnknown = false;
                        }
                        else if (a1==b2 && a2==b1) {
                            switchHapOrder[hapPair] = true;
                            phaseUnknown = false;
                        }
                    }
                }
            }
        }
        return switchHapOrder;
    }
}
