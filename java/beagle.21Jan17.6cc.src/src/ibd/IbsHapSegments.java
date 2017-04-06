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
package ibd;

import blbutil.IndexMap;
import blbutil.IntList;
import blbutil.IntPair;
import haplotype.HapPairs;
import haplotype.SampleHapPairs;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

/**
 * <p>Class {@code IbsHapSegments} identifies IBS haplotype segments in
 * a list of sample halotype pairs.
 * </p>
 * <p>Instances of {@code IbsHapSegments} are immutable.
 * </p>
 *
 * Reference: Gusev A, Lowe JK, Stoffel M, Daly MJ, Altshuler D, Breslow JL,
 *            Friedman JM, Pe'er I (2008) Whole population, genomewide mapping
 *            of hidden relatedness.  Genome Research 2009;19(2):318-26.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class IbsHapSegments {

    private static final int INIT_LIST_SIZE=500;

    private final SampleHapPairs haps;
    private final double[] pos;
    private final double minLength;
    private final int[] windowStarts;
    private final int[][][] idSets;

    private final int DEL = -67;

    /**
     * Constructs a new {@code IbsHapSegments} object from the specified data.
     * @param haps the sample haplotype pairs
     * @param pos an array of non-decreasing marker positions whose {@code j}-th
     * element is the position of marker {@code haps.marker(j)}
     * @param minLength the minimum length of a reported IBS segment
     *
     * @throws IllegalArgumentException if
     * {@code haps.nMarkers() != pos.length}
     * @throws IllegalArgumentException if {@code pos[j] < pos[j-1]}
     * for any {@code j} satisfing {@code (0 < j && j < pos.length)}
     * @throws IllegalArgumentException if
     * {@code (Double.isNaN(pos[j])==true || Double.isInfinite(pos[j]) == true)}
     * for any {@code j} satisfying {@code  (0 <= j && j < pos.length)}
     * @throws IllegalArgumentException if {@code  minLength <= 0.0f}
     * @throws NullPointerException if {@code haps == null || pos == null}
     */
    public IbsHapSegments(SampleHapPairs haps, double[] pos, double minLength) {
        checkArguments(haps, pos, minLength);
        this.haps = haps;
        this.pos = pos.clone();
        this.minLength = minLength;
        this.windowStarts = windowStarts(pos, minLength);
        this.idSets = idSets(haps, windowStarts);
    }

    /**
     * Constructs a new {@code IbsHapSegments} object with marker positions
     * defined to be marker indices.
     * @param haps the sample haplotype pairs
     * @param minMarkers the minimum number of shared markers in a reported
     * IBS segment
     * @throws NullPointerException if {@code haps == null}
     */
    public IbsHapSegments(SampleHapPairs haps, int minMarkers) {
        this(haps, pos(haps.nMarkers()), minMarkers);
    }

   private static double[] pos(int nMarkers) {
       double[] pos = new double[nMarkers];
       for (int j=0; j<pos.length; ++j) {
           pos[j] = j;
       }
       return pos;
   }

    private static void checkArguments(SampleHapPairs haps, double[] pos,
            double minLength) {
        if (minLength <= 0.0f) {
            throw new IllegalArgumentException("minLength: " + minLength);
        }
        if (haps.nMarkers()!= pos.length) {
            throw new IllegalArgumentException("haps.nMarkers()!= pos.length");
        }
        if (Double.isNaN(pos[0]) || Double.isInfinite(pos[0]) ) {
            throw new IllegalArgumentException(String.valueOf(pos[0]));
        }
        for (int j=1; j<pos.length; ++j) {
            if (Double.isNaN(pos[j]) || Double.isInfinite(pos[j]) ) {
                throw new IllegalArgumentException(String.valueOf(pos[j]));
            }
            if (pos[j] < pos[j-1]) {
                String s = "positions are not non-decreasing";
                throw new IllegalArgumentException(s);
            }
        }
    }

    private static int[] windowStarts(double[] pos, double minIbsLength) {
        double step = minIbsLength/2.0f;
        IntList indices = new IntList(pos.length/10);
        int index = 0;
        do {
            indices.add(index);
            double nextPos = pos[index] + step;
            index = nextIndex(pos, index, nextPos);
        } while (index < pos.length);
        return indices.toArray();
    }

    private static int nextIndex(double[] pos, int start, double targetPos) {
        int nextIndex = Arrays.binarySearch(pos, start, pos.length, targetPos);
        return (nextIndex<0) ? -nextIndex-1 : nextIndex;
    }

    private static int[][][] idSets(SampleHapPairs haps, int[] starts) {
        return IntStream.range(0, starts.length)
                .parallel()
                .mapToObj(j -> intPair(j, starts, haps.nMarkers()))
                .map(ip -> hapDictionary(haps, ip))
                .map(m -> fillIdSet(m, haps.nHaps()))
                .toArray(int[][][]::new);
    }

    private static IntPair intPair(int index, int[] starts, int nMarkers) {
        if (index+1 < starts.length) {
            return new IntPair(starts[index], starts[index+1]);
        }
        else {
            return new IntPair(starts[index], nMarkers);
        }
    }

    private static Map<Haplotype, IntList> hapDictionary(SampleHapPairs haps,
            IntPair ip) {
        Map<Haplotype, IntList> map = new HashMap<>();
        IntStream.range(0, haps.nHaps())  // forEach does not allow parallelization
                .mapToObj(h -> new Haplotype(haps, h, ip.first(), ip.second()))
                .forEach((h) -> {
                    IntList list = map.get(h);
                    if (list==null) {
                        list = new IntList(10);
                        map.put(h, list);
                    }
                    list.add(h.hapIndex());
                });
        return map;
    }

    private static int[][] fillIdSet(Map<Haplotype, IntList> hapMap, int nHaps) {
        int[][] value = new int[nHaps][];
        for (Haplotype key : hapMap.keySet()) {
            int[] ia = hapMap.get(key).toArray();
            for (int i : ia) {
                value[i] = ia;
            }
        }
        return value;
    }

    /**
     * Returns the sample haplotype pairs.
     * @return the sample haplotype pairs
     */
    public SampleHapPairs haps() {
        return haps;
    }

    /**
     * Returns an array of non-decreasing marker positions whose {@code j}-th
     * element is the position of marker {@code this.haps().marker(j)}.
     * @return an array of marker positions
     */
    public double[] pos() {
        return pos.clone();
    }

    /**
     * Returns the minimum length of an IBS segment.
     * @return the minimum length of an IBS segment
     */
    public double minIbsLength() {
        return minLength;
    }

    /**
     * Returns the list of haplotype segments for other haplotypes that
     * are IBS with the specified haplotype and have length greater
     * than or equal to {@code this.minIbsLength()}.
     *
     * @param hap the haplotype index
     * @return a list of IBS haplotype segments
     *
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap >= this.haps().nHaps()}
     */
    public List<HapSegment> find(int hap) {
        List<HapSegment> segments = new ArrayList<>(INIT_LIST_SIZE);
        int nil = Integer.MIN_VALUE;
        IndexMap prev = new IndexMap(haps.nHaps()-1, nil);
        IndexMap next = new IndexMap(haps.nHaps()-1, nil);
        int window = 0;
        matches(idSets, hap, window, prev);
        while (++window < idSets.length) {
            matches(idSets, hap, window, next);
            extend(prev, next);
            save(haps, hap, prev, windowStarts[window], segments);
            prev.clear();
            IndexMap tmp = prev;
            prev = next;
            next = tmp;
        }
        save(haps, hap, prev, haps.nMarkers(), segments);
        return segments;
    }

    /**
     * Returns a list of haplotype segments for other haplotypes
     * that are IBS with the specified haplotype and that have length greater
     * than or equal to {@code this.minIbsLength()}. An IBS segment is
     * permitted (but not required) to be excluded from the returned
     * list if both end-points of the IBD segment are interior points of
     * another IBD segment.
     *
     * @param hap the haplotype index
     * @return a list of IBS haplotype segments
     *
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap >= this.haps().nHaps()}
     */
    public List<HapSegment> filteredFind(int hap) {
        List<HapSegment> segments = new ArrayList<>(INIT_LIST_SIZE);
        int nil = Integer.MIN_VALUE;
        IndexMap prev = new IndexMap(haps.nHaps()-1, nil);
        IndexMap next = new IndexMap(haps.nHaps()-1, nil);
        int window = 0;
        matches(idSets, hap, window, prev);
        while (++window < idSets.length) {
            matches(idSets, hap, window, next);
            int minExtendedStartWindow = extend(prev, next);
            filteredSave(haps, hap, prev, minExtendedStartWindow,
                    windowStarts[window], segments);
            prev.clear();
            IndexMap tmp = prev;
            prev = next;
            next = tmp;
        }
        filteredSave(haps, hap, prev, window, haps.nMarkers(), segments);
        return segments;
    }

    private void matches(int[][][] idSets, int hap, int window, IndexMap map) {
        assert map.size()==0;
        int[] hapIbsSet = idSets[window][hap];
        for (int h : hapIbsSet) {
            if (h!=hap) {
                map.put(h, window);
            }
        }
    }

    /* Returns minimum start window index from extended segments */
    private int extend(IndexMap prev, IndexMap next) {
        int nil = next.nil();
        int minStart = Integer.MAX_VALUE;
        for (int i=0, n=next.size(); i<n; ++i) {
            int hap = next.enumeratedKey(i);
            int prevStart = prev.get(hap);
            if (prevStart != nil) {
                next.put(hap, prevStart);
                prev.put(hap, DEL);
                if (prevStart < minStart) {
                    minStart = prevStart;
                }
            }
        }
        return minStart;
    }

    private void save(HapPairs haps, int hap1,
            IndexMap prev, int prevExclEnd, List<HapSegment> segments) {
        for (int i=0, n=prev.size(); i<n; ++i) {
            int hap2 = prev.enumeratedKey(i);
            int startWindow = prev.enumeratedValue(i);
            if (startWindow != DEL) {
                int start = start(haps, hap1, hap2, windowStarts[startWindow]);
                int inclEnd = inclusiveEnd(haps, hap1, hap2, prevExclEnd);
                if ( (pos[inclEnd] - pos[start]) >= minLength) {
                    segments.add( new HapSegment(hap2, start, inclEnd) );
                }
            }
        }
    }

    private void filteredSave(HapPairs haps, int hap1, IndexMap prev,
            int minExtendedStartWindow, int prevExclEnd, List<HapSegment> segments) {
        for (int i=0, n=prev.size(); i<n; ++i) {
            int hap2 = prev.enumeratedKey(i);
            int startWindow = prev.enumeratedValue(i);
            if (startWindow != DEL && startWindow <= minExtendedStartWindow) {
                int start = start(haps, hap1, hap2, windowStarts[startWindow]);
                int inclEnd = inclusiveEnd(haps, hap1, hap2, prevExclEnd);
                if ( (pos[inclEnd] - pos[start]) >= minLength) {
                    segments.add( new HapSegment(hap2, start, inclEnd) );
                }
            }
        }
    }

    private int start(HapPairs haps, int hap1, int hap2, int start) {
        while (start>0
                && haps.allele(start-1, hap1)==haps.allele(start-1, hap2)) {
            --start;
        }
        return start;
    }

    private int inclusiveEnd(HapPairs haps, int hap1, int hap2, int end) {
        while (end<haps.nMarkers()
                && haps.allele(end, hap1)==haps.allele(end, hap2)) {
            ++end;
        }
        return end-1;
    }
}
