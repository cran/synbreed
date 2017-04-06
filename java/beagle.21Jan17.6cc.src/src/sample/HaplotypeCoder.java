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
package sample;

import blbutil.ByteIndexArray;
import blbutil.ShiftedByteIndexArray;
import blbutil.CharIndexArray;
import blbutil.IntArray;
import blbutil.IntList;
import blbutil.WrappedIntArray;
import haplotype.SampleHapPairs;
import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * <p>Class {@code HaplotypeCoder} indexes the observed allele sequences
 * in reference and target haplotype pairs for a list of consecutive markers.
 * </p>
 * <p>Instances of class {@code HaplotypeCoder} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class HaplotypeCoder {

    private final int nRefHaps;
    private final int nHaps;
    private final SampleHapPairs refHapPairs;
    private final SampleHapPairs targetHapPairs;

    /**
     * Constructs a new {@code HaplotypeCoder} instance from the specified
     * data.
     * @param refHapPairs the reference haplotype pairs
     * @param targetHapPairs the target haplotype pairs
     * @throws IllegalArgumentException if
     * {@code refHapPairs.markers().equals(targetHapPairs.markers()) == false}
     * @throws NullPointerException if
     * {@code refHapPairs == null || targetHapPairs == null}
     */
    public HaplotypeCoder(SampleHapPairs refHapPairs,
            SampleHapPairs targetHapPairs) {
        if (refHapPairs.markers().equals(targetHapPairs.markers())==false) {
            throw new IllegalArgumentException("inconsistent markers");
        }
        this.nRefHaps = refHapPairs.nHaps();
        this.nHaps = nRefHaps + targetHapPairs.nHaps();
        this.refHapPairs = refHapPairs;
        this.targetHapPairs = targetHapPairs;
    }

    /**
     * Returns the reference haplotype pairs used to construct this.
     * @return the reference haplotype pairs used to construct this
     */
    public SampleHapPairs refHapPairs() {
        return refHapPairs;
    }

    /**
     * Returns the target haplotype pairs used to construct this.
     * @return the target haplotype pairs used to construct this
     */
    public SampleHapPairs targetHapPairs() {
        return targetHapPairs;
    }

    /**
     * Returns a two element array whose first element maps each reference
     * haplotype index to the index of the allele sequence carried by that
     * reference haplotype, and whose second element maps each target haplotype
     * index to the index of the allele sequence carried by that target
     * haplotype. The size of the first element of the returned array is
     * {@code this.refHapPairs().nHaps()}, and the size of the second
     * element of the returned array is {@code this.targetHapPairs().nHaps()}
     *
     * @param start the first marker index (inclusive)
     * @param end the last marker index (exclusive)
     * @return the haplotype indices for the reference and target samples.
     * @throws IllegalArgumentException if {@code start > end}
     * @throws IndexOutOfBoundsException if
     * {@code start < 0 || end >= this.refHapPairs.nMarkers()}
     */
    public IntArray[] run(int start, int end) {
        if (start >= end) {
            throw new IllegalArgumentException("start > end");
        }
        IntArray[] val = new IntArray[2];
        int[] haps = IntStream.range(0, nHaps).toArray();
        int[] alleles = new int[nHaps];
        IntList lastEnds = new IntList(1);
        lastEnds.add(nHaps);
        for (int m=start; m<end; ++m) {
            lastEnds = partition(m, alleles, haps, lastEnds);
        }
        setAllelesToHapIndices(alleles, lastEnds, haps);
        int nAlleles = lastEnds.size();
        if (nAlleles <= 128) {
            val[0] = new ByteIndexArray(alleles, 0, nRefHaps);
            val[1] = new ByteIndexArray(alleles, nRefHaps, nHaps);
        }
        else if (nAlleles <= 256) {
            val[0] = new ShiftedByteIndexArray(alleles, 0, nRefHaps);
            val[1] = new ShiftedByteIndexArray(alleles, nRefHaps, nHaps);
        }
        else if (nAlleles <= 65535) {
            val[0] = new CharIndexArray(alleles, 0, nRefHaps);
            val[1] = new CharIndexArray(alleles, nRefHaps, nHaps);
        }
        else {
            val[0] = new WrappedIntArray(Arrays.copyOfRange(alleles, 0, nRefHaps));
            val[1] = new WrappedIntArray(Arrays.copyOfRange(alleles, nRefHaps, nHaps));
        }
        return val;
    }

    private IntList partition(int marker, int[] alleles, int[] haps,
            IntList lastEnds) {
        IntList nextEnds = new IntList( (4*lastEnds.size())/3 + 1 );
        setAlleles(marker, alleles);
        int nAlleles = refHapPairs.marker(marker).nAlleles();
        int lastAllele = nAlleles - 1;

        int start = 0;
        for (int j=0, n=lastEnds.size(); j<n; ++j) {
            int end = lastEnds.get(j);
            for (int al=0; al<lastAllele; ++al) {
                int nextStart = partition(alleles, haps, start, end, al);
                if (nextStart > start) {
                    nextEnds.add(nextStart);
                    start = nextStart;
                }
            }
            if (end > start) {
                nextEnds.add(end);
            }
            start = end;
        }
        return nextEnds;
    }

    private void setAlleles(int marker, int[] alleles) {
        for (int j=0; j<nRefHaps; ++j) {
            alleles[j] = refHapPairs.allele(marker, j);
        }
        for (int j = nRefHaps; j<nHaps; ++j) {
            alleles[j] = targetHapPairs.allele(marker, j - nRefHaps);
        }
    }

    /* Returns the start index of second partitioned set */
    private int partition(int[] alleles, int[] haps, int start, int end,
            int splitAllele) {
        int nextStart = end;
        while (start < nextStart) {
            int allele = alleles[haps[start]];
            if (allele == splitAllele) {
                ++start;
            }
            else {
                --nextStart;
                int tmp = haps[nextStart];
                haps[nextStart] = haps[start];
                haps[start] = tmp;
            }
        }
        return nextStart;
    }

    private static void setAllelesToHapIndices(int[] alleles, IntList ends,
            int[] haps) {
        int start = 0;
        for (int j=0, n=ends.size(); j<n; ++j) {
            int end = ends.get(j);
            for (int k=start; k<end; ++k) {
                alleles[haps[k]] = j;
            }
            start = end;
        }
    }
}
