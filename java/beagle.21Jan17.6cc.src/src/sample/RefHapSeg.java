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

import blbutil.IntArray;
import blbutil.IntList;
import haplotype.SampleHapPairs;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * <p>Class {@code RefHapSeg} represents a chromosome segment of
 * reference haplotypes.
 * </p>
 * <p>Instances of class {@code RefHapSeg} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class RefHapSeg {

    private final SampleHapPairs refHapPairs;
    private final int start;    // inclusive
    private final int end;      // exclusive
    private final IntArray hapToSeq;
    private final int[] seqToHap;

    /**
     * Constructs a new {@code RefHapSegs} instance from the specified data.
     * @param refHapPairs the reference haplotype pairs
     * @param start the starting marker index (inclusive)
     * @param end the ending marker index (exclusive)
     * @throws IllegalArgumentException if
     * {@code start < 0 || start >= end || end > refHapPairs.nMarkers()}
     * @throws NullPointerException if {@code refHapPairs == null}
     */
    public RefHapSeg(SampleHapPairs refHapPairs, int start, int end) {
        if (start < 0 || start >= end || end > refHapPairs.nMarkers()) {
            throw new IllegalArgumentException();
        }
        HapSegData hapSegData = new HapSegData(refHapPairs, start, end);
        this.refHapPairs = refHapPairs;
        this.start = start;
        this.end = end;
        this.hapToSeq = hapSegData.hap2Seq();
        this.seqToHap = hapSegData.seq2Hap();
    }

    /**
     * Returns the reference haplotype pairs.
     * @return the reference haplotype pairs
     */
    public SampleHapPairs refHapPairs() {
        return refHapPairs;
    }

    /**
     * Return the number of reference allele sequences in this segment.
     * @return the number of reference allele sequences in this segment
     */
    public int nSeq() {
        return seqToHap.length;
    }

    /**
     * Return the index of the reference allele sequence in this segment
     * for the specified reference haplotype.
     * @param hap a haplotype index
     * @return the index of the reference allele sequence in this segment
     * for the specified reference haplotype
     *
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap >= this.refHapPairs().nHaps()}
     */
    public int seq(int hap) {
        return hapToSeq.get(hap);
    }

    /**
     * Return the specified reference haplotype allele.
     * @param marker index of a marker in this segment
     * @param seq index of a reference allele sequence in this segment
     * @return the specified reference haplotype allele
     *
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= (this.end() - this.start())}
     * @throws IndexOutOfBoundsException if
     * {@code seq < 0 || seg >= this.nSeq(segment)}
     */
    public int allele(int marker, int seq) {
        int refIndex = start + marker;
        if (marker < 0 || refIndex >= end) {
            throw new IllegalArgumentException(String.valueOf(marker));
        }
        return refHapPairs.allele(refIndex, seqToHap[seq]);
    }

    /**
     * Returns the starting marker index (inclusive) of this segment.
     * @return the starting marker index (inclusive) of this segment
     */
    public int start() {
        return start;
    }

    /**
     * Returns the ending marker index (exclusive) of this segment.
     * @return the ending marker index (exclusive) of this segment
     */
    public int end() {
        return end;
    }

    private static class HapSegData {
        private final SampleHapPairs refHapPairs;
        private final int[] hap2seq;
        private final IntList seq2Cnt;

        private final List<IntList> seq2AlleleMap;
        private final IntList seq2NonMajorCnt;
        private final IntList nonMajorSeq;

        public HapSegData(SampleHapPairs refHapPairs, int start, int end) {
            this.refHapPairs = refHapPairs;
            this.hap2seq = new int[refHapPairs.nHaps()];
            this.seq2Cnt = new IntList(200);
            this.seq2Cnt.add(hap2seq.length);
            this.seq2AlleleMap = new ArrayList<>(200);
            this.seq2AlleleMap.add(new IntList(4));
            this.seq2NonMajorCnt = new IntList(200);
            this.nonMajorSeq =  new IntList(20);
            for (int m=start; m<end; ++m) {
                if (refHapPairs.storesNonMajorIndices(m)) {
                    lowMafUpdate(m);
                }
                else {
                    highMafUpdate(m);
                }
            }
        }

        public IntArray hap2Seq() {
            return IntArray.create(hap2seq, 0, seq2AlleleMap.size()-1);
        }

        public int[] seq2Hap() {
            int[] seqToHap = new int[seq2AlleleMap.size()];
            Arrays.fill(seqToHap, -1);
            for (int h=0; h<hap2seq.length; ++h) {
                int seq = hap2seq[h];
                if (seqToHap[seq] == -1) {
                    seqToHap[seq] = h;
                }
            }
            return seqToHap;
        }

        private void lowMafUpdate(int marker) {
            setSeqToAlleleMap(marker); // major allele not in any allele map
            int nAlleles = refHapPairs.nAlleles(marker);
            int majorAllele = refHapPairs.majorAllele(marker);
            for (int al=0; al<nAlleles; ++al) {
                if (al!=majorAllele) {
                    int nCopies = refHapPairs.alleleCount(marker, al);
                    for (int c=0; c<nCopies; ++c) {
                        int h = refHapPairs.hapIndex(marker, al, c);
                        int seq = hap2seq[h];
                        IntList list = seq2AlleleMap.get(seq);
                        int index = indexOfAllele(list, al);
                        if (index < list.size()) {
                            updateHap2Seq(h, list.get(index+1));
                        }
                    }
                }
            }
        }

        private void setSeqToAlleleMap(int marker) {
            int nAlleles = refHapPairs.nAlleles(marker);
            int majorAllele = refHapPairs.majorAllele(marker);

            seq2NonMajorCnt.clear();
            for (int j=0, n=seq2Cnt.size(); j<n; ++j) {
                seq2AlleleMap.get(j).clear();
                seq2NonMajorCnt.add(0);
            }

            for (int al=0; al<nAlleles; ++al) {
                if (al!=majorAllele) {
                    updateNonMajorSeq(marker, al, seq2NonMajorCnt, nonMajorSeq);
                    updateAlleleToSeqList(al, seq2NonMajorCnt, nonMajorSeq);
                    for (int j=0; j<nonMajorSeq.size(); ++j) {
                        seq2NonMajorCnt.set(nonMajorSeq.get(j), 0);
                    }
                    nonMajorSeq.clear();
                }
            }
        }

        private void updateNonMajorSeq(int marker, int allele,
                IntList nonMajorSeqCnts, IntList nonMajorSeq) {
            int nCopies = refHapPairs.alleleCount(marker, allele);
            for (int c=0; c<nCopies; ++c) {
                int h = refHapPairs.hapIndex(marker, allele, c);
                int seq = hap2seq[h];
                if (nonMajorSeqCnts.getAndIncrement(seq) == 0) {
                    nonMajorSeq.add(seq);
                }
            }
        }

        private void updateAlleleToSeqList(int allele, IntList nonMajorSeqCnts,
                IntList nonMajorSeq) {
            for (int j=0; j<nonMajorSeq.size(); ++j) {
                int seq = nonMajorSeq.get(j);
                if (nonMajorSeqCnts.get(seq) < seq2Cnt.get(seq)) {
                    IntList list = seq2AlleleMap.get(seq);
                    mapAlleleToNewSeq(list, allele);
                }
            }
        }

        private void highMafUpdate(int marker) {
            for (int j=0, n=seq2AlleleMap.size(); j<n; ++j) {
                seq2AlleleMap.get(j).clear();
            }
            for (int h=0; h<hap2seq.length; ++h) {
                int seq = hap2seq[h];
                int allele = refHapPairs.allele(marker, h);
                IntList list = seq2AlleleMap.get(seq);
                if (list.isEmpty()) {
                    list.add(allele);
                    list.add(seq);
                }
                else {
                    int index = indexOfAllele(list, allele);
                    if (index==list.size()) {
                        mapAlleleToNewSeq(list, allele);
                    }
                    updateHap2Seq(h, list.get(index+1));
                }
            }
        }

        private int indexOfAllele(IntList list, int allele) {
            int index=0;
            while (index < list.size() && list.get(index)!=allele) {
                index+=2;
            }
            return index;
        }

        private void mapAlleleToNewSeq(IntList list, int allele) {
            list.add(allele);
            list.add(seq2AlleleMap.size());
            seq2AlleleMap.add(new IntList(4));
            seq2Cnt.add(0);
        }

        private void updateHap2Seq(int h, int seq) {
            seq2Cnt.decrementAndGet(hap2seq[h]);
            hap2seq[h] = seq;
            seq2Cnt.incrementAndGet(hap2seq[h]);
        }
    }
}
