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

import beagleutil.Phase;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import vcf.Markers;
import vcf.VcfRecord;

/**
 * Class {@code ConsensusPhasing} determines a consensus phasing from
 * multiple estimated haplotype pairs for an individual.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class ConsensusPhasing {

    private ConsensusPhasing() {
        // private constructor prevents instantiation
    }

    /**
     * Returns a list of consensus haplotype pairs (one pair per individual)
     * for the individuals with haplotype pairs in the the specified list.
     *
     * @param haps a list of haplotype pairs.  The list may contain multiple
     * haplotype pairs for each individual.
     * @return a list of consensus haplotype pairs.
     *
     * @throws IllegalArgumentException if there exist
     * {@code 0 <= j < k < haps.size()} such that
     * haps.get(j).markers().equals(haps.get(k).markers()==false}
     * @throws NullPointerException if {@code haps==null}.
     */
    public static List<HapPair> consensusHaps(List<HapPair> haps) {
        List<HapPair> copy = new ArrayList<>(haps);
        if (copy.isEmpty()) {
            return copy;
        }
        checkMarkers(copy);
        Random random = new Random(copy.size());

        Collections.sort(copy, BasicSampleHapPairs.hapsComparator());
        List<HapPair> consensus = new ArrayList<>(copy.size()/20);
        int start = 0;
        while (start < copy.size()) {
            int end = start+1;
            while (end < copy.size()
                    && copy.get(end).idIndex()==copy.get(start).idIndex()) {
                ++end;
            }
            if (end-start==1) {
                consensus.add(copy.get(start));
            }
            else {
                consensus.add(consensus(copy.subList(start, end), random));
            }
            start = end;
        }
        return consensus;
    }

    private static void checkMarkers(List<HapPair> hapList) {
        Markers markers = hapList.get(0).markers();
        for (int j=1; j<hapList.size(); ++j) {
            if (markers.equals(hapList.get(j).markers())==false) {
                throw new IllegalArgumentException("inconsistent markers");
            }
        }
    }

    private static HapPair consensus(List<HapPair> hapList, Random rand) {
        int idIndex = hapList.get(0).idIndex();
        Markers markers = hapList.get(0).markers();
        int nMarkers = markers.nMarkers();
        Phase[] lastPhaseArray = new Phase[hapList.size()];
        Phase[] phaseArray = new Phase[hapList.size()];
        byte[] alleles1 = new byte[nMarkers];
        byte[] alleles2 = new byte[nMarkers];

        Phase phase = null;
        for (int m=0; m<nMarkers; ++m) {
            int consensusGT = consensusGT(hapList, markers, m, rand);
            int sumGt = 0;
            int a1 = consensusGT;
            int a2 = 0;
            while (a1>a2) {
                ++a2;
                sumGt += a2;
                a1 = consensusGT - sumGt;
            }
            if (a1!=a2) {
                assert a1<a2;
                setPhaseArray(hapList, m, a1, a2, phaseArray);
                if (phase==null) {
                    phase = Phase.IDENTICAL;
                }
                else {
                    Phase relPhase = relPhase(lastPhaseArray, phaseArray, rand);
                    phase = (phase==relPhase) ? Phase.IDENTICAL : Phase.OPPOSITE;
                    if (phase==Phase.OPPOSITE) {
                        int tmp = a1;
                        a1 = a2;
                        a2 = tmp;
                    }
                }
                Phase[] tmp = phaseArray;
                phaseArray = lastPhaseArray;
                lastPhaseArray = tmp;
            }
            alleles1[m] = (byte) a1;
            alleles2[m] = (byte) a2;
        }
        return new BitHapPair(markers, idIndex, alleles1, alleles2);
    }

    private static int consensusGT(List<HapPair> hapList, Markers markers,
            int marker, Random random) {
        int[] gtCounts  = gtCounts(hapList, markers, marker);
        int pivot = random.nextInt(gtCounts.length);
        int bestGt = pivot;
        for (int gt=pivot+1; gt<gtCounts.length; ++gt) {
            if (gtCounts[gt] > gtCounts[bestGt]) {
                bestGt = gt;
            }
        }
        for (int gt=0; gt<pivot; ++gt) {
            if (gtCounts[gt] > gtCounts[bestGt]) {
                bestGt = gt;
            }
        }
        return bestGt;
    }

    private static int[] gtCounts(List<HapPair> hapList, Markers markers,
            int marker) {
        int nGt = markers.marker(marker).nGenotypes();
        int[] gtCounts = new int[nGt];
        for (int j=0, n=hapList.size(); j<n; ++j) {
            HapPair hp = hapList.get(j);
            byte a1 = hp.allele1(marker);
            byte a2 = hp.allele2(marker);
            int gt = VcfRecord.gtIndex(a1, a2);
            ++gtCounts[gt];
        }
        return gtCounts;
    }

    private static void setPhaseArray(List<HapPair> hapList, int marker,
            int a1, int a2, Phase[] phaseArray) {
        assert a1<a2;
        assert phaseArray.length == hapList.size();
        for (int j=0; j<phaseArray.length; ++j) {
            byte b1 = hapList.get(j).allele1(marker);
            byte b2 = hapList.get(j).allele2(marker);
            if (a1==b1 && a2==b2) {
                phaseArray[j] =Phase.IDENTICAL;
            }
            else if (a1==b2 && a2==b1) {
                phaseArray[j] = Phase.OPPOSITE;
            }
            else {
                phaseArray[j] = Phase.INCONSISTENT;
            }
        }
    }

    private static Phase relPhase(Phase[] rp1, Phase[] rp2, Random rand) {
        assert rp1.length == rp2.length;
        int identCnt = 0;
        int oppCnt = 0;
        for (int j=0; j<rp1.length; ++j) {
            if (rp1[j]==Phase.IDENTICAL) {
                switch (rp2[j]) {
                    case IDENTICAL: ++identCnt; break;
                    case OPPOSITE: ++oppCnt; break;
                }
            }
            else if(rp1[j] == Phase.OPPOSITE) {
                switch (rp2[j]) {
                    case IDENTICAL: ++oppCnt; break;
                    case OPPOSITE: ++identCnt; break;
                }
            }
        }
        if (identCnt > oppCnt) {
            return Phase.IDENTICAL;
        }
        else if (oppCnt > identCnt) {
            return Phase.OPPOSITE;
        }
        else {
            return rand.nextBoolean() ? Phase.IDENTICAL : Phase.OPPOSITE;
        }
    }
}
