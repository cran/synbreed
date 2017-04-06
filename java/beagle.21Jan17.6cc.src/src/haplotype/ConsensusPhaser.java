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
import beagleutil.Samples;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;
import vcf.Markers;
import vcf.VcfRecord;

/**
 * Class {@code ConsensusPhaser} contains a static method for
 * calculating a consensus phasing from multiple estimated haplotype pairs
 * for an individual.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class ConsensusPhaser {

    private ConsensusPhaser() {
        // private constructor prevents instantiation
    }

    /**
     * Returns a list of consensus haplotype pairs (one pair per individual)
     * sorted in order of increasing sample index. The specified list of
     * haplotype pairs may contain multiple haplotype pairs for each individual.
     *
     * @param hapPairs a list of haplotype pairs
     * @return a list of consensus haplotype pairs
     *
     * @throws IllegalArgumentException if
     * {@code (hapPairs.get(j).markers().equals(hapPairs.get(k).markers() == false)}
     * for any {@code j, k} satisfying {@code 0 <= j < k < hapPairs.size()}
     * @throws IllegalArgumentException if
     * {@code (hapPairs.get(j).samples().equals(hapPairs.get(k).samples() == false)}
     * for any {@code j, k} satisfying {@code 0 <= j < k < hapPairs.size()}
     * @throws NullPointerException if {@code hapPairs == null}
     */
    public static List<HapPair> run(List<HapPair> hapPairs) {
        List<HapPair> copy = new ArrayList<>(hapPairs);
        if (copy.isEmpty()) {
            return copy;
        }
        checkMarkers(copy);
        Random random = new Random(copy.size());

        Collections.sort(copy, hapsComparator());
        List<HapPair> consensus = new ArrayList<>(copy.size()/20);
        int start = 0;
        while (start < copy.size()) {
            int end = start+1;
            while (end < copy.size()
                    && copy.get(end).sampleIndex()==copy.get(start).sampleIndex()) {
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
        Samples samples = hapList.get(0).samples();
        for (int j=1; j<hapList.size(); ++j) {
            if (markers.equals(hapList.get(j).markers())==false) {
                throw new IllegalArgumentException("inconsistent markers");
            }
            if (samples.equals(hapList.get(j).samples())==false) {
                throw new IllegalArgumentException("inconsistent samples");
            }
        }
    }

    private static HapPair consensus(List<HapPair> hapList, Random rand) {
        HapPair firstHP = hapList.get(0);
        int sampleIndex = firstHP.sampleIndex();
        Samples samples = firstHP.samples();
        Markers markers = firstHP.markers();
        int nMarkers = markers.nMarkers();
        Phase lastConsensus = null;
        Phase[] lastPhase = new Phase[hapList.size()];
        Phase[] currentPhase = new Phase[hapList.size()];
        int[] alleles1 = new int[nMarkers];
        int[] alleles2 = new int[nMarkers];

        for (int m=0; m<nMarkers; ++m) {
            int hp = hapPairWithConsensusGT(hapList, markers, m, rand);
            // retrieve actual allele order to match input phased data
            int a1 = hapList.get(hp).allele1(m);
            int a2 = hapList.get(hp).allele2(m);
            if (a1!=a2) {
                storePhase(hapList, m, a1, a2, currentPhase);
                Phase consensus;
                if (lastConsensus != null) {
                    Phase relPhase = relPhase(lastPhase, currentPhase, rand);
                    if (relPhase == Phase.IDENTICAL) {
                        consensus = lastConsensus;
                    }
                    else {
                        assert relPhase == Phase.OPPOSITE;
                        consensus = flip(lastConsensus);
                    }
                    if ( (consensus == Phase.IDENTICAL && a1 > a2)
                            || (consensus == Phase.OPPOSITE && a1 < a2)) {
                        int tmp = a1;
                        a1 = a2;
                        a2 = tmp;
                    }
                }
                lastConsensus = a1 < a2 ? Phase.IDENTICAL : Phase.OPPOSITE;
                Phase[] tmp = currentPhase;
                currentPhase = lastPhase;
                lastPhase = tmp;
            }
            alleles1[m] = a1;
            alleles2[m] = a2;
        }
        return new BitHapPair(markers, samples, sampleIndex, alleles1, alleles2);
    }

    private static Phase flip(Phase phase) {
        if (phase == Phase.IDENTICAL) {
            return Phase.OPPOSITE;
        }
        else if (phase == Phase.OPPOSITE) {
            return Phase.IDENTICAL;
        }
        else {
            throw new IllegalArgumentException(phase.toString());
        }
    }

    private static int hapPairWithConsensusGT(List<HapPair> hapList,
            Markers markers, int marker, Random random) {
        int consensusGT = consensusGT(hapList, markers, marker, random);
        for (int j=0, n=hapList.size(); j<n; ++j) {
            HapPair hp = hapList.get(j);
            int a1 = hp.allele1(marker);
            int a2 = hp.allele2(marker);
            if (VcfRecord.gtIndex(a1, a2) == consensusGT) {
                return j;
            }
        }
        assert false;
        throw new IllegalArgumentException("no sample with consensus GT");
    }

    private static int consensusGT(List<HapPair> hapList, Markers markers,
            int marker, Random random) {
        int[] gtCounts  = gtCounts(hapList, markers, marker);
        int start = random.nextInt(gtCounts.length);
        int bestGt = start;
        for (int j=1; j<gtCounts.length; ++j) {
            int gt = start + j;
            if (gt >= gtCounts.length) {
                gt -= gtCounts.length;
            }
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
            int a1 = hp.allele1(marker);
            int a2 = hp.allele2(marker);
            int gt = VcfRecord.gtIndex(a1, a2);
            ++gtCounts[gt];
        }
        return gtCounts;
    }

    private static void storePhase(List<HapPair> hapList, int marker,
            int a1, int a2, Phase[] phaseArray) {
        assert phaseArray.length == hapList.size();
        for (int j=0; j<phaseArray.length; ++j) {
            int b1 = hapList.get(j).allele1(marker);
            int b2 = hapList.get(j).allele2(marker);
            if ( (a1==b1 && a2==b2) || (a1==b2 && a2==b1) ) {
                phaseArray[j] = (b1 < b2) ? Phase.IDENTICAL : Phase.OPPOSITE;
            }
            else {
                phaseArray[j] = Phase.INCONSISTENT;
            }
        }
    }

    private static Phase relPhase(Phase[] ph1, Phase[] ph2, Random rand) {
        assert ph1.length == ph2.length;
        int identCnt = 0;
        int oppCnt = 0;
        for (int j=0; j<ph1.length; ++j) {
            if (ph1[j]==Phase.IDENTICAL) {
                switch (ph2[j]) {
                    case IDENTICAL: ++identCnt; break;
                    case OPPOSITE: ++oppCnt; break;
                }
            }
            else if(ph1[j] == Phase.OPPOSITE) {
                switch (ph2[j]) {
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

    /**
     * Returns a {@code Comparator<HapPairInterface>}
     * whose {@code compare(hp1, hp2)} method returns -1, 0, or 1
     * depending on whether {@code hp1.idIndex()} is less than, equal,
     * or greater than {@code hp2.idIndex()}.
     * @return a {@code Comparator<HapPairInterface>}
     * whose {@code compare(hp1, hp2)} method compares two
     * haplotype pairs for order.
     */
    private static Comparator<HapPair> hapsComparator() {
        return (HapPair hp1, HapPair hp2) -> {
            int i1 = hp1.sampleIndex();
            int i2 = hp2.sampleIndex();
            if (i1==i2) {
                return 0;
            }
            else {
                return (i1 < i2) ? -1 : 1;
            }
        } ;
    }
}
