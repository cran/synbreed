/*
 * Copyright 2013 Brian L. Browning
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package sample;

import beagleutil.CenteredIntIntervalTree;
import beagleutil.IntIntervalTree;
import blbutil.IndexSet;
import dag.Dag;
import dag.MergeableDag;
import haplotype.HapPairs;
import haplotype.SampleHapPairs;
import ibd.HapSegment;
import ibd.IbsHapSegments;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * <p>Class {@code RestrictedDag} is a wrapper for a {@code Dag}
 * object that stores segments of identity by descent.
 * </p>
 * <p>Instances of class {@code RestrictedDag} are immutable.
 * </p>
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class RestrictedDag {

    private static final int END_FILTER = 1;

    private final SampleHapPairs haps;
    private final Dag dag;
    private final int[][] hapStates;
    private final IbsHapSegments hapSegments;
    private final double[] pos;
    private final double ibdExtend;

    /**
     * Constructs a {@code RestrictedDag} instance.
     * @param haps the sample haplotypes
     * @param weights an array of length {@code hapPairs.nHaps()}
     * whose {@code j}-th element is the weight for the {@code j}-th haplotype
     * @param nInitLevels the number of initial levels to read
     * @param scale a parameter that multiplicatively scales the node
     * similarity threshold
     * @param ibdLength the minimum length of an IBD segment
     * @param ibdExtend the length by which an IBD segment will be extended
     *
     * @throws IllegalArgumentException if {@code hapPairs.nMarkers() == 0}
     * @throws IllegalArgumentException if
     * {@code weights.length != 2*haps.nSamples()}
     * @throws IllegalArgumentException if
     * {@code (weights[j] <= 0 || Float.isFinite(weights[j]) == false)}
     * for any {@code j} satisfying {@code (0 <= j && j < weights.length)}
     * @throws IllegalArgumentException if {@code nInitLevels < 1}
     * @throws IllegalArgumentException if
     * {@code Double.isFinite(scale) == false || scale <= 0}
     * @throws IllegalArgumentException if
     * {@code ibdLength < 0 || ibdExtend < 0}
     * @throws NullPointerException if
     * {@code hapPairs == null || weights == null}
     */
    public RestrictedDag(SampleHapPairs haps, float[] weights, int nInitLevels,
            float scale, double ibdLength, double ibdExtend) {
        if (ibdLength <= 0d) {
            throw new IllegalArgumentException(String.valueOf(ibdLength));
        }
        if (ibdExtend <= 0d) {
            throw new IllegalArgumentException(String.valueOf(ibdExtend));
        }
        this.haps = haps;
        this.ibdExtend = ibdExtend;
        this.dag = MergeableDag.dag(haps, weights, scale, nInitLevels);
        this.pos = pos(dag);
        this.hapStates = hapStates(dag, haps);
        this.hapSegments = new IbsHapSegments(haps, pos, ibdLength);
    }

    private static double[] pos(Dag dag) {
        double[] pos = dag.posArray();
        double scaleFactor = 0.2;
        for (int j=0; j<pos.length; ++j) {
            pos[j] *= scaleFactor;
        }
        return pos;
    }

    /**
     * Returns a int[][] array whose (j,k)-th element is the edge state
     * at the j-th marker that is traversed by the k-th haplotype.
     * @param dag  the DAG constructed by the specified haplotypes
     * @param haps the haplotype pairs used to construct the specified DAG
     * @return a int[][] array whose (j,k)-th element is the edge state
     * at the j-th marker that is traversed by the k-th haplotype
     * @throws IllegalArgumentException if the specified haplotypes
     * pairs are inconsistent with the specified DAG
     * @throws NullPointerException if {@code dag == null || haps == null}
     */
    private static int[][] hapStates(Dag dag, HapPairs haps) {
        if (dag.nLevels()!=haps.nMarkers()) {
            throw new IllegalArgumentException("dag.nMarkers()!=haps.nMarkers()");
        }
        int nHaps = haps.nHaps();
        int[][] states = new int[dag.nLevels()][nHaps];
        for (int h=0; h<nHaps; ++h) {
            int node = 0;
            int symbol = haps.allele(0, h);
            states[0][h] = dag.outEdgeBySymbol(0, node, symbol);
            for (int j=1; j<states.length; ++j) {
                node = dag.childNode(j-1, states[j-1][h]);
                symbol = haps.allele(j, h);
                states[j][h] = dag.outEdgeBySymbol(j, node, symbol);
                assert states[j][h] != -1;
            }
        }
        return states;
    }

    /**
     * Returns the haplotypes used to construct {@code this}.
     * @return the haplotypes used to construct {@code this}
     */
    public SampleHapPairs sampleHaps() {
        return haps;
    }

    /**
     * Returns the DAG.
     * @return the DAG
     */
    public Dag dag() {
        return dag;
    }

    /**
     * Returns the permitted states for the specified sample.
     * @param sample the sample index
     * @return the permitted states for the specified sample
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= haps.nSamples()}
     */
    public DiploidStates singleStates(int sample) {
        if (sample < 0 || sample >= haps.nSamples()) {
            throw new IndexOutOfBoundsException("sample: " + sample);
        }
        int hap1 = 2*sample;
        int hap2 = 2*sample + 1;
        List<HapSegment> hapSegs1 = ibsSegs(hap1);
        List<HapSegment> hapSegs2 = ibsSegs(hap2);

        return new SinglePermittedStates(haps.nMarkers(), haps.nHaps(),
                sample, hapSegs1, hapSegs2);
    }

    private List<HapSegment> ibsSegs(int hap) {
        List<HapSegment> hapSegs = hapSegments.filteredFind(hap);
        containmentFilter(hapSegs, END_FILTER);
        return hapSegs;
    }

    /* filter if minimum requirements are not met */
    private void containmentFilter(List<HapSegment> ibdSegments,
            int minEndDiff) {
        assert minEndDiff >= 0;
        Collections.sort(ibdSegments, modStartComparator());
        if (ibdSegments.isEmpty()==false) {
            List<HapSegment> list = new LinkedList<>();
            List<HapSegment> filtered = new ArrayList<>(ibdSegments.size()/5);
            for (int k=0, m=ibdSegments.size(); k<m; ++k) {
                HapSegment hs = ibdSegments.get(k);
                boolean exclude = false;
                Iterator<HapSegment> it = list.iterator();
                while (it.hasNext() && exclude==false) {
                    HapSegment cover = it.next();
                    int cStart = cover.start();
                    int cEnd = cover.end();
                    if (cEnd <= hs.start() ) {
                        it.remove();
                    }
                    else {
                        if ( (hs.start() - cStart) >= minEndDiff
                                && (cEnd - hs.end()) >= minEndDiff) {
                            exclude = true;
                        }
                    }
                }
                if (exclude==false) {
                    list.add(hs);
                    filtered.add(hs);
                }
            }
            ibdSegments.clear();
            ibdSegments.addAll(filtered);
        }
    }


    // for a given start field, the end field is sorted in in reverse order.
    private static Comparator<HapSegment> modStartComparator() {
        return (HapSegment hs1, HapSegment hs2) -> {
            if (hs1.start() != hs2.start()) {
                return (hs1.start() < hs2.start()) ? -1 : 1;
            }
            else if (hs1.end() != hs2.end()) {
                return (hs1.end() > hs2.end()) ? -1 : 1;
            }
            if (hs1.hap() != hs2.hap()) {
                return (hs1.hap() < hs2.hap()) ? -1 : 1;
            }
            return 0;
        };
    }

    private int modifyStart(HapSegment targetHS,
            IntIntervalTree<HapSegment> tree) {
        int maxStart = extStartIndex(targetHS.start(), ibdExtend, pos);
        int minEnd = targetHS.end();
        List<HapSegment> list = new ArrayList<>(10);
        tree.intersectAll(maxStart, minEnd, list);
        return list.isEmpty() ? maxStart : targetHS.start();
    }

    private int modifyEnd(HapSegment targetHS,
            IntIntervalTree<HapSegment> tree) {
        int maxStart = targetHS.start();
        int minEnd = extEndIndex(targetHS.end(), ibdExtend, pos);
        List<HapSegment> list = new ArrayList<>(10);
        tree.intersectAll(maxStart, minEnd, list);
        return list.isEmpty() ? minEnd : targetHS.end();
    }

    private static int extStartIndex(int start, double extension, double[] pos) {
        double target = pos[start] - extension;
        int x = Arrays.binarySearch(pos, target);
        return (x<0) ? -x-1 : x;
    }

    private static int extEndIndex(int end, double extension, double[] pos) {
        double target = pos[end] + extension;
        int x = Arrays.binarySearch(pos, target);
        return (x<0) ? -x-2 : x; // end is inclusive
    }

    private class SinglePermittedStates implements DiploidStates {

        private final int nMarkers;
        private final IndexSet indices1;
        private final IndexSet indices2;
        private final IntIntervalTree<HapSegment> tree1;
        private final IntIntervalTree<HapSegment> tree2;

        private int marker = -1;
        private int i1 = 0;
        private int i2 = 0;
        private int edge1 = -1;
        private int edge2 = -1;
        private boolean rev = false;

        private SinglePermittedStates(int nMarkers, int nHaps, int sample,
                List<HapSegment> list1, List<HapSegment> list2) {
            int hap1 = 2*sample;
            int hap2 = 2*sample + 1;
            this.nMarkers = nMarkers;
            this.indices1 = new IndexSet(nHaps);
            this.indices2 = new IndexSet(nHaps);

            List<HapSegment> extList1 = extendSegment(hap1, list1);
            List<HapSegment> extList2 = extendSegment(hap2, list2);

            tree1 = getTree(nMarkers, extList1);
            tree2 = getTree(nMarkers, extList2);
        }

        private List<HapSegment> extendSegment(int hap, List<HapSegment> ibdSegs) {
            List<HapSegment> extendedSegs = new ArrayList<>(ibdSegs.size());
            IntIntervalTree<HapSegment> tree = getTree(haps.nMarkers(), ibdSegs);

            // permit states traversed by hap
            int lastMarker = haps.nMarkers()-1;
            extendedSegs.add( new HapSegment(hap, 0, lastMarker) );

            // permit states traversed by IBS haps
            if (ibdSegs.isEmpty()==false) {
                for (int k=0, n=ibdSegs.size(); k<n; ++k) {
                    HapSegment targetHS = ibdSegs.get(k);
                    int start = modifyStart(targetHS, tree);
                    int end  = modifyEnd(targetHS, tree);
                    extendedSegs.add(new HapSegment(targetHS.hap(), start, end));
                }
            }
            return extendedSegs;
        }

        private IntIntervalTree<HapSegment> getTree(int nMarkers,
                Collection<HapSegment> c) {
            IntIntervalTree<HapSegment> tree
                    = new CenteredIntIntervalTree<>(0, nMarkers-1);
            c.stream().forEach((hs) -> {
                tree.add(hs);
            });
            return tree;
        }

        private void convertToIndices(int marker,
                IntIntervalTree<HapSegment> tree, IndexSet set) {
            set.clear();
            Collection<HapSegment> c = new ArrayList<>(30);
            tree.intersect(marker, c);
            c.stream().forEach((hs) -> {
                set.add(hapStates[marker][hs.hap()]);
            });
        }

        @Override
        public int nMarkers() {
            return nMarkers;
        }

        @Override
        public int marker() {
            return marker;
        }

        @Override
        public void setMarker(int marker) {
            this.marker = marker;
            this.i1 = 0;
            this.i2 = 0;
            this.edge1 = -1;
            this.edge2 = -1;
            this.rev = false;
            convertToIndices(marker, tree1, indices1);
            convertToIndices(marker, tree2, indices2);
        }

        @Override
        public boolean hasNext() {
            return i1<indices1.size();
        }

        @Override
        public void next() {
            if (hasNext()==false) {
                throw new NoSuchElementException();
            }
            if (rev) {
                int tmp = edge1;
                edge1 = edge2;
                edge2 = tmp;
                ++i2;
                if (i2==indices2.size()) {
                    ++i1;
                    i2 = 0;
                }
                rev = false;
            }
            else {
                edge1 = indices1.enumeratedValue(i1);
                edge2 = indices2.enumeratedValue(i2);
                if (indices1.contains(edge2)==false
                        || indices2.contains(edge1)==false) {
                    rev = true;
                }
                else {
                    ++i2;
                    if (i2==indices2.size()) {
                        ++i1;
                        i2 = 0;
                    }
                }
            }
        }

        @Override
        public int edge1() {
            return edge1;
        }

        @Override
        public int edge2() {
            return edge2;
        }
    }
}
