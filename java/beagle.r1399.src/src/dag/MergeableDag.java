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
package dag;

import blbutil.FileIterator;
import blbutil.Utilities;
import haplotype.HapPairs;
import haplotype.HapsMarker;
import haplotype.HapsMarkerIterator;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;
import vcf.Marker;
import vcf.Markers;

/**
 * <p>Class {@code MergeableDag} constructs a Directed Acyclic Graph (DAG)
 * from sequence data.
 * </p>
 *
 * References:
 * <br>
 * Ron D, Singer Y, and Tishby N (1998) On the Learnability and
 * usage of acyclic probabilistic finite automata.  Journal of Computer
 * and SystemSciences 56:133-152.
 * <br>
 * Browning S (2006) Multi-locus association mapping using variable length
 * Markov chains.  Am J Hum Genet. 78:903-913.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class MergeableDag {

    private static final Score MAX_SCORE =
            new Score(-1, -1, Float.POSITIVE_INFINITY, false);
    private static final float MAX_THRESHOLD_RATIO = 1.4f;

    private final FileIterator<HapsMarker> it;
    private final float[] weights;
    private final int minWindow;
    private final int maxWindow;
    private final float scale;

    private final Dag dag;

    /**
     * Constructs and returns a new {@code Dag} instance.
     * @param hapPairs the sequence data.
     * @param weights an array whose {@code j}-th element is the
     * weight for the {@code j}-th haplotype.
     * @param maxWindow maximum window size used when constructing the DAG.
     * @param scale a parameter that multiplicatively scales the node
     * similarity threshold.
     * @return a new {@code Dag} instance.
     *
     * @throws IllegalArgumentException if {@code hapPairs.nMarkers()==0}
     * @throws IllegalArgumentException if any element of the
     * {@code weights} array is non-positive, NaN or infinite
     * @throws IllegalArgumentException if {@code maxWindow<1}
     * @throws IllegalArgumentException if
     * {@code Double.isInfinite(scale) || Double.isNaN(scale)
     *      || scale<=0}
     * @throws NullPointerException if {@code hapPairs==null || weights==null}
     */
    public static Dag dag(HapPairs hapPairs, float[] weights, int maxWindow,
            float scale) {
        try (FileIterator<HapsMarker> tmpIt = new HapsMarkerIterator(hapPairs)) {
            MergeableDag md = new MergeableDag(tmpIt, hapPairs.markers(),
                    weights, maxWindow, scale);
            return md.dag();
        }
    }

    /**
     * Constructs and returns a new {@code Dag} instance.
     * @param it an iterator that returns per-marker phased genotype data.
     * @param weights an array whose {@code j}-th element is the
     * weight for the {@code j}-th haplotype.
     * @param maxWindow maximum window size used when constructing the DAG.
     * @param scale a parameter that multiplicatively scales the node
     * similarity threshold.
     * @return a new {@code Dag} instance.
     *
     * @throws IllegalArgumentException if {@code it.hasNext()==false}
     * @throws IllegalArgumentException if any element of the
     * {@code weights} array is non-positive, NaN or infinite
     * @throws IllegalArgumentException if {@code maxWindow<1}
     * @throws IllegalArgumentException if
     * {@code Double.isInfinite(scale) || Double.isNaN(scale)
     *      || scale<=0}
     * @throws NullPointerException if {@code it==null || weights==null}
     */
    public static Dag dag(FileIterator<HapsMarker> it, float[] weights,
            int maxWindow, float scale) {
        Markers markers = null;
        return new MergeableDag(it, markers, weights, maxWindow, scale).dag();
    }

    /**
     * Constructs a new {@code MergeableDag} instance.
     * @param it an iterator that returns per-marker phased genotype data.
     * @param markers the list of markers that correspond to the
     * {@code HapsMarker} instances returned by {@code it}, or {@code null}
     * if the list of markers should be constructed from the
     * {@code HapsMarker} instances returned by {@code it}.
     * @param weights an array whose {@code j}-th element is the
     * weight for the {@code j}-th haplotype.
     * @param maxWindow maximum window size used when constructing the DAG.
     * @param scale a parameter that multiplicatively scales the node
     * similarity threshold.
     *
     * @throws IllegalArgumentException if {@code it.hasNext()==false}
     * @throws IllegalArgumentException if any element of the
     * {@code weights} array is non-positive, NaN or infinite
     * @throws IllegalArgumentException if {@code maxWindow<1}
     * @throws IllegalArgumentException if
     * {@code Double.isInfinite(scale) || Double.isNaN(scale)
     *      || scale<=0}
     * @throws NullPointerException if {@code it==null || weights==null}
     */
    private MergeableDag(FileIterator<HapsMarker> it, Markers markers,
            float[] weights, int maxWindow, float scale) {
        checkParameters(it, weights, maxWindow, scale);
        this.it = it;
        this.weights = weights.clone();
        this.minWindow = maxWindow/12 + 1;
        this.maxWindow = maxWindow;
        this.scale = scale;

        List<DagLevel> mergedLevels = new ArrayList<>(25000);
        MergeableDagLevel currentLevel = new MergeableDagLevel(it.next(),
                weights);
        readLevels(it, weights, currentLevel, minWindow);
        while (currentLevel.next() != null) {
            currentLevel = currentLevel.next();
            mergeParentNodes(currentLevel);
            MergeableDagLevel previousLevel = currentLevel.setPreviousToNull();
            mergedLevels.add(previousLevel.toDagLevel());
        }
        mergedLevels.add(currentLevel.toDagLevel());

        if (markers==null) {
            markers = markers(mergedLevels);
        }
        DagLevel[] levels =  mergedLevels.toArray(new DagLevel[0]);
        this.dag = new ImmutableDag(markers, levels);
    }

    private static Markers markers(List<DagLevel> list) {
        Marker[] ma = new Marker[list.size()];
        for (int j=0; j<ma.length; ++j) {
            ma[j] = list.get(j).marker();
        }
        return new Markers(ma);
    }

    /**
     * Returns the constructed DAG.
     * @return the constructed DAG.
     */
    public Dag dag() {
        return dag;
    }

    private static void checkParameters(FileIterator<HapsMarker> it,
            float[] weights, int maxWindow, double scale) {
        if (it.hasNext()==false) {
            throw new IllegalArgumentException("it.hasNext()==false");
        }
        for (int j=0; j<weights.length; ++j) {
            float f = weights[j];
            if (f <= 0.0 || Float.isNaN(f) || Float.isInfinite(f)) {
                throw new IllegalArgumentException("weights[" + j + "]=" + f);
            }
        }
        if (maxWindow<1) {
            throw new IllegalArgumentException("maxWindow: " + maxWindow);
        }
        if (Double.isInfinite(scale) || Double.isNaN(scale) || scale <=0) {
            throw new IllegalArgumentException("scale: " + scale);
        }
    }

    private static void readLevels(FileIterator<HapsMarker> it, float[] weights,
            MergeableDagLevel leafLevel, int maxLevelsToRead) {
        for (int j=0; j<maxLevelsToRead && it.hasNext(); ++j) {
            MergeableDagLevel newLevel = new MergeableDagLevel(leafLevel,
                    it.next(), weights);
            leafLevel.setNextLevel(newLevel);
            leafLevel = newLevel;
        }
    }

    private void mergeParentNodes(MergeableDagLevel level) {
        List<Score> scores = new LinkedList<>();
        Score min = getPairwiseScores(level, scores);
        while (min.isMergeable()) {
            int retainedNode = min.nodeA();
            int removedNode = min.nodeB();
            if (level.hasSibling(retainedNode)==false) {
                // Ensure that no-sibling nodes are always removed
                retainedNode = min.nodeB();
                removedNode = min.nodeA();
                assert level.hasSibling(retainedNode);
            }
            else if (level.hasSibling(removedNode)
                    && level.nodeCount(min.nodeA())<level.nodeCount(min.nodeB())) {
                removedNode = min.nodeB();
                retainedNode = min.nodeA();
            }
            level.mergeParentNodes(retainedNode, removedNode);

            min = MAX_SCORE;
            ListIterator<Score> scoreIt = scores.listIterator();
            while (scoreIt.hasNext()) {
                Score s = scoreIt.next();
                if (s.nodeA()==removedNode || s.nodeB()==removedNode) {
                    scoreIt.remove();
                }
                else {
                    if (s.nodeA()==retainedNode || s.nodeB()==retainedNode) {
                        s = score(level, s.nodeA(), s.nodeB());
                        if (s!=null) {
                            if (s.score()<min.score() && s.isMergeable()) {
                                min = s;
                            }
                            scoreIt.set(s);
                        }
                        else {
                            scoreIt.remove();
                        }
                    }
                    else if (s.score()<min.score() && s.isMergeable()) {
                        min = s;
                    }
                }
            }
        }
    }

    private Score getPairwiseScores(MergeableDagLevel level,
            List<Score> scores) {
        if (level.next()==null) {
            readLevels(it, weights, level, minWindow);
        }
        Score minScore = MAX_SCORE;
        int[] nodeArray = level.parentNodeArray();
        boolean[] hasSibling = hasSibling(level, nodeArray);
        for (int j=0; j<nodeArray.length; ++j) {
            int nodeA = nodeArray[j];
            for (int k=j+1; k<nodeArray.length; ++k) {
                int nodeB = nodeArray[k];
                if (hasSibling[j] || hasSibling[k]) {
                    Score s = score(level, nodeA, nodeB);
                    if (s!=null) {
                        if (s.score()<minScore.score() && s.isMergeable()) {
                            minScore = s;
                        }
                        scores.add(s);
                    }
                }

            }
        }
        return minScore;
    }

    private boolean[] hasSibling(MergeableDagLevel level,
            int[] parentNodeArray) {
        boolean[] hasSibling = new boolean[parentNodeArray.length];
        for (int j=0; j<parentNodeArray.length; ++j) {
            hasSibling[j] = level.hasSibling(parentNodeArray[j]);
        }
        return hasSibling;
    }

    private Score score(MergeableDagLevel level, int nodeA, int nodeB) {
        float maxDiff = 0.0f;
        float nodeCntA = level.nodeCount(nodeA);
        float nodeCntB = level.nodeCount(nodeB);
        float threshold = (float) (scale*Math.sqrt((1.0/nodeCntA)+(1.0/nodeCntB)));
        maxDiff = similar(level.previous(), level, nodeA, nodeB,
                nodeCntA, nodeCntB, level.markerIndex(),
                nodeCntA, nodeCntB, maxDiff, threshold);
        if (maxDiff > MAX_THRESHOLD_RATIO*threshold) {
            return null;
        }
        else {
            boolean isMergeable = (maxDiff < threshold);
            return new Score(nodeA, nodeB, maxDiff, isMergeable);
        }
    }

    private MergeableDagLevel nextLevel(MergeableDagLevel prevLevel,
            int baseMarker, float propA, float propB, float maxDiff,
            float threshold) {
        float t1 = 0.7f * threshold;
        float t2 = 0.5f * threshold;
        int depth = prevLevel.markerIndex() - baseMarker;
        if (depth<maxWindow && ( (maxDiff>t1 && Math.min(propA, propB)>t2)
                                 || depth<minWindow) ) {
//        float maxProp = Math.max(propA, propB);
//        float minProp = Math.min(propA, propB);
//        float estMaxDiff = maxProp - 0.5f*minProp;
//        int depth = prevLevel.markerIndex() - baseMarker;
//        if ( (depth<maxWindow && threshold<0.4 && estMaxDiff>1.1*maxDiff
//                && estMaxDiff>threshold) || depth<minWindow) {
            MergeableDagLevel newLeaf =
                    new MergeableDagLevel(prevLevel, it.next(), this.weights);
            prevLevel.setNextLevel(newLeaf);
            return newLeaf;
        }
        else {
            return null;
        }
    }

    /*
     * Returns a similarity-score (lower scores correspond to higher
     * similarity).
     *
     * @param marker marker DAG marker containing the specified nodes
     * @param nodeA a parent node at the specified DAG marker in tree A.
     * @param nodeB a parent node at the specified DAG marker in tree B.
     * @param nodeCntA the count of the parent node in tree A.
     * @param nodeCntB the count of the parent node in tree B.
     * @param baseMarker the marker index at the root of trees A and B
     * @param nA the node count of the root of tree A.
     * @param nB the node count of the root of tree B
     * @param maxDiff the current maximum difference in proportions in
     * the counts of corresponding tree branches.
     * @param threshold the maximum permitted node similarity.
     *
     * @return a similarity-score. Lower scores correspond to greater
     * similarity.
     */
    private float similar(MergeableDagLevel prevLevel, MergeableDagLevel level,
            int nodeA, int nodeB, float nodeCntA, float nodeCntB,
            int baseMarker, float nA, float nB, float maxDiff, float threshold) {
        float propA = nodeCntA / nA;
        float propB = nodeCntB / nB;
        float diff = Math.abs(propA - propB);
        if (diff >= threshold) {
            return diff;
        }
        else if (propA <= maxDiff && propB <= maxDiff) {
            return maxDiff;
        }
        else if (diff > maxDiff) {
            maxDiff = diff;
        }
        if (level==null && it.hasNext()) {
            level = nextLevel(prevLevel, baseMarker, propA, propB, maxDiff,
                    threshold);
        }
        if (nodeA == -1 || nodeB == -1 || level==null) {
            return maxDiff;
        }
        for (int j=0, n=level.nAlleles(); j<n; ++j) {
            int edgeA = level.outEdge(nodeA, j);
            int edgeB = level.outEdge(nodeB, j);
            int childA = (edgeA != -1) ? level.childNode(edgeA) : -1;
            int childB = (edgeB != -1) ? level.childNode(edgeB) : -1;
            nodeCntA = (edgeA != -1) ? level.edgeCount(edgeA) : 0.0f;
            nodeCntB = (edgeB != -1) ? level.edgeCount(edgeB) : 0.0f;
            float childMaxDiff = similar(level, level.next(), childA, childB,
                    nodeCntA, nodeCntB, baseMarker, nA, nB, maxDiff, threshold);
            if (childMaxDiff > maxDiff) {
                if (childMaxDiff >= threshold) {
                    return childMaxDiff;
                }
                else {
                    maxDiff = childMaxDiff;
                }
            }
        }
        return maxDiff;
    }

    /**
     * Returns a string description of {@code this}.  The
     * exact details of the representation are unspecified and subject
     * to change.
     * @return a string description of {@code this}.
     */
    @Override
    public String toString() {
        return dag.toString();
    }
}
