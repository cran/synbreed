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

import blbutil.FileIt;
import haplotype.HapPairs;
import vcf.HapsMarker;
import haplotype.HapsMarkerIterator;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.ListIterator;

/**
 * <p>Class {@code MergeableDag} contains a static, thread-safe factory
 * method that constructs a Directed Acyclic Graph (DAG) from sequence data.
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
    private static final int MIN_DEPTH = 10;

    private final float scale;
    private final Dag dag;

    private float nUnmergedAtLeaf = 0f;

  /**
     * Constructs and returns a new {@code Dag} instance from the
     * specified data.
     * @param hapPairs the sequence data
     * @param weights an array whose {@code j}-th element is the
     * weight for the {@code j}-th haplotype
     * @param scale a parameter that multiplicatively scales the node
     * similarity threshold
     * @param nInitLevels the number of initial levels to read
     * @return a new {@code Dag} instance
     *
     * @throws IllegalArgumentException if {@code hapPairs.nMarkers() == 0}
     * @throws IllegalArgumentException if
     * {@code (weights[j] <= 0 || Float.isFinite(weights[j]) == false)}
     * for any {@code j} satisfying {@code (0 <= j && j < weights.length)}
     * @throws IllegalArgumentException if
     * {@code Double.isFinite(scale) == false || scale <= 0}
     * @throws IllegalArgumentException if {@code nInitLevels < 1}
     * @throws NullPointerException if
     * {@code hapPairs == null || weights == null}
     */
    public static Dag dag(HapPairs hapPairs, float[] weights, float scale,
            int nInitLevels) {
        MergeableDag md = new MergeableDag(hapPairs, weights, scale, nInitLevels);
        return md.dag();
    }

  /**
     * Constructs a new {@code MergeableDag} instance from the specified data.
     * @param hapPairs the sequence data
     * @param weights an array whose {@code j}-th element is the
     * weight for the {@code j}-th haplotype
     * @param scale a parameter that multiplicatively scales the node
     * similarity threshold
     * @param nInitLevels the number of initial levels to read
     *
     * @throws IllegalArgumentException if {@code hapPairs.nMarkers() == 0}
     * @throws IllegalArgumentException if
     * {@code (weights[j] <= 0 || Float.isFinite(weights[j]) == false)}
     * for any {@code j} satisfying {@code (0 <= j && j < weights.length)}
     * @throws IllegalArgumentException if
     * {@code Double.isFinite(scale) == false || scale <= 0}
     * @throws IllegalArgumentException if {@code nInitLevels < 1}
     * @throws NullPointerException if
     * {@code hapPairs == null || weights == null}
     */
    private MergeableDag(HapPairs hapPairs, float[] weights, float scale,
            int nInitLevels) {
        checkParameters(hapPairs, weights, scale, nInitLevels);
        this.scale = scale;
        List<DagLevel> mergedLevels = new ArrayList<>(hapPairs.nMarkers());
        float maxUnmerged = maxUnmergedAtLeaf(hapPairs, weights);
        int lastReadDepth = nInitLevels;
        try (FileIt<HapsMarker> it = new HapsMarkerIterator(hapPairs)) {
            MergeableDagLevel current = readFirstLevel(it, weights);
            MergeableDagLevel leaf = readLevels(it, nInitLevels, current);
            while (current.next() != null) {
                nUnmergedAtLeaf = 0f;
                current = current.next();
                mergeParentNodes(current);
                MergeableDagLevel previousLevel = current.setPreviousToNull();
                mergedLevels.add(previousLevel.toDagLevel());

                if (it.hasNext()) {
                    float ratio = (nUnmergedAtLeaf / maxUnmerged);
                    int depth = (leaf.index() - current.index());
                    int readDepth = nextReadDepth(ratio, depth, lastReadDepth);
                    if (readDepth>depth) {
                        leaf = readLevels(it, (readDepth - depth), leaf);
                        lastReadDepth = readDepth;
                    }
                }
            }
            mergedLevels.add(current.toDagLevel());
            DagLevel[] levels =  mergedLevels.toArray(new DagLevel[0]);
            this.dag = new ImmutableDag(hapPairs.markers(), levels);
        }
    }

    private static void checkParameters(HapPairs hapPairs, float[] weights,
            float scale, int nInitLevels) {
        if (nInitLevels < 1) {
            throw new IllegalArgumentException(String.valueOf(nInitLevels));
        }
        if (hapPairs.nMarkers()==0) {
            throw new IllegalArgumentException("hapPairs.nMarkers()==0");
        }
        if (weights!=null) {
            for (int j=0; j<weights.length; ++j) {
                float f = weights[j];
                if (f <= 0.0 || Float.isFinite(f) == false) {
                    throw new IllegalArgumentException(String.valueOf(f));
                }
            }
        }
        if (scale <= 0 || Float.isFinite(scale) == false) {
            throw new IllegalArgumentException(String.valueOf(scale));
        }
    }

    private static MergeableDagLevel readFirstLevel(FileIt<HapsMarker> it,
            float[] weights) {
        HapsMarker hapsMarker = it.next();
        if (weights==null) {
                return new MergeableDagLevel(hapsMarker);
        }
        else {
            return new MergeableDagLevel(hapsMarker, weights);
        }
    }

    private static float maxUnmergedAtLeaf(HapPairs hapPairs, float[] weights) {
        float maxPropUnmerged = 0.01f;
        float sum = (weights==null ? hapPairs.nHaps() : sum(weights));
        return maxPropUnmerged * sum;
    }

    private static int nextReadDepth(float unmergedRatio, int depth,
            int lastDepth) {
        if (unmergedRatio <= 1) {
            return MIN_DEPTH;
        }
        else if (depth < (0.85*lastDepth) ) {
            return 1 + (int) Math.round(0.95*lastDepth);
        }
        else if ( (unmergedRatio>2) && (depth > (0.95*lastDepth)) ) {
            return (int) Math.ceil((1 + unmergedRatio/20)*lastDepth);
        }
        else {
            return lastDepth;
        }
    }

    private static MergeableDagLevel readLevels(FileIt<HapsMarker> it,
            int nLevelsToRead, MergeableDagLevel leafLevel) {
        for (int j=0; it.hasNext() && j<nLevelsToRead; ++j) {
            HapsMarker hapsMarker = it.next();
            MergeableDagLevel newLeaf =
                    new MergeableDagLevel(leafLevel, hapsMarker);
            leafLevel.setNextLevel(newLeaf);
            leafLevel = newLeaf;
        }
        return leafLevel;
    }

    private static float sum(float[] fa) {
        float sum = 0f;
        for (float f : fa) {
            sum += f;
        }
        return sum;
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
        Score minScore = MAX_SCORE;
        SortedNodes parentNodes = sortedParents(level);
        int[] parents = parentNodes.sorted;
        int nParentsWithSibs = parentNodes.nWithSibs;
        for (int j=0; j<nParentsWithSibs; ++j) {
            int nodeA = parents[j];
            for (int k=j+1; k<parents.length; ++k) {
                int nodeB = parents[k];
                Score s = score(level, nodeA, nodeB);
                if (s!=null) {
                    if (s.score()<minScore.score() && s.isMergeable()) {
                        minScore = s;
                    }
                    scores.add(s);
                }
            }
        }
        return minScore;
    }

    private static class SortedNodes {

        public int[] sorted;
        public int nWithSibs;

        public SortedNodes(int[] sorted, int nWithSibs) {
            this.sorted = sorted;
            this.nWithSibs = nWithSibs;
        }
    }

    private SortedNodes sortedParents(MergeableDagLevel level) {
        int[] nodeArray = level.parentNodeArray();
        int index1 = 0;
        int index2 = nodeArray.length-1;
        while (index1<index2) {
            if (level.hasSibling(nodeArray[index1])) {
                ++index1;
            }
            else {
                int tmp = nodeArray[index1];
                nodeArray[index1] = nodeArray[index2];
                nodeArray[index2--] = tmp;
            }
        }
        if (level.hasSibling(nodeArray[index1])) {
            ++index1;
        }
        return new SortedNodes(nodeArray, index1);
    }

    private Score score(MergeableDagLevel level, int nodeA, int nodeB) {
        float maxDiff = 0.0f;
        float nodeCntA = level.nodeCount(nodeA);
        float nodeCntB = level.nodeCount(nodeB);
        float threshold = (float) (scale*Math.sqrt((1.0/nodeCntA)+(1.0/nodeCntB)));
        maxDiff = similar(level, nodeA, nodeB,
                nodeCntA, nodeCntB, level.index(),
                nodeCntA, nodeCntB, maxDiff, threshold);
        if (maxDiff > MAX_THRESHOLD_RATIO*threshold) {
            return null;
        }
        else {
            boolean isMergeable = (maxDiff < threshold);
            return new Score(nodeA, nodeB, maxDiff, isMergeable);
        }
    }

    /*
     * Returns a similarity-score (lower scores correspond to higher
     * similarity).
     *
     * @param marker marker DAG marker containing the specified nodes
     * @param nodeA a parent node at the specified DAG level in tree A
     * @param nodeB a parent node at the specified DAG level in tree B
     * @param nodeCntA the count of the parent node in tree A
     * @param nodeCntB the count of the parent node in tree B
     * @param baseMarker the marker index at the root of trees A and B
     * @param nA the node count of the root of tree A
     * @param nB the node count of the root of tree B
     * @param maxDiff the current maximum difference in proportions in
     * the counts of corresponding tree branches
     * @param threshold the maximum permitted node similarity
     * @return a similarity-score
     */
    private float similar(MergeableDagLevel level,
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
        if (nodeA == -1 ^ nodeB == -1) {
            return maxDiff;
        }
        else if (level==null) {
            nUnmergedAtLeaf += (nodeCntA + nodeCntB);
            return maxDiff;
        }
        for (int j=0, n=level.nAlleles(); j<n; ++j) {
            int edgeA = level.outEdge(nodeA, j);
            int edgeB = level.outEdge(nodeB, j);
            int childA = (edgeA != -1) ? level.childNode(edgeA) : -1;
            int childB = (edgeB != -1) ? level.childNode(edgeB) : -1;
            nodeCntA = (edgeA != -1) ? level.edgeCount(edgeA) : 0.0f;
            nodeCntB = (edgeB != -1) ? level.edgeCount(edgeB) : 0.0f;
            float childMaxDiff = similar(level.next(), childA, childB,
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
     * Returns the constructed DAG.
     * @return the constructed DAG
     */
    public Dag dag() {
        return dag;
    }

    /**
     * Returns a string description of {@code this}.  The
     * exact details of the representation are unspecified and subject
     * to change.
     * @return a string description of {@code this}
     */
    @Override
    public String toString() {
        return dag.toString();
    }
}
