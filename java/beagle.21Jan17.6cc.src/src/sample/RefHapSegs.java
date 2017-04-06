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

import haplotype.SampleHapPairs;
import java.util.stream.IntStream;

/**
 * <p>Class {@code RefHapSegs} represents reference haplotypes that span
 * segments determined by non-overlapping clusters of markers.
 * </p>
 * <p>Instances of class {@code RefHapSegs} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class RefHapSegs {

    private final int[] segStart;
    private final int[] segEnd;
    private final SampleHapPairs refHapPairs;
    private final RefHapSeg[] refHapSegs;

    /**
     * Constructs a new {@code RefHapSegs} instance from the specified data.
     * @param refHapPairs the reference haplotype pairs
     * @param segStart an array whose {@code j}-th element is the
     * starting reference marker index (inclusive) for the {@code j}-th
     * segment of markers
     * @param segEnd an array whose {@code j}-th element is the
     * ending reference marker index (exclusive) for the {@code j}-th
     * segment of markers
     * @throws IllegalArgumentException if
     * {@code segStart.length != segEnd.length}
     * @throws IllegalArgumentException if
     * {@code (segStart[j] < 0 || segStart[j] >= segEnd[j]
     * || segEnd[j] >= refHapPairs.nMarkers())} for some {@code j} satisfying
     * {@code 0 <= j && j < segStart.length}
     * @throws NullPointerException if
     * {@code refHapPairs == null || segStart == null || segEnd == null}
     */
    public RefHapSegs(SampleHapPairs refHapPairs, int[] segStart, int[] segEnd) {
        int nMarkers = refHapPairs.nMarkers();
        checkClusters(segStart, segEnd, nMarkers);
        this.segStart = segStart.clone();
        this.segEnd = segEnd.clone();
        this.refHapPairs = refHapPairs;
        this.refHapSegs = IntStream.range(0, this.segStart.length)
                .parallel()
                .mapToObj(j -> new RefHapSeg(refHapPairs, segStart[j], segEnd[j]))
                .toArray(RefHapSeg[]::new);
    }

    private void checkClusters(int[] starts, int[] ends, int nMarkers) {
        if (starts.length != ends.length) {
            throw new IllegalArgumentException("inconsistent data");
        }
        for (int j=0; j<starts.length; ++j) {
            if (starts[j] < 0 || starts[j] >= ends[j] || ends[j] > nMarkers) {
                throw new IllegalArgumentException("inconsistent data");
            }
        }
    }

    /**
     * Returns the reference haplotype pairs.
     * @return the reference haplotype pairs
     */
    public SampleHapPairs refHapPairs() {
        return refHapPairs;
    }

    /**
     * Return the number of distinct reference allele sequences in the
     * specified chromosome segment.
     * @param segment index of a chromosome segment determined by
     * the marker clusters
     * @return the number of distinct reference allele sequences in the
     * specified chromosome segment
     * @throws IndexOutOfBoundsException if
     * {@code segment < 0 || segment > this.nClusters()}
     */
    public int nSeq(int segment) {
        return refHapSegs[segment].nSeq();
    }

    /**
     * Return the number of markers in the specified chromosome segment.
     * @param segment index of a chromosome segment determined by
     * the marker clusters
     * @return the number of markers in the specified chromosome segment
     * @throws IndexOutOfBoundsException if
     * {@code segment < 0 || segment > this.nClusters()}
     */
    public int nMarkers(int segment) {
        return refHapSegs[segment].end() - refHapSegs[segment].start();
    }

    /**
     * Return the index of the allele sequence in the specified chromosome
     * segment for the specified reference haplotype.
     *
     * @param segment index of a chromosome segment determined by
     * the marker clusters
     * @param hap a haplotype index
     *
     * @return the index of the allele sequence in the specified chromosome
     * segment for the specified reference haplotype
     *
     * @throws IndexOutOfBoundsException if
     * {@code segment < 0 || segment > this.nClusters()}
     * @throws IndexOutOfBoundsException if
     * {@code hap < 0 || hap >= this.refHapPairs().nHaps()}
     */
    public int seq(int segment, int hap) {
        return refHapSegs[segment].seq(hap);
    }

    /**
     * Return the specified reference haplotype allele.
     *
     * @param segment index of a chromosome segment determined by
     * the marker clusters
     * @param marker index of a marker in the specified interval
     * @param seq index of a reference allele sequence in the specified
     * interval
     * @return the specified reference haplotype allele
     *
     * @throws IndexOutOfBoundsException if
     * {@code segment < 0 || segment > this.nClusters()}
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers(interval)}
     * @throws IndexOutOfBoundsException if
     * {@code seq < 0 || seg >= this.nSeq(segment)}
     */
    public int allele(int segment, int marker, int seq) {
        return refHapSegs[segment].allele(marker, seq);
    }

    /**
     * Returns the number of marker segments..
     * @return the number of marker segments
     */
    public int nSegs() {
        return segStart.length;
    }

    /**
     * Returns the index of the first marker (inclusive) in the specified
     * marker segment.
     * @param segment an index of a marker segment
     * @return the index of the first marker (inclusive) in the specified
     * marker segment
     * @throws IndexOutOfBoundsException if
     * {@code segment < 0 || segment >= this.nSegs()}
     */
    public int segStart(int segment) {
        return segStart[segment];
    }

    /**
     * Returns the index of the last marker (exclusive) in the specified
     * marker segment.
     * @param segment an index of a marker segment
     * @return the index of the last marker (exclusive) in the specified
     * marker segment
     * @throws IndexOutOfBoundsException if
     * {@code segment 0  || segment >= this.nSegs()}
     */
    public int segEnd(int segment) {
        return segEnd[segment];
    }
}
