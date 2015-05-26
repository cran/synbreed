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
package vcf;

import beagleutil.ChromIds;
import blbutil.Const;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * <p>Class {@code Markers} represent a list of markers in chromosome order.
 * </p>
 * <p>Instances of class {@code Markers} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class Markers {

    private final Set<Marker> markerSet;

    private final Marker[] fwdMarkers;
    private final int[] fwdSumAlleles;
    private final int[] fwdSumGenotypes;
    private final int[] fwdSumHaplotypeBits;
    private final int fwdHashCode;

    private final Marker[] bwdMarkers;
    private final int[] bwdSumAlleles;
    private final int[] bwdSumGenotypes;
    private final int[] bwdSumHaplotypeBits;
    private final int bwdHashCode;

    /**
     * Construct a new {@code Markers} instance that represents the
     * specified list of markers.
     * @param markers a list of markers in chromosome order.
     *
     * @throws IllegalArgumentException if markers on a chromosome are not
     * in chromosome order.
     * @throws IllegalArgumentException if there are duplicate markers.
     * @throws IllegalArgumentException if the markers on a chromosome
     * do not form a contiguous set of entries within the array.
     *
     * @throws NullPointerException if
     * {@code markers==null || markers[j]==null} for any
     * {@code 0<=j<markers.length}.
     */
    public Markers(Marker[] markers) {
        checkMarkerPosOrder(markers);
        this.fwdMarkers = markers.clone();
        this.bwdMarkers = reverse(this.fwdMarkers);
        this.markerSet = markerSet(fwdMarkers);

        this.fwdSumAlleles = cumSumAlleles(fwdMarkers);
        this.fwdSumGenotypes = cumSumGenotypes(fwdMarkers);
        this.fwdSumHaplotypeBits = cumSumHaplotypeBits(fwdMarkers);
        this.fwdHashCode = Arrays.deepHashCode(fwdMarkers);

        this.bwdSumAlleles = cumSumAlleles(bwdMarkers);
        this.bwdSumGenotypes = cumSumGenotypes(bwdMarkers);
        this.bwdSumHaplotypeBits = cumSumHaplotypeBits(bwdMarkers);
        this.bwdHashCode = Arrays.deepHashCode(bwdMarkers);
    }

    /*
     * Constructs a new {@code Markers} instance that has the same
     * markers as the specified {@code Markers} object, but with
     * order reversed.
     */
    private Markers(Markers markers) {
        this.markerSet = markers.markerSet;
        this.fwdMarkers = markers.bwdMarkers;
        this.bwdMarkers = markers.fwdMarkers;

        this.fwdSumAlleles = markers.bwdSumAlleles;
        this.fwdSumGenotypes = markers.bwdSumGenotypes;
        this.fwdSumHaplotypeBits = markers.bwdSumHaplotypeBits;
        this.fwdHashCode = markers.bwdHashCode;

        this.bwdSumAlleles = markers.fwdSumAlleles;
        this.bwdSumGenotypes = markers.fwdSumGenotypes;
        this.bwdSumHaplotypeBits = markers.fwdSumHaplotypeBits;
        this.bwdHashCode = markers.fwdHashCode;
    }

    private static void checkMarkerPosOrder(Marker[] markers) {
        if (markers.length < 2) {
            return;
        }
        Set<Integer> chromIndices = new HashSet<>();
        chromIndices.add(markers[0].chromIndex());
        chromIndices.add(markers[1].chromIndex());
        for (int j=2; j<markers.length; ++j) {
            int chr0 = markers[j-2].chromIndex();
            int chr1 = markers[j-1].chromIndex();
            int chr2 = markers[j].chromIndex();
            if (chr0 == chr1 && chr1==chr2) {
                int pos0 = markers[j-2].pos();
                int pos1 = markers[j-1].pos();
                int pos2 = markers[j].pos();
                if ( (pos1<pos0 && pos1<pos2) || (pos1>pos0 && pos1>pos2) ) {
                    String s = "markers not in chromosomal order: "
                            + Const.nl + markers[j-2]
                            + Const.nl + markers[j-1]
                            + Const.nl + markers[j];
                    throw new IllegalArgumentException(s);
                }
            }
            else if (chr1!=chr2) {
                if (chromIndices.contains(chr2)) {
                    String s = "markers on chromosome are not contiguous: "
                            + ChromIds.instance().id(chr2);
                    throw new IllegalArgumentException(s);
                }
                chromIndices.add(chr2);
            }
        }
    }

    private static Marker[] reverse(Marker[] markers) {
        int lastIndex = markers.length - 1;
        Marker[] rev = new Marker[markers.length];
        for (int j=0; j<markers.length; ++j) {
            rev[j] = markers[lastIndex - j];
        }
        return rev;
    }

    private static Set<Marker> markerSet(Marker[] markers) {
        Set<Marker> markerSet = new HashSet<>(markers.length);
        for (Marker m : markers) {
            if (markerSet.add(m)==false) {
                throw new IllegalArgumentException("Duplicate marker: " + m);
            }
        }
        return markerSet;
    }

    private static int[] cumSumAlleles(Marker[] markers) {
        int[] ia = new int[markers.length + 1];
        for (int j=1; j<ia.length; ++j) {
            ia[j] = ia[j-1] + markers[j-1].nAlleles();
        }
        return ia;
    }

    private static int[] cumSumGenotypes(Marker[] markers) {
        int[] ia = new int[markers.length + 1];
        for (int j=1; j<ia.length; ++j) {
            ia[j] = ia[j-1] + markers[j-1].nGenotypes();
        }
        return ia;
    }

    private static int[] cumSumHaplotypeBits(Marker[] markers) {
        int[] ia = new int[markers.length + 1];
        for (int j=1; j<ia.length; ++j) {
            int nAllelesM1 = markers[j-1].nAlleles() - 1;
            int nStorageBits = Integer.SIZE
                    - Integer.numberOfLeadingZeros(nAllelesM1);
            ia[j] = ia[j-1] + nStorageBits;
        }
        return ia;
    }

    /**
     * Returns a hash code value for the object.
     * @return a hash code value for the object.
     */
    @Override
    public int hashCode() {
        return fwdHashCode;
    }

    /**
     * Returns {@code true} if {@code this} and the specified object
     * represent an identical list of markers, and returns
     * {@code false} otherwise.  Markers are considered equal if
     * the {@code Marker.equals()} method returns {@code true}.
     *
     * @param obj the object to be tested for equality with {@code this}.
     *
     * @return {@code true} if {@code this} and the specified object
     * represent an identical list of markers, and returns
     * {@code false} otherwise.
     */
    @Override
    public boolean equals(Object obj) {
        if (this==obj) {
            return true;
        }
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Markers other = (Markers) obj;
        return Arrays.deepEquals(this.fwdMarkers, other.fwdMarkers);
    }

    /**
     * Constructs and returns a new {@code Markers} instance that is
     * obtained by reversing the order of markers in {@code this}.
     * @return a new {@code Markers} instance that is obtained by
     * reversing the order of markers in {@code this}.
     */
    public Markers reverse() {
        return new Markers(this);
    }

    /**
     * Returns the number of markers.
     * @return the number of markers.
     */
    public int nMarkers() {
        return fwdMarkers.length;
    }

    /**
     * Returns the specified marker.
     * @param marker a marker index.
     * @return the specified marker.
     * @throws IndexOutOfBoundsException if
     * {@code marker < 0 || marker >= this.nMarkers()}.
     */
    public Marker marker(int marker) {
        return fwdMarkers[marker];
    }

    /**
     * Returns the list of markers.
     * @return the list of markers.
     */
    public Marker[] markers() {
        return fwdMarkers.clone();
    }

    /**
     * Returns {@code true} if the specified marker is an element of
     * the list of markers represented by {@code this}, and returns
     * {@code false} otherwise.  Returns {@code false} if
     * {@code marker==null}.
     *
     * @param marker a marker.
     *
     * @return {@code true} if the specified marker is an element of
     * the list of markers represented by {@code this}, and returns
     * {@code false} otherwise.
     */
    public boolean contains(Marker marker) {
        return markerSet.contains(marker);
    }

    /**
     * Returns a {@code Markers} instance that represents
     * the specified range of marker indices.
     * @param start the starting marker index (inclusive).
     * @param end the ending marker index (exclusive).
     * @return a {@code Markers} instance that represents
     * the specified range of marker indices.
     *
     * @throws IndexOutOfBoundsException if
     * {@code start<0 || end>this.nMarkers()}.
     * @throws IllegalArgumentException if {@code start>=end}.
     */
    public Markers restrict(int start, int end) {
        if (end > fwdMarkers.length) {
            throw new IndexOutOfBoundsException("end > this.nMarkers(): " + end);
        }
        return new Markers(Arrays.copyOfRange(fwdMarkers, start, end));
    }

    /**
     * Returns the sum of the number of alleles for
     * the markers with index less than the specified index.
     * @param marker a marker index.
     * @return the sum of the number of alleles for
     * the markers with index less than the specified index.
     * @throws IndexOutOfBoundsException if
     * {@code marker<0 || marker>this.nMarkers()}
     */
    public int sumAlleles(int marker) {
        return fwdSumAlleles[marker];
    }

    /**
     * Returns {@code this.sumAlleles(this.nMarkers())}.
     * @return {@code this.sumAlleles(this.nMarkers())}.
     */
    public int sumAlleles() {
        return fwdSumAlleles[fwdMarkers.length];
    }

    /**
     * Returns the sum of the number of possible genotypes for the markers
     * with index less than the specified index.
     * @param marker a marker index.
     * @return the sum of the number of possible genotypes for the markers
     * with index less than the specified index.
     * @throws IndexOutOfBoundsException if
     * {@code marker<0 || marker>this.nMarkers()}
     */
    public int sumGenotypes(int marker) {
        return fwdSumGenotypes[marker];
    }

    /**
     * Returns {@code this.sumGenotypes(this.nMarkers())}.
     * @return {@code this.sumGenotypes(this.nMarkers())}.
     */
    public int sumGenotypes() {
        return fwdSumGenotypes[fwdMarkers.length];
    }

    /**
     * Returns the number of bits requires to store a haplotype for the
     * markers with index less than the specified index.
     * @param marker a marker index.
     * @return the number of bits requires to store a haplotype for the
     * markers with index less than the specified index.
     * @throws IndexOutOfBoundsException if
     * {@code marker<0 || marker>this.nMarkers()}
     */
    public int sumHaplotypeBits(int marker) {
        return fwdSumHaplotypeBits[marker];
    }

    /**
     * Returns {@code this.sumHaplotypeBits(this.nMarkers())}.
     * @return {@code this.sumHaplotypeBits(this.nMarkers())}.
     */
    public int sumHaplotypeBits() {
        return fwdSumHaplotypeBits[fwdMarkers.length];
    }

    /**
     * Returns a string representation of {@code this}.
     * The exact details of the representation are unspecified and
     * subject to change.
     * @return a string representation of {@code this}.
     */
    @Override
    public String toString() {
        return Arrays.toString(fwdMarkers);
    }
}
