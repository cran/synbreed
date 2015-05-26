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

import beagleutil.Samples;
import haplotype.HapPair;
import java.util.List;

/**
 * Interface {@code Data} represents a sliding marker window
 * for reference sample data and non-reference sample data.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface Data {

    /**
     * Returns {@code true} if the sliding marker window is the last
     * marker window for the chromosome and returns {@code false}
     * otherwise.
     * @return {@code true} if the sliding marker window is the last
     * marker window for the chromosome and {@code false} otherwise.
     */
    public boolean lastWindowOnChrom();

    /**
     * Returns {@code true} if the sliding marker window can advance
     * and returns {@code false} otherwise.
     * @return {@code true} if the sliding marker window can advance
     * and returns {@code false} otherwise.
     */
    boolean canAdvanceWindow();

    /**
     * Advances the sliding marker window.
     *
     * @param overlap the requested number of markers of overlap between the
     * marker window immediately before the method invocation and the
     * marker window immediately after the method returns.
     * If the marker window immediately before method invocation contains
     * the last markers on a chromosome, then there will be 0 markers
     * in the overlap.
     * @param windowSize the requested number of the markers in the window
     * immediately after the method returns.
     *
     * @throws IllegalArgumentException if a format error in the input data
     * is encountered.
     * @throws IllegalArgumentException if
     * {@code overlap<0 || overlap>=windowSize}.
     * @throws IllegalStateException if
     * {@code this.canAdvanceWindow()==false}.
     */
    void advanceWindow(int overlap, int windowSize);

    /**
     * Returns the window index.  The window index is the number of previous
     * invocations of the {@code advanceWindow()} method.
     * @return the window index
     */
    public int window();

    /**
     * Returns the number of markers overlap between the current marker window
     * and the previous marker window.  Returns 0 if the current marker window
     * is the first marker window.
     *
     * @return the number of markers overlap between the current marker window
     * and the previous marker window.
     */
    public int overlap();

     /**
     * Returns the number of markers in the non-reference data
     * in the overlap between the current marker window and the previous
     * marker window.  Returns 0 if the current marker window is the first
     * marker window.
     *
     * @return the number of markers in the non-reference data
     * in the overlap between the current marker window and the previous
     * marker window.
     */
    int nonRefOverlap();

    /**
     * Returns the number of markers in the union of the current marker window
     * and all previous marker windows.
     * @return the number of markers in the union of the current marker window
     * and all previous marker windows.
     */
    int cumMarkerCnt();

    /**
     * Returns the list of markers in the current marker window.
     * @return the list of markers in the current marker window.
     */
     Markers markers();

    /**
     * Returns the list of markers with non-reference data in the current
     * marker window.
     * @return the list of markers with non-reference data in the current
     * marker window.
     */
    Markers nonRefMarkers();

    /**
     * Returns the number of markers in the current marker window.
     * @return the number of markers in the current marker window.
     */
     int nMarkers();

     /**
     * Returns the number of markers with non-reference data in
     * the current marker window.
     * @return the number of markers with non-reference data in
     * the current marker window.
     */
    int nNonRefMarkers();

     /**
      * Returns the marker index corresponding to the
      * specified marker in the non-reference data.  Indices are with
      * respect to the current marker window.
      * @param nonRefMarker index of a marker in the non-reference data.
      * @return the marker index corresponding to
      * the specified marker in the non-reference data.
      * @throws IndexOutOfBoundsException if
      * {@code nonRefIndex<0 || nonRefIndex>=this.nNonRefMarkers()}.
      */
     int markerIndex(int nonRefMarker);

     /**
      * Returns the marker index in the non-reference data corresponding to
      * the specified marker, or returns -1 if no corresponding marker exists.
      * Indices are with respect to the current marker window.
      * @param refMarker index of a marker in the reference data.
      * @return the marker index in the non-reference data corresponding to
      * the specified marker, or returns -1 if no corresponding marker exists.
      * @throws IndexOutOfBoundsException if
      * {@code refIndex<0 || refIndex>=this.nRefMarkers()}.
      */
     int nonRefMarkerIndex(int refMarker);

    /**
     * Returns the number of reference samples.
     * @return the number of reference samples.
     */
    int nRefSamples();

    /**
     * Returns the list of reference samples, or {@code null} if
     * there are no reference samples.
     * @return the list of reference samples, or {@code null} if
     * there are no reference samples.
     */
    Samples refSamples();

    /**
     * Returns the number of non-reference samples.
     * @return the number of non-reference samples.
     */
    int nNonRefSamples();

    /**
     * Returns the list of non-reference samples.
     * @return the list of non-reference samples.
     */
    Samples nonRefSamples();

    /**
     * Returns the emission probabilities for the reference samples,
     * or returns {@code null} if there are no reference samples.
     * @return the emission probabilities for the reference samples,
     * or returns {@code null} if there are no reference samples.
     */
    GL refEmissions();

    /**
     * Returns the genotype emission probabilities for the
     * non-reference samples at the markers in the non-reference data.
     * @return the genotype emission probabilities for the
     * non-reference samples at the markers in the non-reference data.
     */
    GL nonRefEmissions();

    /**
     * Returns a list of reference haplotype pairs that are restricted
     * to the markers with non-reference data.  The returned list will
     * be empty if there are no reference samples.
     * @return a list of reference haplotype pairs that are restricted
     * to the markers with the non-reference data.
     */
    List<HapPair> restrictedRefHaps();

    /**
     * Returns a list of reference haplotype pairs.  The returned list will
     * be empty if there are no reference samples.
     * @return a list of reference haplotype pairs.
     */
    List<HapPair> refHaps();

    /**
     * Closes any I/O resources controlled by this
     * {@code Data} object.
     */
    void close();
}
