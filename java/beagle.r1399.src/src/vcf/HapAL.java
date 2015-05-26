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
import haplotype.SampleHapPairs;

/**
 * <p>Class {@code HapAL} represents allele emission probabilities
 * for a set of haplotypes.
 * </p>
 * <p>Instances of class {@code HapAL} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class HapAL implements AL {

    private final Markers markers;
    private final SampleHapPairs partialHapPairs;
    private final int[] markers2HapMarkers;
    private final float errProb;
    private final float noErrProb;

    /**
     * Constructs a new {@code HapAL} instance.  Allele emission probabilities
     * are determined by the called alleles at a marker and the specified
     * error probability if the marker is present in the sample haplotype pairs.
     * Allele emission probabilities are constant if a marker is absent in
     * the sample haplotype pairs.
     *
     * @param markers a list of markers.
     * @param partialHapPairs a list of sample haplotype pairs whose markers
     * are a sublist of the specified markers.
     * @param err the probability that the emitted allele differs from the
     * haplotype's actual allele in the case where the marker is present in
     * the sample haplotype pairs.
     * @throws IllegalArgumentException if {@code partialHapsPairs.markers()}
     * is not a sublist of {@code markers}.
     * @throws IllegalArgumentException if
     * {@code Float.isNaN(err) || err<0.0f || err>1.0f}.
     * @throws NullPointerException if
     * {@code markers==null || partialHapPairs==null}.
     */
    public HapAL(Markers markers, SampleHapPairs partialHapPairs,
            float err) {
        if (Float.isNaN(err) || err < 0.0 || err > 1.0) {
            throw new IllegalArgumentException("err: " + err);
        }
        this.markers2HapMarkers = markers2HapMarkers(markers, partialHapPairs);
        this.markers = markers;
        this.partialHapPairs = partialHapPairs;
        this.errProb = err;
        this.noErrProb = 1.0f - err;
    }

    private static int[] markers2HapMarkers(Markers markers,
            SampleHapPairs hapPairs) {
        int[] markers2HapMarkers = new int[markers.nMarkers()];
        int hapIndex = 0;
        for (int j=0; j<markers2HapMarkers.length; ++j) {
            Marker m = markers.marker(j);
            if (hapIndex<hapPairs.nMarkers()
                    && m.equals(hapPairs.marker(hapIndex))) {
                markers2HapMarkers[j] = hapIndex++;
            }
            else {
                markers2HapMarkers[j] = -1;
            }
        }
        if (hapIndex != hapPairs.nMarkers()) {
            String s = "partialHapPairs.markers() is not a sublist of markers";
            throw new IllegalArgumentException(s);
        }
        return markers2HapMarkers;
    }

    @Override
    public float al(int marker, int haplotype, byte allele) {
        if (allele<0 || allele >= markers.marker(marker).nAlleles()) {
            String s = "marker=" + marker + " allele: " + allele;
            throw new IllegalArgumentException(s);
        }
        int index = markers2HapMarkers[marker];
        if (index==-1) {
            return 1.0f;
        }
        else {
            byte a = partialHapPairs.allele(index, haplotype);
            return (a==allele) ? noErrProb : errProb;
        }
    }

    @Override
    public int nMarkers() {
        return markers.nMarkers();
    }

    @Override
    public Marker marker(int markerIndex) {
        return markers.marker(markerIndex);
    }

    @Override
    public Markers markers() {
        return markers;
    }

    @Override
    public int nSamples() {
        return partialHapPairs.nSamples();
    }

    @Override
    public Samples samples() {
        return partialHapPairs.samples();
    }

    @Override
    public String toString() {
        StringBuilder sb  = new StringBuilder();
        sb.append("[HapAL: nMarkers=");
        sb.append(nMarkers());
        sb.append(" nSamples=");
        sb.append(nSamples());
        sb.append(']');
        return sb.toString();
    }
}
