/*
 * Copyright (C) 2015 Brian L. Browning
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
package main;

import vcf.Markers;
import vcf.Marker;
import beagleutil.Samples;
import java.util.Arrays;
import java.util.Comparator;

/**
 * <p>Class {@code BasicAlleleProbs} stores per-haplotype allele probabilities
 * for a list of samples.
 * </p>
 * <p>Instances of class {@code BasicAlleleProbs} are immutable.
 * </p>
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class BasicAlleleProbs implements AlleleProbs {

    private final Markers markers;
    private final Samples samples;
    private final HapAlleleProbs[] alleleProbs;

    /**
     * Construct a new {@code BasicAlleleProbs} instance from the specified
     * data.
     * @param alProbs allele probabilities for each haplotype
     * @throws IllegalArgumentException if
     * {@code alProbs[j].markers().equals(alProbs[k].markers) == false}
     * for any {@code j, k} satisfying
     * {@code (0 <= j && j < k && k < alProbs.length)}
     * @throws IllegalArgumentException if
     * {@code alProbs[j].samples().equals(alProbs[k].samples) == false}
     * for any {@code j, k} satisfying
     * {@code (0 <= j && j < k && k < alProbs.length)}
     * @throws IllegalArgumentException if
     * {@code alProbs.length == 0 || alProbs.length != alProbs[0].nMarkers()}
     * @throws NullPointerException if
     * {@code alProbs == null || alProbs[j] == null} for any {@code j} satisfying
     * {@code (0 <= j && j < alProbs.length)}
     */
    public BasicAlleleProbs(HapAlleleProbs[] alProbs) {
        if (alProbs.length==0) {
            throw new IllegalArgumentException("alProbs.length==0");
        }
        this.alleleProbs = alProbs.clone();
        Arrays.sort(alleleProbs, comparator());
        this.markers = alleleProbs[0].markers();
        this.samples = alleleProbs[0].samples();
        for (int j=1; j<alleleProbs.length; ++j) {
            if (markers.equals(alleleProbs[j].markers())==false) {
                throw new IllegalArgumentException("inconsistent markers");
            }
            if (samples.equals(alleleProbs[j].samples())==false) {
                throw new IllegalArgumentException("inconsistent samples");
            }
        }
    }

    private static Comparator<HapAlleleProbs> comparator() {
        return (HapAlleleProbs t, HapAlleleProbs t1) -> {
            if (t.hapIndex() != t1.hapIndex()) {
                return t.hapIndex() < t1.hapIndex() ? -1 : 1;
            }
            else {
                return 0;
            }
        } ;
    }

    @Override
    public float alProb1(int marker, int sample, int allele) {
        assert alleleProbs[2*sample].hapIndex()/2 == sample;
        return alleleProbs[2*sample].allele(marker, allele);
    }

    @Override
    public float alProb2(int marker, int sample, int allele) {
        assert alleleProbs[2*sample + 1].hapIndex()/2 == sample;
        return alleleProbs[2*sample + 1].allele(marker, allele);
    }

    @Override
    public float gtProb(int marker, int sample, int allele1, int allele2) {
        return alProb1(marker, sample, allele1)*alProb2(marker, sample, allele2);
    }

    @Override
    public int allele1(int marker, int sample) {
        return alleleProbs[2*sample].alleleWithMaxProb(marker);
    }

    @Override
    public int allele2(int marker, int sample) {
        return alleleProbs[2*sample + 1].alleleWithMaxProb(marker);
    }

    @Override
    public int nMarkers() {
        return markers.nMarkers();
    }

    @Override
    public Markers markers() {
        return markers;
    }

    @Override
    public Marker marker(int marker) {
        return markers.marker(marker);
    }

    @Override
    public int nSamples() {
        return samples.nSamples();
    }

    @Override
    public Samples samples() {
        return samples;
    }
}
