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

import vcf.Marker;
import vcf.Markers;

/**
 * <p>Interface {@code GeneticMap} represents a genetic map for one or more
 * chromosomes.
 * </p>
 * <p>Instances of class {@code GeneticMap} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface GeneticMap {

    /**
     * Returns the base position corresponding to the specified genetic map
     * position. If the genetic position is not a map position then the base
     * position is estimated from the nearest genetic map positions using
     * linear interpolation.
     *
     * @param chrom the chromosome index
     * @param geneticPosition the genetic position on the chromosome
     * @return the base position corresponding to the specified genetic map
     * position
     * @throws IllegalArgumentException if the calculated base position
     * exceeds {@code Integer.MAX_VALUE}
     * @throws IllegalArgumentException if this genetic map has no
     * map positions for the specified chromosome
     * @throws IndexOutOfBoundsException if
     * {@code chrom < 0 || chrom >= ChromIds.instance().size()}
     */
    int basePos(int chrom, double geneticPosition);

    /**
     * Returns the genetic map position of the specified marker. The
     * genetic map position is estimated using linear interpolation.
     *
     * @param marker a genetic marker
     * @return the genetic map position of the specified marker
     * @throws IllegalArgumentException if this genetic map has no
     * map positions for the specified chromosome
     * @throws NullPointerException if {@code marker == null}
     */
    double genPos(Marker marker);

    /**
     * Returns the genetic map position of the specified genome coordinate.
     * The genetic map position is estimated using linear interpolation.
     *
     * @param chrom the chromosome index
     * @param basePosition the base coordinate on the chromosome
     * @return the genetic map position of the specified genome coordinate
     * @throws IllegalArgumentException if this genetic map has no
     * map positions for the specified chromosome
     * @throws IndexOutOfBoundsException if
     * {@code chrom < 0 || chrom >= ChromIds.instance().size()}
     */
    double genPos(int chrom, int basePosition);

    /**
     * Returns a string representation of this genetic map. The exact details
     * of the representation are unspecified and subject to change.
     *
     * @return a string representation of this genetic map
     */
    @Override
    String toString();

    /**
     * Returns the an array of length {@code hapPairs.nMarkers()} whose
     * whose {@code j}-th element is the genetic map position
     * of the {@code j}-th marker.
     * @param markers the list of markers
     * @return an array of genetic map positions
     * @throws NullPointerException if {@code markers == null}
     */
    default double[] genPos(Markers markers) {
        double[] genPos = new double[markers.nMarkers()];
        for (int j=0; j<genPos.length; ++j) {
            genPos[j] = this.genPos(markers.marker(j));
        }
        return genPos;
    }

    /**
     * Returns the an array of length {@code hapPairs.nMarkers()} whose
     * whose {@code j}-th element for {@code j > 0} is the
     * probability of recombination between marker {@code j - 1}
     * and marker {@code j}, and whose initial element is {@code 0}.
     * Any inter-marker genetic distances less than {@code 1e-7} cM are
     * increased to {@code 1e-7} cM.
     * @param markers the list of markers
     * @param nHaps the number of haplotypes in the sample
     * @param ne the effective population size
     * @return an array of inter-marker recombination probabilities
     * @throws IllegalArgumentException if {@code nHaps < 1}
     * @throws IllegalArgumentException if {@code ne < 1f}
     * @throws NullPointerException if {@code markers == null}
     */
    default float[] pRecomb(Markers markers, int nHaps, float ne) {
        if (nHaps < 1) {
            throw new IllegalArgumentException(String.valueOf(nHaps));
        }
        if (ne < 1f) {
            throw new IllegalArgumentException(String.valueOf(ne));
        }
        double MIN_CM_DIST = 1e-7;
        int chrom = markers.marker(0).chromIndex();
        float[] pRecomb = new float[markers.nMarkers()];
        double c = -(0.04*ne/nHaps);    // 0.04 = 4/(100 cM/M)
        double lastGenPos = this.genPos(chrom, markers.marker(0).pos());
        pRecomb[0] = 0f;
        for (int j=1; j<pRecomb.length; ++j) {
            double genPos = this.genPos(markers.marker(j));
            double genDist = Math.max(Math.abs(genPos - lastGenPos), MIN_CM_DIST);
            pRecomb[j] = (float) -Math.expm1(c*genDist);
            lastGenPos = genPos;
        }
        return pRecomb;
    }

}
