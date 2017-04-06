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

import beagleutil.ChromIds;
import vcf.Marker;

/**
 * <p>Class {@code PositionMap} represents a genetic map obtained by
 * multiplying chromosome position by a scale factor.
 * </p>
 * <p>Instances of class {@code PositionMap} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class PositionMap implements GeneticMap {

    private final double scaleFactor;

    /**
     * Returns the scale factor that is multiplied by the chromosome position
     * to obtain the corresponding genetic map position
     * @return the scale factor.
     */
    public double scaleFactor() {
        return scaleFactor;
    }

    /**
     * Constructs a new {@code PositionMap} instance.
     * @param scaleFactor the factor that is multiplied by
     * a base position to obtain the corresponding genetic map position
     * @throws IllegalArgumentException if
     * {@code scaleFactor <= 0d || Double.isFinite(scaleFactor) == false}
     */
    public PositionMap(double scaleFactor) {
        if (Double.isFinite(scaleFactor) == false || scaleFactor <= 0d) {
            throw new IllegalArgumentException(String.valueOf(scaleFactor));
        }
        this.scaleFactor = scaleFactor;
    }

    @Override
    public int basePos(int chrom, double geneticPosition) {
        if (chrom < 0 || chrom >= ChromIds.instance().size()) {
            throw new IndexOutOfBoundsException(String.valueOf(chrom));
        }
        long pos = Math.round(geneticPosition / scaleFactor);
        if (pos > Integer.MAX_VALUE) {
            throw new IllegalArgumentException(String.valueOf(pos));
        }
        return (int) pos;
    }

    @Override
    public double genPos(Marker marker) {
        return scaleFactor*marker.pos();
    }

    @Override
    public double genPos(int chrom, int basePosition) {
        if (chrom < 0 || chrom >= ChromIds.instance().size()) {
            throw new IndexOutOfBoundsException(String.valueOf(chrom));
        }
        return scaleFactor*basePosition;
    }
}
