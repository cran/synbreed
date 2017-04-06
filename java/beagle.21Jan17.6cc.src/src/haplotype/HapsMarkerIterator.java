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
package haplotype;

import vcf.HapsMarker;
import blbutil.FileIt;
import java.io.File;
import java.util.NoSuchElementException;
import vcf.Marker;

/**
 * <p>Class {@code HapsMarkerIterator} represents a file iterator whose
 * {@code next()} method returns {@code HapsMarker} objects.
 * </p>
 * <p>Instances of class {@code HapsMarkerIterator} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class HapsMarkerIterator implements FileIt<HapsMarker> {

    private final HapPairs haps;
    private int nextIndex = 0;

    /**
     * Constructs a new {@code HapsMarkerIterator} instance that iterates
     * through the markers of the specified {@code HapPairs} object.
     *
     * @param haps the haplotype pairs
     * @throws NullPointerException if {@code haps == null}
     */
    public HapsMarkerIterator(HapPairs haps) {
        if (haps==null) {
            throw new NullPointerException("haps==nullt");
        }
        this.haps = haps;
    }

    @Override
    public File file() {
        return null;
    }

    @Override
    public void close() {
        nextIndex = haps.nMarkers();
    }

    /**
     * Returns {@code true} if the iteration has more elements, and returns
     * {@code false} otherwise.
     * @return {@code true} if the iteration has more elements
     */
    @Override
    public boolean hasNext() {
        return nextIndex < haps.nMarkers();
    }

    /**
     * Returns the next element in the iteration.
     * @return the next element in the iteration
     * @throws NoSuchElementException if the iteration has no more elements
     */
    @Override
    public HapsMarker next() {
        if (hasNext() == false) {
            throw new NoSuchElementException("hasNext()==false");
        }
        final int index = nextIndex++;
        return new HapsMarker() {
            @Override
            public int allele(int haplotype) {
                return haps.allele(index, haplotype);
            }

            @Override
            public int allele1(int hapPair) {
                return haps.allele(index, 2 * hapPair);
            }

            @Override
            public int allele2(int hapPair) {
                return haps.allele(index, 2 * hapPair + 1);
            }

            @Override
            public Marker marker() {
                return haps.marker(index);
            }

            @Override
            public int nHaps() {
                return haps.nHaps();
            }

            @Override
            public int nHapPairs() {
                return haps.nHapPairs();
            }
        };
    }

    /**
     * The {@code remove} method is not supported by this iterator.
     * @throws UnsupportedOperationException if this method is invoked
     */
    @Override
    public void remove() {
        throw new UnsupportedOperationException("Not supported.");
    }
}
