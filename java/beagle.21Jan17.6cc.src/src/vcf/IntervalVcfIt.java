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

import blbutil.SampleFileIt;
import beagleutil.ChromInterval;
import beagleutil.Samples;
import java.io.File;
import java.util.NoSuchElementException;

/**
 * <p>Class {@code IntervalVcfIterator} is a sample file iterator whose
 * {@code next()} method returns a marker container.
 * </p>
 *
 * @param <E> the type parameter
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class IntervalVcfIt<E extends MarkerContainer>
        implements SampleFileIt<E> {

    private final SampleFileIt<E> it;
    private final ChromInterval interval;
    private E next;

    /**
     * Constructs a new {@code IntervalVcfIterator} instance.
     * @param it an iterator whose {@code next()} method returns a marker
     * container
     * @param interval a chromosome interval
     * @throws NullPointerException if {@code it == null || interval == null}
     */
    public IntervalVcfIt(SampleFileIt<E> it, ChromInterval interval) {
        if (it==null) {
            throw new IllegalArgumentException("it==null");
        }
        if (interval==null) {
            throw new IllegalArgumentException("interval==null");
        }
        this.it = it;
        this.interval = interval;
        this.next = readFirstRecord(it, interval);
    }

    @Override
    public File file() {
        return it.file();
    }

    @Override
    public Samples samples() {
        return it.samples();
    }

    /**
     * Returns {@code true} if the iteration has more elements.
     * @return {@code true} if the iteration has more elements.
     */
    @Override
    public boolean hasNext() {
        return (next != null);
    }

    /**
     * Returns the next element in the iteration.
     * @return the next element in the iteration.
     * @throws NoSuchElementException if the iteration has no more elements.
     */
    @Override
    public E next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        E current = next;
        this.next = readNextRecord(it, interval);
        return current;
    }

    private E readFirstRecord(SampleFileIt<E> it, ChromInterval interval) {
        E nextRecord = null;
        while (nextRecord==null && it.hasNext()) {
            E candidate = it.next();
            if (inInterval(interval, candidate.marker())) {
                nextRecord = candidate;
            }
        }
        return nextRecord;
    }

    private E readNextRecord(SampleFileIt<E> it, ChromInterval interval) {
        E nextRecord = null;
        if (it.hasNext()) {
            E candidate = it.next();
            if (inInterval(interval, candidate.marker())) {
                nextRecord = candidate;
            }
        }
        return nextRecord;
    }

    private static boolean inInterval(ChromInterval interval, Marker marker) {
        return (marker.chromIndex() == interval.chromIndex()
                && interval.start() <= marker.pos()
                && marker.pos() <= interval.end());
    }

    /**
     * The {@code remove} method is not supported by this iterator.
     * @throws UnsupportedOperationException if this method is invoked
     */
    @Override
    public void remove() {
        throw new UnsupportedOperationException(this.getClass().toString());
    }

    @Override
    public void close() {
        it.close();
    }
}
