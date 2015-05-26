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

import blbutil.SampleFileIterator;
import beagleutil.ChromInterval;
import beagleutil.Samples;
import java.io.File;
import java.util.NoSuchElementException;

/**
 * Class {@code IntervalVcfIterator} is an iterator whose
 * {@code next()} method returns VCF records contained within a
 * chromosome interval.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class IntervalVcfIterator implements SampleFileIterator<VcfRecord> {

    private final SampleFileIterator<VcfRecord> it;
    private final ChromInterval interval;
    private VcfRecord next;

    /**
     * Constructs a new {@code FilteredVcfIterator} instance.
     * @param it an iterator whose {@code next()} method returns
     * VCF records.
     * @param interval a chromosome interval
     * @throws NullPointerException if {@code it==null || interval==null}.
     */
    public IntervalVcfIterator(SampleFileIterator<VcfRecord> it,
            ChromInterval interval) {
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
    public VcfRecord next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        VcfRecord current = next;
        this.next = readNextRecord(it, interval);
        return current;
    }

    private static VcfRecord readFirstRecord(SampleFileIterator<VcfRecord> it,
            ChromInterval interval) {
        VcfRecord nextRecord = null;
        while (nextRecord==null && it.hasNext()) {
            VcfRecord candidate = it.next();
            if (inInterval(interval, candidate.marker())) {
                nextRecord = candidate;
            }
        }
        return nextRecord;
    }

    private static VcfRecord readNextRecord(SampleFileIterator<VcfRecord> it,
            ChromInterval interval) {
        VcfRecord nextRecord = null;
        if (it.hasNext()) {
            VcfRecord candidate = it.next();
            if (inInterval(interval, candidate.marker())) {
                nextRecord = candidate;
            }
        }
        return nextRecord;
    }

    private static boolean inInterval(ChromInterval interval, Marker marker) {
        return (marker.chromIndex()==interval.chromIndex()
                && interval.start() <= marker.pos()
                && marker.pos() <= interval.end());
    }

    /**
     * The {@code remove} method is not supported by this iterator.
     * @throws UnsupportedOperationException if this method is invoked.
     */
    @Override
    public void remove() {
        throw new UnsupportedOperationException("IntervalVcfIterator.remove()");
    }

    @Override
    public void close() {
        it.close();
    }
}
