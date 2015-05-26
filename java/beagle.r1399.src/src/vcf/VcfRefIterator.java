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
import blbutil.SampleFileIterator;
import java.io.File;
import java.util.NoSuchElementException;

/**
 * Class {@code VcfRefIterator} is an iterator whose {@code next()}
 * method returns {@code VcfEmission} objects.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class VcfRefIterator implements SampleFileIterator<VcfEmission> {

    private final SampleFileIterator<VcfRecord> it;

    /**
     * Constructs a new {@code VcfRefIterator} instance.
     * @param it an iterator whose {@code next()} method returns
     * VCF records.
     * @throws IllegalArgumentException if there are no samples in the input
     * data.
     * @throws NullPointerException if {@code it==null}.
     */
    public VcfRefIterator(SampleFileIterator<VcfRecord> it) {
        if (it.samples().nSamples()==0) {
            String s = "missing sample data";
            throw new IllegalArgumentException(s);
        }
        this.it = it;
    }

    @Override
    public Samples samples() {
        return it.samples();
    }

    @Override
    public File file() {
        return it.file();
    }

    @Override
    public void close() {
        it.close();
    }

    /**
     * Returns {@code true} if the iteration has more elements and
     * returns {@code false} otherwise.
     * @return {@code true} if the iteration has more elements and
     * {@code false} otherwise.
     */
    @Override
    public boolean hasNext() {
        return it.hasNext();
    }

    /**
     * Returns the next element in the iteration.
     * @return the next element in the iteration.
     * @throws NoSuchElementException if the iteration has no more elements.
     *
     * @throws IllegalArgumentException if the next VCF record
     * does not have a "GT" format field, has a missing allele, or
     * has an unphased genotype.
     */
    @Override
    public VcfEmission next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        return new BitSetRefGT(it.next());
    }

    /**
     * The {@code remove} method is not supported by this iterator.
     * @throws UnsupportedOperationException if this method is invoked.
     */
    @Override
    public void remove() {
        throw new UnsupportedOperationException("VcfRefIterator.remove()");
    }
}
