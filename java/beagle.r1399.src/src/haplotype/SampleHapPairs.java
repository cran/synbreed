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

import beagleutil.Samples;

/**
 * <p>Interface {@code SampleHapPairs} represents a list of haplotype pairs
 * for which each haplotype pair correspond to a distinct sample.
 * Each pair of haplotypes is guaranteed to have non-missing alleles at each
 * marker.
 * </p>
 * All instances of {@code SampleHaplotypesInterface} are
 * required to be immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface SampleHapPairs extends HapPairs {

    /**
     * Returns the samples.  The {@code k}-th sample corresponds to
     * the {@code k}-th haplotype pair.
     * @return the samples;
     */
    Samples samples();

    /**
     * Returns the number of samples.  The returned value is equal to
     * {@code this.samples().nSamples()}.
     * @return the number of samples.
     */
    int nSamples();
}
