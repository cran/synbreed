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

/**
 * Interface {@code VcfEmissionFactory} is used to create a
 * {@code VcfEmission} instance from a {@code VcfRecord} instances.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface VcfEmissionFactory {

    /**
     * Constructs and returns a new {@code VcfEmission} instance corresponding
     * to the specified VCF record.
     * @param rec a VCF record.
     * @return a new {@code VcfEmission} instance corresponding to the
     * specified VCF record.
     * @throws NullPointerException if {@code rec==null}.
     */
    VcfEmission create(VcfRecord rec);

}
