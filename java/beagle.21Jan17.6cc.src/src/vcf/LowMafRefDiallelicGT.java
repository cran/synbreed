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
package vcf;

import beagleutil.Samples;
import java.util.Arrays;

/**
 * <p>Class {@code LowMafRefDiallelicGT} represent represents phased,
 * non-missing genotypes for a list of reference samples at a single diallelic
 * marker.  Genotype emission probabilities are determined by the sample
 * genotypes.
 * </p>
 * <p>
 * Class {@code LowMafRefDiallelicGT} stores the minor allele indices.
 * </p>
 * <p>Instances of class {@code LowMemRefDiallelicGT} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class LowMafRefDiallelicGT implements VcfEmission {

    private final Marker marker;
    private final Samples samples;
    private final int nHaps;
    private final int majorAllele;
    private final int minorAllele;
    private final int[] minorAlleles;

    /**
     * Constructs a new {@code LowMafRefDiallelicGT} instance with phased
     * non-missing genotypes from the specified marker, samples, and
     * haplotype indices.
     * @param marker the marker
     * @param samples the samples
     * @param minorAllele the minor allele
     * @param minorIndices an array whose elements are indices of haplotypes
     * that carry the minor allele
     *
     * @throws IllegalArgumentException {@code marker.nAlleles() != 2}
     * @throws IllegalArgumentException
     * {@code minorAllele < 0 || minorAllele > 1}
     * @throws IllegalArgumentException if any element in
     * {@code minorIndices} is negative or greater than or equal to
     * {@code 2*samples.nSamples()}
     * @throws IllegalArgumentException if any two elements in
     * {@code minorIndices} are equal
     * @throws NullPointerException if
     * {@code marker == null || samples == null || minorIndices == null}
     */
    public LowMafRefDiallelicGT(Marker marker, Samples samples, int minorAllele,
            int[] minorIndices) {
        int[] sortedIndices = checkAndSortIndices(marker, samples, minorAllele,
                minorIndices);
        this.marker = marker;
        this.samples = samples;
        this.nHaps = 2*samples.nSamples();
        this.majorAllele = 1 - minorAllele;
        this.minorAllele = minorAllele;
        this.minorAlleles = sortedIndices;
    }

    private static int[] checkAndSortIndices(Marker marker, Samples samples,
            int minorAllele, int[] minorIndices) {
        if (marker.nAlleles() != 2 || minorAllele < 0 || minorAllele > 1) {
            throw new IllegalArgumentException("ERROR: inconsistent data");
        }
        int[] sorted = minorIndices.clone();
        Arrays.sort(sorted);
        for (int j=1; j<sorted.length; ++j) {
            if (sorted[j] == sorted[j-1]) {
                throw new IllegalArgumentException("ERROR: inconsistent data");
            }
        }
        int nHaps = 2*samples.nSamples();
        if (sorted.length>0
                && (sorted[0]<0 || sorted[sorted.length-1] >= nHaps)) {
            throw new IllegalArgumentException("ERROR: inconsistent data");
        }
        return sorted;
    }

    @Override
    public Samples samples() {
        return samples;
    }

    @Override
    public int nSamples() {
        return samples.nSamples();
    }

    @Override
    public int nHaps() {
        return nHaps;
    }

    @Override
    public int nHapPairs() {
        return samples.nSamples();
    }

    @Override
    public Marker marker() {
        return marker;
    }

    @Override
    public boolean isRefData() {
        return true;
    }

    @Override
    public float gl(int sample, int allele1, int allele2) {
        if (allele1 != 0 &&  allele1 != 1) {
            throw new IndexOutOfBoundsException(String.valueOf(allele1));
        }
        if (allele2 != 0 &&  allele2 != 1) {
            throw new IndexOutOfBoundsException(String.valueOf(allele2));
        }
        boolean matches = (allele1==allele1(sample) && allele2==allele2(sample));
        return matches ? 1.0f : 0.0f;
    }

    @Override
    public boolean isPhased(int sample) {
        if (sample >= samples.nSamples()) {
            throw new IndexOutOfBoundsException(String.valueOf(sample));
        }
        return true;
    }

    @Override
    public int allele1(int sample) {
        return allele(2*sample);
    }

    @Override
    public int allele2(int sample) {
        return allele(2*sample + 1);
    }

    @Override
    public int allele(int hap) {
        if (hap < 0 || hap >= nHaps) {
            throw new IndexOutOfBoundsException(String.valueOf(hap));
        }
        if (Arrays.binarySearch(minorAlleles, hap) >= 0) {
            return minorAllele;
        }
        else {
            return majorAllele;
        }
    }

    @Override
    public int nAlleles() {
        return this.marker().nAlleles();
    }

    @Override
    public boolean storesNonMajorIndices() {
        return true;
    }

    @Override
    public int majorAllele() {
        return majorAllele;
    }

    @Override
    public int alleleCount(int allele) {
        if (allele==majorAllele) {
            throw new IllegalArgumentException("major allele");
        }
        else {
            return minorAlleles.length;
        }
    }

    @Override
    public int hapIndex(int allele, int copy) {
        if (allele==majorAllele) {
            throw new IllegalArgumentException("major allele");
        }
        else {
            return minorAlleles[copy];
        }
    }

    /**
     * Returns the data represented by {@code this} as a VCF
     * record with a GT format field. The returned VCF record
     * will have missing QUAL and INFO fields, will have "PASS"
     * in the filter field, and will have a GT format field.
     * @return the data represented by {@code this} as a VCF
     * record with a GT format field
     */
    @Override
    public String toString() {
        return toVcfRec();
    }
}
