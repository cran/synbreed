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
 * <p>Class {@code LowMafRefGT} represent represents phased, non-missing
 * genotypes for a list of reference samples at a single marker.
 * Genotype emission probabilities are determined by the sample
 * genotypes.
 * </p>
 * <p>
 * Class {@code LowMafRefGT} stores the non-major allele indices.
 * </p>
 * <p>Instances of class {@code LowMemRefGT} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class LowMafRefGT implements VcfEmission {

    private final Marker marker;
    private final Samples samples;
    private final int nHaps;
    private final int[][] hapIndices;
    private final int majorAllele;

    /**
     * Constructs a new {@code LowMafRefGT} instance with phased,
     * non-missing genotypes from the specified marker, samples, and haplotype
     * indices. If a haplotype index is duplicated in the specified
     * {@code hapIndices} array, the haplotype will be assigned the allele
     * with the smallest index.
     *
     * @param marker the marker
     * @param samples the samples
     * @param hapIndices an array whose {@code j}-th element is {@code null}
     * if {@code j} is the unique (or first) major allele, and is an array of
     * indices of haplotypes that carry the {@code j}-th allele otherwise
     *
     * @throws IllegalArgumentException if {@code hapIndices[j] == null} and
     * {@code j} is not the unique or first major allele
     * @throws IllegalArgumentException if any haplotype index in
     * {@code hapIndices} is negative or greater than or equal to
     * {@code 2*samples.nSamples()}
     * @throws IllegalArgumentException if
     * {@code marker.nAlleles() != hapIndices.length}
     * @throws NullPointerException if any parameter is {@code null}
     */
    public LowMafRefGT(Marker marker, Samples samples, int[][] hapIndices) {
        int[][] sortedCopy = copyAndSortIndices(hapIndices);
        checkSortedIndices(marker, samples, sortedCopy);
        this.marker = marker;
        this.samples = samples;
        this.nHaps = 2*samples.nSamples();
        this.hapIndices = sortedCopy;
        int nullIndex = -1;
        for (int j=0; j<hapIndices.length; ++j) {
            if (sortedCopy[j] == null) {
                nullIndex = j;
                break;
            }
        }
        assert nullIndex != -1;
        this.majorAllele = nullIndex;
    }

    private static int[][] copyAndSortIndices(int[][] hapIndices) {
        int[][] sortedCopy = new int[hapIndices.length][];
        for (int j=0; j<hapIndices.length; ++j) {
            if (hapIndices[j] != null) {
                sortedCopy[j] = hapIndices[j].clone();
                Arrays.sort(sortedCopy[j]);
            }
        }
        return sortedCopy;
    }

    private static void checkSortedIndices(Marker marker, Samples samples,
            int[][] hapIndices) {
        if (marker.nAlleles() != hapIndices.length) {
            throw new IllegalArgumentException("ERROR: inconsistent data");
        }
        int nHaps = 2*samples.nSamples();
        checkAlleleCounts(hapIndices, nHaps);
        for (int[] ia : hapIndices) {
            if (ia != null) {
                if (ia.length > 0 && (ia[0] < 0 || ia[ia.length-1] >= nHaps)) {
                    throw new IndexOutOfBoundsException(Arrays.toString(ia));
                }
            }
        }
    }

    private static void checkAlleleCounts(int[][] hapIndices, int nHaps) {
        int nMajorAlleles = nHaps;
        int maxIndex = -1;
        int nullIndex = -1;
        for (int j=0; j<hapIndices.length; ++j) {
            if (hapIndices[j] == null) {
                if (nullIndex != -1) {
                    throw new IllegalArgumentException("ERROR: major allele error");
                }
                nullIndex = j;
            }
            else {
                if (maxIndex == -1
                        || hapIndices[j].length > hapIndices[maxIndex].length) {
                    maxIndex = j;
                }
                nMajorAlleles -= hapIndices[j].length;
            }
        }
        boolean majorAlleleError =
                maxIndex != -1
                && hapIndices[maxIndex].length == nMajorAlleles
                && maxIndex < nullIndex;
        if (nullIndex == -1 || majorAlleleError) {
            throw new IllegalArgumentException("ERROR: major allele error");
        }
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
        if (allele1 < 0 || allele1 >= hapIndices.length) {
            throw new IndexOutOfBoundsException(String.valueOf(allele1));
        }
        if (allele2 < 0 || allele2 >= hapIndices.length) {
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
        for (int j=0; j<hapIndices.length; ++j) {
            if (j != majorAllele) {
                if (Arrays.binarySearch(hapIndices[j], hap) >= 0) {
                    return j;
                }
            }
        }
        return majorAllele;
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
        if (hapIndices[allele]==null) {
            throw new IllegalArgumentException("major allele");
        }
        else {
            return hapIndices[allele].length;
        }
    }

    @Override
    public int hapIndex(int allele, int copy) {
        if (hapIndices[allele]==null) {
            throw new IllegalArgumentException("major allele");
        }
        else {
            return hapIndices[allele][copy];
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
