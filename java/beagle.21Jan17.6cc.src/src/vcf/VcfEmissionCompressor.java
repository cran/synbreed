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
import blbutil.IntArray;
import blbutil.IntList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * <p>Class {@code VcfEmissionCompressor} compresses a sequence of
 * {@code VcfEmission} objects which contain reference genotype data.
 * Reference genotype data does not contain any unphased or missing genotypes.
 * Compression is performed by storing the list of distinct allele sequences
 * and the allele sequence carried by each haplotype.
 * </p>
 * <p>Class {@code VcfEmissionCompressor} is not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class VcfEmissionCompressor {

    private final Samples samples;
    private final int capacity;
    private final int[] hapToSeq;
    private final List<Marker> markers;
    private final List<IntList> alleleToSeqList;
    private final List<IntList> sequences;
    private final IntList copiedSeqToSrcSeq;

    private int nSeq;

    /**
     * Constructs a new {@code VcfEmissionCompressor} for the specified
     * samples.
     * @param samples the list of samples whose data will be compressed
     * @param capacity the maximum number of allele sequences that is
     * permitted to exist in the list of compressed {@code VcfEmission} objects
     * @throws IllegalArgumentException if {@code capacity < 0}
     * @throws NullPointerException if {@code samples == null}
     */
    public VcfEmissionCompressor(Samples samples, int capacity) {
        if (samples == null) {
            throw new NullPointerException("samples==null");
        }
        if (capacity < 0) {
            throw new IllegalArgumentException(String.valueOf(capacity));
        }
        this.samples = samples;
        this.capacity = capacity;
        this.hapToSeq = new int[2*samples.nSamples()];
        this.markers = new ArrayList<>(100);
        this.sequences = new ArrayList<>(100);
        this.alleleToSeqList = new ArrayList<>(100);
        this.copiedSeqToSrcSeq = new IntList(50);
        clear();
    }

    /**
     * Returns the list of samples whose phased genotype data will be compressed.
     * @return the list of samples whose phased genotype data will be compressed
     */
    public Samples samples() {
        return samples;
    }

    /**
     * Returns the maximum number of allele sequences that is
     * permitted to exist in the list of compressed {@code VcfEmission} objects.
     * @return the maximum number of allele sequences that is
     * permitted to exist in the list of compressed {@code VcfEmission} objects
     */
    public int capacity() {
        return capacity;
    }

    /**
     * Attempts to add the specified {@code VcfEmission} object to the list of
     * compressed {@code VcfEmission} objects, and returns {@code true}
     * if the {@code VcfEmission} object was added.
     *
     * @param em reference genotypes for a marker
     * @return {@code true} if the specified {@code VcfEmission} object was
     * added to the list of compressed markers, and {@code false}
     * others
     * @throws IllegalArgumentException if
     * {@code em.samples().equals(this.samples()) == false}
     * @throws IllegalArgumentException if {@code em.isRefData() == false}
     * @throws NullPointerException if {@code em == null}
     */
    public boolean addToCompessedList(VcfEmission em) {
        checkEmission(em);
        if (inconsistentChrom(em) || em.marker().nAlleles() > capacity) {
            return false;
        }
        boolean success = true;
        int startNSeq = nSeq;
        for (int j=0; j<nSeq; ++j) {
            alleleToSeqList.get(j).clear();
        }
        copiedSeqToSrcSeq.clear();
        for (int h=0; h<hapToSeq.length && success; ++h) {
            success = addHaplotype(em, h);
        }
        if (success) {
            markers.add(em.marker());
        }
        else {
            rollBackChanges(startNSeq);
        }
        return success;
    }

    private void checkEmission(VcfEmission em) {
        if (em.samples().equals(samples)==false) {
            throw new IllegalArgumentException("inconsistent samples");
        }
        if (em.isRefData()==false) {
            throw new IllegalArgumentException("unphased data");
        }
    }

    private boolean inconsistentChrom(VcfEmission em) {
        return (markers.isEmpty()==false
                && em.marker().chromIndex() != markers.get(0).chromIndex());
    }

    private boolean addHaplotype(VcfEmission em, int hap) {
        int seq = hapToSeq[hap];
        int allele = em.allele(hap);
        IntList alleleToSeq = alleleToSeqList.get(seq);
        if (alleleToSeq.isEmpty()) {
            alleleToSeq.add(allele);
            alleleToSeq.add(seq);
            sequences.get(seq).add(allele);
        }
        else {
            int index=0;
            while (index < alleleToSeq.size()
                    && alleleToSeq.get(index)!=allele) {
                index+=2;
            }
            if (index==alleleToSeq.size()) {
                if (nSeq == capacity) {
                    return false;
                }
                else {
                    addCopyOfSequence(seq);
                    sequences.get(nSeq - 1).add(allele);
                    alleleToSeq.add(allele);
                    alleleToSeq.add(nSeq - 1);
                    copiedSeqToSrcSeq.add(seq);
                }
            }
            hapToSeq[hap] = alleleToSeq.get(index+1);
        }
        return true;
    }

    private void rollBackChanges(int startNSeq) {
        alleleToSeqList.subList(startNSeq, nSeq).clear();
        sequences.subList(startNSeq, nSeq).clear();
        nSeq = startNSeq;
        for (int h=0; h<hapToSeq.length; ++h) {
            if (hapToSeq[h] >= startNSeq) {
                hapToSeq[h] = copiedSeqToSrcSeq.get(hapToSeq[h] - startNSeq);
            }
        }
    }

    /**
     * Returns the size of the list of compressed {@code VcfEmission} objects.
     * @return the size of the list of compressed {@code VcfEmission} objects
     */
    public int size() {
        return markers.size();
    }

    /**
     * Returns the number of distinct allele sequences in the list of
     * compressed {@code VcfEmission} objects.
     * @return the number of distinct allele sequences
     */
    public int nSeq() {
        return nSeq;
    }

    /**
     * Returns the specified marker.
     * @param index an index in the list of compressed {@code VcfEmission}
     * objects
     * @return the specified marker
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public Marker marker(int index) {
        return markers.get(index);
    }

    /**
     * Returns an array of length {@code 2*this.samples().nSamples()}, whose
     * {@code j}-th element is the index of the allele sequence carried by the
     * {@code j}-th haplotype in the list of compressed {@code VcfEmission}
     * objects.
     * @return an array mapping haplotype indices to allele sequence indices
     */
    public IntArray hapToSeq() {
        if (nSeq==0) {
            return IntArray.create(hapToSeq, 0, 0);
        }
        else {
            return IntArray.create(hapToSeq, 0, nSeq-1);
        }
    }

    /**
     * Returns an array of length {@code this.nSeq()} whose {@code j}-th
     * element is the allele carried by the {@code j}-th distinct allele
     * sequence at {@code this.marker(index)}.
     * @param index an index in the list of compressed {@code VcfEmission}
     * objects
     * @return an array mapping allele sequence indices to allele indices.
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.size()}
     */
    public IntArray seqToAllele(int index) {
        int maxAllele = markers.get(index).nAlleles() - 1;
        int[] seq2allele = new int[nSeq];
        for (int j=0; j<seq2allele.length; ++j) {
            seq2allele[j] = sequences.get(j).get(index);
        }
        return IntArray.create(seq2allele, 0, maxAllele);
    }

    /**
     * Returns the list of compressed {@code VcfEmission} objects.
     *
     * @return the list of compressed {@code VcfEmission} objects
     */
    public List<VcfEmission> getCompressedList() {
        List<VcfEmission> list = new ArrayList<>(markers.size());
        IntArray hap2seq = hapToSeq();
        for (int j=0, n=markers.size(); j<n; ++j) {
                Marker marker = markers.get(j);
                IntArray seq2allele = seqToAllele(j);
                list.add(new SeqCodedRefGT(marker, samples, hap2seq, seq2allele));
        }
        return list;
    }

    /**
     * Clears the list of compressed {@code VcfEmission} objects.
     */
    public final void clear() {
        Arrays.fill(hapToSeq, 0);
        markers.clear();
        sequences.clear();
        alleleToSeqList.clear();
        copiedSeqToSrcSeq.clear();
        nSeq = 0;
        addEmptySequence();
    }

    private void addEmptySequence() {
        alleleToSeqList.add(new IntList(4));
        sequences.add(new IntList(100));
        ++nSeq;
    }

    private void addCopyOfSequence(int seq) {
        addEmptySequence();
        IntList srcSeq = sequences.get(seq);
        IntList destSeq = sequences.get(nSeq - 1);
        assert destSeq.isEmpty();
        for (int j=0, n=markers.size(); j<n; ++j) {
            destSeq.add(srcSeq.get(j));
        }
    }
}