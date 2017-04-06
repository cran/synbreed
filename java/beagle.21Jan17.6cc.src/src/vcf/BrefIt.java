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

import beagleutil.ChromIds;
import beagleutil.Samples;
import blbutil.Const;
import blbutil.Filter;
import blbutil.IntArray;
import blbutil.SampleFileIt;
import blbutil.Utilities;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayDeque;
import java.util.Deque;
import java.util.NoSuchElementException;
import java.util.zip.GZIPInputStream;

/**
 * <p>Class {@code BrefIt} represents  an iterator whose {@code next()}
 * method returns an object storing data from a VCF record with phased,
 * non-missing genotypes.
 * </p>
 * <p>Instances of class {@code BrefIt} are not thread-safe.
 * </p>
 * <p>Methods of this class will terminate the Java Virtual Machine with
 * an error message if an I/O error or file format error is detected.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class BrefIt implements SampleFileIt<VcfEmission> {

    private static final int BUFFER_SIZE = 1<<22;
    private static final int SHIFT = 128;
    private static final String[] EMPTY_STRING_ARRAY = new String[0];
    private static final String err = "Error reading file.";

    private final File file;
    private final Filter<Marker> markerFilter;
    private final DataInputStream is;
    private final long initNumber;
    private final String version;
    private final Samples samples;
    private final int nHaps;

    private final Deque<VcfEmission> emBuffer;

    /**
     * Constructs a new {@code BrefIt} instance.
     * @param brefFile a bref file
     *
     * @throws IllegalArgumentException if a format error is detected in a
     * line of the specified bref file
     * @throws NullPointerException if {@code file == null}
     */
    public BrefIt(File brefFile) {
        this(brefFile, Filter.acceptAllFilter());
    }

    /**
     * Constructs a new {@code BrefIt} instance.
     * @param brefFile a bref file
     * @param markerFilter a marker filter or {@code null}
     *
     * @throws IllegalArgumentException if a format error is detected in a
     * line of the specified bref file
     * @throws NullPointerException if {@code file == null}
     */
    public BrefIt(File brefFile, Filter<Marker> markerFilter) {
        if (markerFilter == null) {
            markerFilter = Filter.acceptAllFilter();
        }
        this.file = brefFile;
        this.markerFilter = markerFilter;
        this.is = dataInputStream(brefFile);
        this.initNumber = readInitialNumber(is);
        this.version = readVersion(is);
        this.samples = readSamples(is);
        this.nHaps = 2*samples.nSamples();
        this.emBuffer = new ArrayDeque<>(500);
        fillBuffer();
    }

    private static DataInputStream dataInputStream(File file) {
        InputStream is = null;
        try {
            is = new FileInputStream(file);
            is = new GZIPInputStream(is, BUFFER_SIZE);
        } catch (FileNotFoundException ex) {
            Utilities.exit("File not found: " + file, ex);
        }
        catch (IOException ex) {
            Utilities.exit("Error opening: " + file, ex);
        }
        return new DataInputStream(is);
    }

    private static long readInitialNumber(DataInputStream is) {
        try {
            long initialNumber = is.readInt();
            if (initialNumber != Bref.INITIAL_NUMBER) {
                String s = "ERROR: unrecognized input file.  Was file created "
                        + Const.nl
                        + "with a different version of the bref program?";
                Utilities.exit(s);
            }
            return initialNumber;
        } catch (IOException ex) {
            Utilities.exit(err, ex);
        }
        Utilities.exit(err);
        return -1;
    }

    private static String readVersion(DataInputStream is) {
        try {
            return is.readUTF();
        } catch (IOException ex) {
            Utilities.exit(err, ex);
        }
        Utilities.exit(err);
        return null;
    }

    private static Samples readSamples(DataInputStream is) {
        try {
            int length = is.readInt();
            String[] ids = readStringArray(is, length);
            return Samples.fromIds(ids);
        } catch (IOException ex) {
            Utilities.exit(err, ex);
        }
        Utilities.exit(err);
        return null;
    }

    /**
     * Returns {@code true} if the iteration has more elements, and returns
     * {@code false} otherwise.
     * @return {@code true} if the iteration has more elements
     */
    @Override
    public boolean hasNext() {
        return !emBuffer.isEmpty();
    }

    /**
     * Returns the next element in the iteration.
     * @return the next element in the iteration
     * @throws NoSuchElementException if the iteration has no more elements
     */
    @Override
    public VcfEmission next() {
        if (hasNext()==false) {
            throw new NoSuchElementException();
        }
        VcfEmission first = emBuffer.removeFirst();
        if (emBuffer.isEmpty()) {
            fillBuffer();
        }
        return first;
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
        try {
            is.close();
        } catch (IOException ex) {
            Utilities.exit("Error closing file", ex);
        }
        emBuffer.clear();
    }

    @Override
    public File file() {
        return file;
    }

    @Override
    public Samples samples() {
        return samples;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(80);
        sb.append(this.getClass().toString());
        sb.append(" : ");
        sb.append(file);
        return sb.toString();
    }

    private void fillBuffer() {
        assert emBuffer.isEmpty();
        try {
            int nRecords = is.readInt();
            if (nRecords>0) {
                String chrom  = is.readUTF();
                int chromIndex = ChromIds.instance().getIndex(chrom);
                int nSeq = is.readByte() + SHIFT;
                IntArray hapToSeq = readHapToSeq(nSeq);
                for (int j=0; j<nRecords; ++j) {
                    Marker marker = readMarker(chromIndex);
                    byte flag = is.readByte();
                    switch (flag) {
                        case 0:
                            VcfEmission em = readSeqCodedRecord(marker,
                                    samples, hapToSeq, nSeq);
                            if (markerFilter.accept(marker)) {
                                emBuffer.add(em);
                            }
                            break;
                        case 1:
                            em = readLowMafRecord(marker, samples);
                            if (markerFilter.accept(marker)) {
                                emBuffer.add(em);
                            }
                            break;
                        default:
                            Utilities.exit("Error reading file.");
                    }
                }
            }
        } catch (IOException ex) {
            Utilities.exit("Error reading file", ex);
        }
    }

    private VcfEmission readSeqCodedRecord(Marker marker, Samples samples,
            IntArray hapToSeq, int nSeq) throws IOException {
        IntArray seqToAllele = readSeqToAllele(nSeq, marker.nAlleles());
        return new SeqCodedRefGT(marker, samples, hapToSeq, seqToAllele);
    }

    private IntArray readHapToSeq(int nSeq) throws IOException {
        int[] hap2seq = new int[nHaps];
        for (int j=0; j<hap2seq.length; ++j) {
            hap2seq[j] = is.readByte() + SHIFT;
            if (hap2seq[j] >= nSeq) {
                throw new IllegalStateException("inconsistent data");
            }
        }
        return IntArray.create(hap2seq, 0, (nSeq-1));
    }

    private IntArray readSeqToAllele(int nSeq, int nAlleles) throws IOException {
        int[] seqToAllele = new int[nSeq];
        for (int j=0; j<seqToAllele.length; ++j) {
            seqToAllele[j] = is.readByte() + SHIFT;
            if (seqToAllele[j] >= nAlleles) {
                throw new IllegalStateException("inconsistent data");
            }
        }
        return IntArray.create(seqToAllele, 0, (nAlleles-1));
    }

    private VcfEmission readLowMafRecord(Marker marker, Samples samples)
            throws IOException {
        int nAlleles = marker.nAlleles();
        int[][] hapIndices = new int[nAlleles][];
        for (int j=0; j<nAlleles; ++j) {
            int length = is.readInt();
            hapIndices[j] = (length == -1) ? null : readIntArray(is, length);
        }
        if (nAlleles==2) {
            int x = (hapIndices[0]==null) ? 1 : 0;
            return new LowMafRefDiallelicGT(marker, samples, x, hapIndices[x]);
        }
        else {
            return new LowMafRefGT(marker, samples, hapIndices);
        }
    }

    private Marker readMarker(int chromIndex) throws IOException {
        int end = -1;
        int pos = is.readInt();
        int length = is.readByte() + SHIFT;
        String[] ids = readStringArray(is, length);
        String[] strAlleles;
        byte alleleCode = is.readByte();
        if (alleleCode == -1) {
            length = is.readInt();
            strAlleles = readStringArray(is, length);
            end = is.readInt();
        }
        else {
            int nAlleles = 1 + (alleleCode & 0b11);
            int permIndex = alleleCode >> 2;
            strAlleles = Bref.alleleString(permIndex, nAlleles);
        }
        return new BasicMarker(chromIndex, pos, ids, strAlleles, end);
    }

    private static int[] readIntArray(DataInputStream is, int length)
            throws IOException {
        int[] ia = new int[length];
        for (int j=0; j<ia.length; ++j) {
            ia[j] = is.readInt();
        }
        return ia;
    }

    /* Returns null if length is negative */
    private static String[] readStringArray(DataInputStream is, int length)
            throws IOException {
        if (length < 0) {
            return null;
        }
        else if (length==0) {
            return EMPTY_STRING_ARRAY;
        }
        else {
            String[] sa = new String[length];
            for (int j=0; j<sa.length; ++j) {
                sa[j] = is.readUTF();
            }
            return sa;
        }
    }
}
