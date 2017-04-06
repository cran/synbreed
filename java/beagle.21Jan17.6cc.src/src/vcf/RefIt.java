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
import beagleutil.Samples;
import blbutil.Const;
import blbutil.FileIt;
import blbutil.Filter;
import blbutil.Utilities;
import java.io.File;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * <p>Class {@code RefIt} represents  an iterator whose {@code next()}
 * method returns an object storing data from a VCF record with
 * phased, non-missing genotypes.
 * </p>
 * <p>Instances of class {@code RefIt} are not thread-safe.
 * </p>
 * <p>Methods of this class will terminate the Java Virtual Machine with
 * an error message if an I/O error or file format error is detected.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class RefIt implements SampleFileIt<VcfEmission> {

    /**
     * The default number of {@code VcfEmission} objects that are
     * stored in a buffer.
     */
    public static final int DEFAULT_EM_BUFFER_SIZE = 1000;

    private static final int MAX_NSEQ = 255;

    private final VcfHeader vcfHeader;
    private final FileIt<String> strIt;
    private final Function<String, VcfEmission> mapper;
    private final Filter<Marker> markerFilter;
    private final Thread fileReaderThread;
    private volatile boolean stopFileReadingThread = false;

    private final BlockingQueue<String[]> stringBuffers;
    private final Deque<VcfEmission> emBuffer;

    private final List<VcfEmission> uncompressedBuffer;
    private final VcfEmissionCompressor emCompressor;

    public static final BiFunction<VcfHeader, String, VcfEmission> toRef
            = (VcfHeader h, String s) -> refEmission(h, s);

    /**
     * Create and returns a new {@code RefIt} instance from the specified
     * iterator.
     * @param strIt an iterator that returns lines of a VCF file
     * @return a new {@code RefIt} instance
     * @throws IllegalArgumentException if a format error is detected in a
     * line of a VCF file returned by {@code strIt}
     * @throws NullPointerException if {@code strIt == null}
     */
    public static RefIt create(FileIt<String> strIt) {
        return RefIt.create(strIt, Filter.acceptAllFilter(),
                Filter.acceptAllFilter(), DEFAULT_EM_BUFFER_SIZE);
    }

    /**
     * Create and returns a new {@code RefIt} instance from the specified
     * objects.
     * @param strIt an iterator that returns lines of a VCF file
     * @param sampleFilter a sample filter or {@code null}
     * @param markerFilter a marker filter or {@code null}
     * @param bufferSize the buffer size
     * @return a new {@code RefIt} instance
     * @throws IllegalArgumentException if a format error is detected in a
     * line of a VCF file returned by {@code strItt}
     * @throws IllegalArgumentException if {@code bufferSize < 1}
     * @throws NullPointerException if {@code strIt == null}
     */
    public static RefIt create(FileIt<String> strIt,
            Filter<String> sampleFilter, Filter<Marker> markerFilter,
            int bufferSize) {
        RefIt refIt = new RefIt(strIt, sampleFilter, markerFilter, bufferSize);
        refIt.start();
        return refIt;
    }

    private RefIt(FileIt<String> strIt, Filter<String> sampleFilter,
            Filter<Marker> markerFilter, int bufferSize) {
        if (bufferSize < 1) {
            throw new IllegalArgumentException(String.valueOf(bufferSize));
        }
        if (markerFilter==null) {
            markerFilter = Filter.acceptAllFilter();
        }
        this.vcfHeader = new VcfHeader(strIt, sampleFilter);
        this.strIt = strIt;
        this.mapper = (String s) -> toRef.apply(vcfHeader, s);
        this.markerFilter = markerFilter;
        this.stringBuffers = new ArrayBlockingQueue<>(1);
        this.emBuffer = new ArrayDeque<>(DEFAULT_EM_BUFFER_SIZE);
        this.uncompressedBuffer = new ArrayList<>();
        this.emCompressor = new VcfEmissionCompressor(vcfHeader.samples(),
                MAX_NSEQ);
        this.fileReaderThread = fileReadingThread();
    }

    private void start() {
        this.fileReaderThread.setDaemon(true);
        this.fileReaderThread.start();
        fillEmissionBuffer();
        if (emBuffer.isEmpty()) {
            noRecordFoundError(strIt);
        }
    }

    private void noRecordFoundError(FileIt<String> it) {
        if (it.hasNext()==false) {
            StringBuilder sb = new StringBuilder(100);
            sb.append("No VCF records found (data source: ");
            sb.append(it.file()==null ? "stdin" : it.file());
            sb.append(")");
            sb.append(Const.nl);
            sb.append("Check that the chromosome identifiers are the same in each input VCF");
            sb.append(Const.nl);
            sb.append("file and in the \'chrom=\' command line argument (if \'chrom=\' is used).");
            throw new IllegalArgumentException(sb.toString());
        }
    }

    private Thread fileReadingThread() {
        Runnable runnable = () -> {
            String line = readLine(strIt);
            int bufferSize = stringBufferSize(line);
            while (line != null && stopFileReadingThread == false) {
                String chromPlusTab = chromFieldPlusTab(line);
                String[] sa = new String[bufferSize];
                int size = 0;
                while (line != null && size < bufferSize
                        && line.startsWith(chromPlusTab)) {
                    sa[size++] = line;
                    line = readLine(strIt);
                }
                if (size < bufferSize) {
                    sa = Arrays.copyOf(sa, size);
                }
                putInBlockingQueue(stringBuffers, sa);
            }
            if (stopFileReadingThread == false) {
                putInBlockingQueue(stringBuffers, new String[0]);    // sentinel
            }
        };
        return new Thread(runnable);
    }

    private static int stringBufferSize(String line) {
        if (line == null) {
            return 0;
        }
        long nBytesPerLine = 2*line.length();
        Runtime rt = Runtime.getRuntime();
        long maxMem = rt.maxMemory();
        if (maxMem == Long.MAX_VALUE) {
            maxMem = 500 * (1 << 30);
        }
        long bufferSize = maxMem / (100*nBytesPerLine);
        if (bufferSize > DEFAULT_EM_BUFFER_SIZE) {
            bufferSize = DEFAULT_EM_BUFFER_SIZE;
        }
        if (bufferSize < DEFAULT_EM_BUFFER_SIZE/20) {
            bufferSize = DEFAULT_EM_BUFFER_SIZE/20;
        }
        return (int) bufferSize;
    }

    private static <E> void putInBlockingQueue(BlockingQueue<E> q, E e) {
        try {
            q.put(e);
        } catch (InterruptedException ex) {
            Utilities.exit("Error: InterruptedException", ex);
        }
    }

    private static <E> E takeFromBlockingQueue(BlockingQueue<E> q) {
        try {
            return q.take();
        } catch (InterruptedException ex) {
            Utilities.exit("Error: InterruptedException", ex);
        }
        assert false;
        return null;
    }

    private static String chromFieldPlusTab(String vcfRecord) {
        int tabIndex = vcfRecord.indexOf(Const.tab);
        if (tabIndex == -1) {
            String s = Const.nl + "ERROR: Missing tab delimiter in VCV Record:"
                    + Const.nl + vcfRecord
                    + Const.nl + "Exiting Program";
            Utilities.exit(s);
        }
        return vcfRecord.substring(0, tabIndex + 1);
    }

    private void fillEmissionBuffer() {
        assert emBuffer.isEmpty();
        int lastLength = -1;
        while (lastLength != 0 && emBuffer.size() < DEFAULT_EM_BUFFER_SIZE) {
            String[] stringBuffer = takeFromBlockingQueue(stringBuffers);
            lastLength = stringBuffer.length;
            if (stringBuffer.length>0) {
                List<VcfEmission> list = Arrays.stream(stringBuffer)
                        .parallel()
                        .map(mapper)
                        .filter(e -> markerFilter.accept(e.marker()))
                        .collect(Collectors.toList());
                for (int j=0, n=list.size(); j<n; ++j) {
                    VcfEmission e = list.get(j);
                    if (e.storesNonMajorIndices()
                            || e.marker().nAlleles() > MAX_NSEQ) {
                        uncompressedBuffer.add(e);
                    }
                    else {
                        boolean success = emCompressor.addToCompessedList(e);
                        if (success == false) {
                            flushToEmBuffer();
                            success = emCompressor.addToCompessedList(e);
                            assert success;
                        }
                        uncompressedBuffer.add(null);
                    }
                }
            }
            else {
                // put sentinel element back
                putInBlockingQueue(stringBuffers, stringBuffer);
            }
        }
        if (lastLength==0) {
            flushToEmBuffer();
        }
    }

    private void flushToEmBuffer() {
        List<VcfEmission> list = emCompressor.getCompressedList();
        emCompressor.clear();
        int index = 0;
        for (int j=0, n=uncompressedBuffer.size(); j<n; ++j) {
            VcfEmission ve = uncompressedBuffer.get(j);
            if (ve==null) {
                uncompressedBuffer.set(j, list.get(index++));
            }
        }
        emBuffer.addAll(uncompressedBuffer);
        uncompressedBuffer.clear();
    }

    private static String readLine(FileIt<String> it) {
        if (it.hasNext()==false) {
            return null;
        }
        String line = it.next();
        while (line.trim().isEmpty() && it.hasNext()) {
            line = it.next();
        }
        return line;
    }

    @Override
    public void close() {
        stopFileReadingThread = true;
        stringBuffers.poll();   // unblock file reading thread
        try {
            fileReaderThread.join();
        } catch (InterruptedException ex) {
            Utilities.exit("Error: InterruptedException", ex);
        }
        strIt.close();
        emBuffer.clear();
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
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        VcfEmission first = emBuffer.removeFirst();
        if (emBuffer.isEmpty()) {
            fillEmissionBuffer();
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
    public File file() {
        return strIt.file();
    }

    @Override
    public Samples samples() {
        return vcfHeader.samples();
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(80);
        sb.append("RefVcfIt from file: ");
        sb.append(strIt.file()==null ? "stdin" : strIt.file().toString());
        return sb.toString();
    }

    private static VcfEmission refEmission(VcfHeader vcfHeader,
            String vcfRecord) {
        VcfEmission ve = refEmission(new VcfRecGTParser(vcfHeader, vcfRecord));
        int nHaps = 2*ve.nSamples();
        int[] alleleCounts = alleleCounts(ve);
        int nonMajorCnt = nHaps - max(alleleCounts);
        if (nonMajorCnt < (1 + nHaps/200)) {
            if (alleleCounts.length == 2) {
                int minorAllele = 1 - majorAllele(alleleCounts);
                int[] minorIndices = minorIndices(ve, minorAllele, nonMajorCnt);
                return new LowMafRefDiallelicGT(ve.marker(), ve.samples(),
                        minorAllele, minorIndices);
            }
            else {
                int[][] hapIndices = hapIndices(ve, alleleCounts);
                return new LowMafRefGT(ve.marker(), ve.samples(), hapIndices);
            }
        }
        else {
            return ve;
        }
    }

    private static int[] minorIndices(VcfEmission ve, int minorAllele, int mac) {
        int[] minorIndices = new int[mac];
        int index = 0;
        for (int h = 0, n = ve.nHaps(); h < n; ++h) {
            if (ve.allele(h) == minorAllele) {
                minorIndices[index++] = h;
            }
        }
        assert index==mac;
        return minorIndices;
    }

    private static int[][] hapIndices(VcfEmission ve, int[] alCnts) {
        int majorAllele = majorAllele(alCnts);
        int[][]  hapIndices = new int[alCnts.length][];
        for (int j=0; j<hapIndices.length; ++j) {
            hapIndices[j] = (j == majorAllele) ? null : new int[alCnts[j]];
        }
        int[] indices = new int[alCnts.length];
        for (int h=0, n=ve.nHaps(); h<n; ++h) {
            int a = ve.allele(h);
            if (a != majorAllele) {
                hapIndices[a][indices[a]++] = h;
            }
        }
        return hapIndices;
    }

    private static int majorAllele(int[] alleleCnts) {
        int major = 0;
        for (int j=1; j<alleleCnts.length; ++j) {
            if (alleleCnts[j] > alleleCnts[major]) {
                major = j;
            }
        }
        return major;
    }

    private static VcfEmission refEmission(VcfRecGTParser gtp) {
        if (gtp.marker().nAlleles() <= Byte.MAX_VALUE + 1) {
            return new ByteArrayRefGT(gtp);
        }
        else {
            return new BitSetRefGT(gtp);
        }
    }

    private static int max(int[] ia) {
        int maxIndex = 0;
        for (int j=1; j<ia.length; ++j) {
            if (ia[j] > ia[maxIndex]) {
                maxIndex = j;
            }
        }
        return ia[maxIndex];
    }

    private static int[] alleleCounts(VcfEmission ve) {
        int[] cnts = new int[ve.marker().nAlleles()];
        for (int h = 0, n = ve.nHaps(); h < n; ++h) {
            ++cnts[ve.allele(h)];
        }
        return cnts;
    }
}
