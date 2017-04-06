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
 * <p>Class {@code VcfIt} represents  an iterator whose {@code next()}
 * method returns an object storing data from a VCF record.
 * </p>
 * <p>Instances of class {@code VcfIt} are not thread-safe.
 * </p>
 * <p>Methods of this class will terminate the Java Virtual Machine with
 * an error message if an I/O error or file format error is detected.
 * </p>
 * @param <E> the type parameter
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class VcfIt<E extends MarkerContainer> implements SampleFileIt<E> {

    private static final float DEFAULT_MAX_LR = Float.MAX_VALUE;

    private final VcfHeader vcfHeader;
    private final FileIt<String> it;
    private final Function<String, E> mapper;
    private final Filter<Marker> markerFilter;
    private final Thread fileReaderThread;
    private volatile boolean stopFileReadingThread = false;

    private final BlockingQueue<String[]> stringBuffers;
    private final Deque<E> emBuffer;

    /**
     * The default number of VCF records stored in a buffer, which is 1000.
     */
    public static final int DEFAULT_BUFFER_SIZE = 1000;

    /**
     * A function mapping a string VCF record with GT or GL format fields
     * to a {@code VcfRecord} object.
     */
    public static final BiFunction<VcfHeader, String, VcfRecord> toGTGLRec
            = (VcfHeader h, String s) -> VcfRecord.fromGTGL(h, s, DEFAULT_MAX_LR);

    /**
     * A function mapping a string VCF record with GL format fields
     * to a {@code VcfRecord} object.
     */
    public static final BiFunction<VcfHeader, String, VcfRecord> toGLRec
            = (VcfHeader h, String s) -> VcfRecord.fromGL(h, s, DEFAULT_MAX_LR);

    /**
     * A function mapping a string VCF record with GT format fields
     * to a {@code VcfEmission} object.
     */
    public static final BiFunction<VcfHeader, String, VcfEmission> toBitSetGT
            = (VcfHeader h, String s) -> new BitSetGT(h, s);

    /**
     * Create and returns a new {@code VcfIt} instance from the specified
     * objects.
     * @param <R> the type returned by the returned {@code VcfIt}
     * @param strIt an iterator that returns lines of a VCF file
     * @param recMapper a function mapping string VCF records to
     * {@code VcfEmission} objects
     * @return a new {@code VcfIt} instance
     * @throws IllegalArgumentException if a format error is detected in a
     * line of a VCF file returned by {@code strIt}
     * @throws NullPointerException if
     * {@code strIt == null || recMapper == null}
     */
    public static <R extends VcfEmission> VcfIt<R> create(
            FileIt<String> strIt, BiFunction<VcfHeader, String, R> recMapper) {
        return VcfIt.create(strIt, Filter.acceptAllFilter(), recMapper);
    }

    /**
     * Create and returns a new {@code VcfIt} instance from the specified
     * objects.
     * @param <R> the type returned by the returned {@code VcfIt}
     * @param strIt an iterator that returns lines of a VCF file
     * @param sampleFilter a sample filter or {@code null}
     * @param recMapper a function mapping string VCF records to
     * {@code VcfEmission} objects
     * @return a new {@code VcfIt} instance
     * @throws IllegalArgumentException if a format error is detected in a
     * line of a VCF file returned by {@code strIt}
     * @throws NullPointerException if
     * {@code strIt == null || recMapper == null}
     */
    public static <R extends VcfEmission> VcfIt<R> create(
            FileIt<String> strIt, Filter<String> sampleFilter,
            BiFunction<VcfHeader, String, R> recMapper) {
        return VcfIt.create(strIt, sampleFilter, Filter.acceptAllFilter(),
                recMapper);
    }

    /**
     * Create and returns a new {@code VcfIt} instance from the specified
     * objects.
     * @param <R> the type returned by the returned {@code VcfIt}
     * @param strIt an iterator that returns lines of a VCF file
     * @param sampleFilter a sample filter or {@code null}
     * @param markerFilter a marker filter or {@code null}
     * @param recMapper a function mapping string VCF records to
     * {@code VcfEmission} objects
     * @return a new {@code VcfIt} instance
     * @throws IllegalArgumentException if a format error is detected in a
     * line of a VCF file returned by {@code strIt}
     * @throws NullPointerException if
     * {@code strIt == null || recMapper == null}
     */
    public static <R extends VcfEmission> VcfIt<R> create(
            FileIt<String> strIt, Filter<String> sampleFilter,
            Filter<Marker> markerFilter,
            BiFunction<VcfHeader, String, R> recMapper) {
        return VcfIt.create(strIt, sampleFilter, markerFilter, recMapper,
                DEFAULT_BUFFER_SIZE);
    }

    /**
     * Create and returns a new {@code VcfIt} instance from the specified
     * objects.
     * @param <R> the type returned by the returned {@code VcfIt}
     * @param strIt an iterator that returns lines of a VCF file
     * @param sampleFilter a sample filter or {@code null}
     * @param markerFilter a marker filter or {@code null}
     * @param recMapper a function mapping string VCF records to
     * {@code VcfEmission} objects
     * @param bufferSize the buffer size
     * @return a new {@code VcfIt} instance
     * @throws IllegalArgumentException if a format error is detected in a
     * line of a VCF file returned by {@code strIt}
     * @throws IllegalArgumentException if {@code bufferSize < 1}
     * @throws NullPointerException if
     * {@code strIt == null || recMapper == null}
     */
    public static <R extends VcfEmission> VcfIt<R> create(
            FileIt<String> strIt, Filter<String> sampleFilter,
            Filter<Marker> markerFilter,
            BiFunction<VcfHeader, String, R> recMapper, int bufferSize) {
        VcfIt<R> vcfIt = new VcfIt<>(strIt, sampleFilter, markerFilter,
                recMapper, bufferSize);
        vcfIt.start();
        return vcfIt;
    }

    private VcfIt(FileIt<String> it, Filter<String> sampleFilter,
            Filter<Marker> markerFilter,
            BiFunction<VcfHeader, String, E> recMapper, int bufferSize) {
        if (bufferSize < 1) {
            throw new IllegalArgumentException(String.valueOf(bufferSize));
        }
        if (markerFilter==null) {
            markerFilter = Filter.acceptAllFilter();
        }
        this.vcfHeader = new VcfHeader(it, sampleFilter);
        this.it = it;
        this.mapper = (String s) -> recMapper.apply(vcfHeader, s);
        this.markerFilter = markerFilter;
        this.stringBuffers = new ArrayBlockingQueue<>(1);
        this.emBuffer = new ArrayDeque<>(bufferSize);
        this.fileReaderThread = fileReadingThread();
    }

    private void start() {
        this.fileReaderThread.setDaemon(true);
        this.fileReaderThread.start();
        fillEmissionBuffer();
        if (emBuffer.isEmpty()) {
            noRecordFoundError(it);
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
            String line = readLine(it);
            int bufferSize = stringBufferSize(line);
            while (line != null && stopFileReadingThread == false) {
                String chromPlusTab = chromFieldPlusTab(line);
                String[] sa = new String[bufferSize];
                int size = 0;
                while (line != null && size < bufferSize
                        && line.startsWith(chromPlusTab)) {
                    sa[size++] = line;
                    line = readLine(it);
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
        if (bufferSize > DEFAULT_BUFFER_SIZE) {
            bufferSize = DEFAULT_BUFFER_SIZE;
        }
        if (bufferSize < DEFAULT_BUFFER_SIZE/20) {
            bufferSize = DEFAULT_BUFFER_SIZE/20;
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
        while (lastLength != 0 && emBuffer.size() < DEFAULT_BUFFER_SIZE) {
            String[] stringBuffer = takeFromBlockingQueue(stringBuffers);
            lastLength = stringBuffer.length;
            if (stringBuffer.length>0) {
                List<E> list = Arrays.stream(stringBuffer)
                        .parallel()
                        .map(mapper)
                        .filter(e -> markerFilter.accept(e.marker()))
                        .collect(Collectors.toList());
                emBuffer.addAll(list);
            }
            else {
                // put sentinel element back
                putInBlockingQueue(stringBuffers, stringBuffer);
            }
        }
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
        stringBuffers.poll();  // unblock file reading thread
        try {
            fileReaderThread.join();
        } catch (InterruptedException ex) {
            Utilities.exit("Error: InterruptedException", ex);
        }
        it.close();
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
     * @throws NoSuchElementException if the iteration has no more elements.
     */
    @Override
    public E next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        E first = emBuffer.removeFirst();
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
        String s = "remove() is not supported by VcfIterator";
        throw new UnsupportedOperationException(s);
    }

    @Override
    public File file() {
        return it.file();
    }

    @Override
    public Samples samples() {
        return vcfHeader.samples();
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(80);
        sb.append(this.getClass().toString());
        sb.append(" : ");
        sb.append(it.file()==null ? "stdin" : it.file().toString());
        return sb.toString();
    }
}
