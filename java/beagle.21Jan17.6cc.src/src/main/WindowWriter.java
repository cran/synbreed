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
package main;

import beagleutil.Samples;
import blbutil.BGZIPOutputStream;
import blbutil.Const;
import blbutil.FileUtil;
import blbutil.IntPair;
import blbutil.Utilities;
import ibd.IbdSegment;
import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.Closeable;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import vcf.VcfWriter;

/**
 * <p>Class {@code WindowWriter} writes VCF and IBD output data.
 * </p>
 * <p>Instances of class {@code WindowWriter} are not thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class WindowWriter implements Closeable {

    private static final DecimalFormat df2 = new DecimalFormat("#.##");

    private boolean appendIbd = false;

    private final Samples samples;
    private final String outPrefix;
    private final File vcfOutFile;
    private final File ibdOutFile;
    private final File hbdOutFile;
    private final Map<IntPair, IbdSegment> ibdBuffer = new HashMap<>();

    /**
     * Constructs a new {@code WindowWriter} object.
     * @param samples the sample whose data will be printed
     * @param outPrefix the output file prefix
     *
     * @throws IllegalArgumentException if {@code outPrefix.length() == 0}
     * @throws NullPointerException if
     * {@code samples == null || outPrefix == null}
     */
    public WindowWriter(Samples samples, String outPrefix) {
        if (samples==null) {
            throw new NullPointerException("samples==null");
        }
        if (outPrefix.length()==0) {
            throw new IllegalArgumentException("outPrefix.length()==0");
        }
        this.samples = samples;
        this.outPrefix = outPrefix;
        this.vcfOutFile = new File(outPrefix + ".vcf.gz");
        this.ibdOutFile = new File(outPrefix + ".ibd.gz");
        this.hbdOutFile = new File(outPrefix + ".hbd.gz");

        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        try (PrintWriter vcfOut=new PrintWriter(
                new BGZIPOutputStream(baos, false))) {
            boolean printGT = true;
            boolean printGP = true;
            boolean printGL = false;
            VcfWriter.writeMetaLines(samples.ids(), Main.program,
                    printGT, printGP, printGL, vcfOut);
        }
        boolean append = false;
        writeBytesToFile(baos.toByteArray(), vcfOutFile, append);
    }

    /**
     * Returns the output file prefix.
     * @return the output file prefix
     */
    public String outPrefix() {
        return outPrefix;
    }

    private static void write(byte[] ba, OutputStream out) {
        try {
            out.write(ba);
        } catch (IOException e) {
            Utilities.exit("Error writing byte", e);
        }
    }

    private static void writeBytesToFile(byte[] ba, File file, boolean append) {
        try {
            try (FileOutputStream fos=new FileOutputStream(file, append)) {
                fos.write(ba);
            }
        } catch (IOException e) {
            Utilities.exit("Error writing to file: " + file, e);
        }
    }

    /**
     * Returns the samples whose data is written by {@code this}.
     * @return the samples whose data is written by {@code this}
     */
    public Samples samples() {
        return samples;
    }

    /**
     * Prints VCF records with GT and GP format fields for markers with
     * index between {@code cd.lastSplice()} (inclusive) and
     * {@code cd.nextSplice()} (exclusive).
     *
     * @param cd the input data for the current marker window
     * @param gv scaled genotype probabilities for the target samples
     *
     * @throws NullPointerException if {@code cd == null || gv == null}
     */
    public void printGV(CurrentData cd, GenotypeValues gv) {
        boolean append = true;
        try (PrintWriter vcfOut = FileUtil.bgzipPrintWriter(vcfOutFile, append)) {
            VcfWriter.appendRecords(gv, cd.prevTargetSpliceStart(),
                    cd.nextTargetSpliceStart(), vcfOut);
        }
    }

    /**
     * Prints the data in {@code alProbs} for markers
     * with index between {@code cd.lastSplice()} (inclusive) and
     * {@code cd.nextSplice()} (exclusive) to the output
     * VCF file: {@code this.outPrefix() + ".vcf.gz"}.
     *
     * @param alProbs the estimated haplotype allele probabilities
     * @param isImputed an array of length {@code alProbs.nMarkers()}
     * whose {@code j}-th element is {@code true} if the corresponding
     * marker is imputed, and {@code false} otherwise
     * @param start the starting marker index (inclusive)
     * @param end the ending marker index (exclusive)
     * @param dose {@code true} if the output FORMAT fields should contain
     * a DS subfield, and {@code false} otherwise
     * @param gprobs {@code true} if the output FORMAT fields should contain
     * a GP subfield, and {@code false} otherwise
     * @param nThreads the number of parallel threads to use
     *
     * @throws IllegalArgumentException if
     * {@code isImputed.length != alProbs.nMarkers()}
     * @throws IllegalArgumentException if {@code nThreads < 1}
     * @throws IndexOutOfBoundsException if
     * {@code start < 0 || end > alProbs.nMarkers() || start > end}
     * @throws NullPointerException if
     * {@code alProbs == null || isImputed == null}
     */
    public void print(AlleleProbs alProbs, boolean[] isImputed,
            int start, int end, boolean dose, boolean gprobs, int nThreads) {
        int step = nMarkersPerStep(alProbs.nSamples(), dose, gprobs);
        int nSteps = nSteps(end-start, step);
        final AtomicInteger atomicInt = new AtomicInteger(0);
        final ConcurrentHashMap<Integer, byte[]> map = new ConcurrentHashMap<>();
        ExecutorService es = Executors.newFixedThreadPool(nThreads);
        for (int j=0; j<nThreads; ++j) {
            es.submit(
                () -> {
                    try {
                        int index = atomicInt.getAndIncrement();
                        while (index < nSteps) {
                            int segStart = start + step*index;
                            int segEnd = Math.min(segStart + step, end);
                            ByteArrayOutputStream baos = new ByteArrayOutputStream();
                            try (PrintWriter vcfOut=new PrintWriter(
                                    new BGZIPOutputStream(baos, false))) {
                                VcfWriter.appendRecords(alProbs, isImputed,
                                        segStart, segEnd, dose, gprobs, vcfOut);
                            }
                            map.put(index, baos.toByteArray());
                            index = atomicInt.getAndIncrement();
                        }
                    }
                    catch (Exception ex) {
                        Utilities.exit("", ex);
                    }
                }
            ) ;
        }
        try {
            es.shutdown();
            es.awaitTermination(Long.MAX_VALUE, TimeUnit.DAYS);
        }
        catch (Throwable e) {
            Utilities.exit("ERROR", e);
        }
        print(map, vcfOutFile);
    }

    private static void print(ConcurrentHashMap<Integer, byte[]> map,
            File outFile)  {
        boolean append = true;
        int index = 0;
        try {
            try (OutputStream fos = new BufferedOutputStream(
                    new FileOutputStream(outFile, append))) {
                byte[] bytes = map.get(index++);
                while (bytes!=null) {
                    write(bytes, fos);
                    bytes = map.get(index++);
                }
            }
        } catch (IOException e) {
            Utilities.exit("Error writing to file: " + outFile, e);
        }
    }

    private static int nMarkersPerStep(int nSamples, boolean dose, boolean gprobs) {
        int nBytesPerStep = 65536*50;
        int bytesPerSample = 4;
        if (dose) {
            bytesPerSample += 2;
        }
        if (gprobs) {
            bytesPerSample += 6;
        }
        return nBytesPerStep / (nSamples*bytesPerSample);
    }

    private static int nSteps(int n, int step) {
        int nSteps = n / step;
        if (nSteps * step != n) {
            ++nSteps;
        }
        return nSteps;
    }

    /**
     * Prints IBD segments that end between the markers
     * with index between {@code cd.lastSplice()} (inclusive) and
     * {@code cd.nextSplice()} (exclusive).
     * IBD segments that end on or after the marker with index
     * {@code cd.nextSplice()} are saved so that they can be merged
     * with IBD segments from the next marker window.
     *
     * <p>It is the the caller's responsibility to ensure that the ordered
     * haplotype pairs between adjacent consecutive markers windows
     * are identical for each sample.
     * </p>
     *
     * @param cd the input data for the current window
     * @param ibdMap a map whose keys are pairs of haplotype indices and whose
     * values are lists of IBD segments involving the haplotype pair key
     *
     * @throws IllegalArgumentException if
     * {@code this.samples().equals(cd.targetSamples()) == false}
     * @throws NullPointerException if {@code cd == null || ibdMap == null}
     */
    public void printIbd(CurrentData cd, Map<IntPair, List<IbdSegment>> ibdMap) {
        if (samples.equals(cd.targetSamples()) == false) {
            throw new IllegalArgumentException("inconsistent samples");
        }
        printIbd(ibdMap, cd.prevTargetSpliceStart(), cd.nextTargetOverlapStart(),
                cd.nextTargetSpliceStart(), cd.nTargetMarkers());
        if (appendIbd==false) {
            appendIbd = true;
        }
    }

    private void printIbd(Map<IntPair, List<IbdSegment>> ibd, int lastSplice,
            int nextOverlap, int nextSplice, int nMarkers) {
        Map<IntPair, IbdSegment> lastBuffer = new HashMap<>(ibdBuffer);
        ibdBuffer.clear();
        try (PrintWriter ibdOut = FileUtil.bgzipPrintWriter(ibdOutFile, appendIbd);
               PrintWriter hbdOut = FileUtil.bgzipPrintWriter(hbdOutFile, appendIbd)) {
            Iterator<IntPair> keyIt = ibd.keySet().iterator();
            while (keyIt.hasNext()) {
                IntPair key = keyIt.next();
                List<IbdSegment> list = ibd.get(key);
                for (IbdSegment seg : list) {
                    if (seg.startIndex()==0) {
                        IbdSegment saved = lastBuffer.get(key);
                        if (saved!=null) {
                            seg = merge(saved, seg);
                        }
                    }
                    int ep1 = seg.endIndex()+1;
                    if (ep1>=lastSplice && (nextSplice==nMarkers || ep1<nextSplice)) {
                        printSegment(samples, seg, ibdOut, hbdOut);
                    }
                    else if (seg.startIndex()<nextOverlap) {
                        ibdBuffer.put(key, seg);
                    }
                }
                keyIt.remove();
            }
        }
    }

    private static IbdSegment merge(IbdSegment a, IbdSegment b) {
        assert a.hapPair().equals(b.hapPair());
        assert a.start().chromIndex()==b.start().chromIndex();
        int newStartIndex = -1;
        float newScore = Math.max(a.score(), b.score());
        return new IbdSegment(a.hapPair(), a.start(), b.end(),
                newScore, newStartIndex, b.endIndex());
    }

    private static void printSegment(Samples samples, IbdSegment tract,
            PrintWriter ibdOut, PrintWriter hbdOut) {
        int h1 = tract.hap1();
        int h2 = tract.hap2();
        int s1 = h1/2;
        int s2 = h2/2;
        PrintWriter out = (s1==s2) ? hbdOut : ibdOut;
        out.print(samples.id(s1));
        out.print(Const.tab);
        out.print((h1 % 2) + 1);
        out.print(Const.tab);
        out.print(samples.id(s2));
        out.print(Const.tab);
        out.print((h2 % 2) + 1);
        out.print(Const.tab);
        out.print(tract.start().chrom());
        out.print(Const.tab);
        out.print(tract.start().pos());
        out.print(Const.tab);
        out.print(tract.end().pos());
        out.print(Const.tab);
        out.println(df2.format(tract.score()));
    }

    @Override
    public void close() {
        boolean append = true;
        try {
            try (FileOutputStream fos = new FileOutputStream(vcfOutFile, append);
                    BufferedOutputStream bos = new BufferedOutputStream(fos);
                    BGZIPOutputStream bgzip = new BGZIPOutputStream(bos, true)) {
                // write empty BGZIP block to bgzip by closing bgzip
            }
        } catch (IOException e) {
            Utilities.exit("Error closing file: " + vcfOutFile, e);
        }
    }
}
