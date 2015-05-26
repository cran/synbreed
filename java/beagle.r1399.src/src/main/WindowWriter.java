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
import blbutil.Const;
import blbutil.FileUtil;
import blbutil.IntPair;
import haplotype.SampleHapPairs;
import ibd.IbdSegment;
import java.io.File;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import vcf.VcfWriter;

/**
 * Class for writing Beagle VCF and IBD output data in overlapping
 * marker windows.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class WindowWriter {

    private static final DecimalFormat df2 = new DecimalFormat("#.##");

    private boolean isClosed = false;
    private boolean appendIbd = false;

    private final Samples samples;
    private final File vcfOutFile;
    private final File ibdOutFile;
    private final File hbdOutFile;
    private final PrintWriter vcfOut;
    private final Map<IntPair, IbdSegment> ibdBuffer = new HashMap<>();

    /**
     * Constructs a {@code WindowWriter} object.
     * @param samples the sample whose data will be printed.
     * @param outPrefix the output file prefix.
     *
     * @throws IllegalArgumentException if {@code outPrefix.length()==0}
     * @throws NullPointerException if {@code samples==null || outPrefix==null}
     */
    public WindowWriter(Samples samples, String outPrefix) {
        if (samples==null) {
            throw new NullPointerException("samples==null");
        }
        if (outPrefix==null) {
            throw new NullPointerException("outPrefix==null");
        }
        if (outPrefix.length()==0) {
            throw new IllegalArgumentException("outPrefix.length()==0");
        }
        this.samples = samples;
        this.vcfOutFile = new File(outPrefix + ".vcf.gz");
        this.ibdOutFile = new File(outPrefix + ".ibd");
        this.hbdOutFile = new File(outPrefix + ".hbd");
        this.vcfOut = FileUtil.bgzipPrintWriter(vcfOutFile);

        boolean printGT = true;
        boolean printGP = true;
        boolean printGL = false;
        VcfWriter.writeMetaLines(samples.ids(), Main.version,
                printGT, printGP, printGL, vcfOut);
    }

    /**
     * Returns the samples.
     * @return the samples.
     */
    public Samples samples() {
        return samples;
    }


    /**
     * Returns {@code true} if {@code this.close()} method has
     * been previously invoked and returns {@code false} otherwise.
     *
     * @return {@code true} if {@code this.close()} method has
     * been previously invoked and {@code false} otherwise.
     */
    public boolean isClosed() {
        return isClosed;
    }

    /**
     * Closes this {@code WindowWriter} for writing.  Any call to the
     * {@code print()} after invoking {@code close()} will
     * throw an {@code IllegalStateException}.
     */
    public void close() {
        vcfOut.close();
        isClosed = true;
    }

    /**
     * Submits the specified window of data to be printed.  The data in
     * {@code hapPair} and {@code gv} between markers with index
     * {@code lastSplice} (inclusive) and {@code nextSplice} (exclusive)
     * will be printed.  HBD segments and IBD segments which terminate
     * in this region will be printed after being merged any corresponding
     * HBD and IBD segments from the previous window that did not
     * definitely terminate in the previous window.
     * HBD segments and IBD segments which could extend beyond this region
     * will be stored and will be merged with any corresponding
     * HBD and IBD segments in the next marker window.
     *
     * <p>It is the the caller's responsibility to ensure that the ordered
     * haplotype pairs in the overlap between adjacent marker windows
     * are identical for each sample.
     * </p>
     *
     * @param hapPairs the haplotype pairs.
     * @param gv the scaled posterior genotype probabilities.  If
     * {@code gv==null}, no posterior genotype probabilities will be printed.
     * @param ibdMap a map whose keys are pairs of haplotype indices and whose
     * values are lists of IBD segments involving the haplotype pair key.
     * If {@code ibdMap==null}, no IBD segments will be printed.
     * @param lastSplice first index after splice point for the last marker
     * window.
     * @param nextOverlap first index of overlap with the next marker window.
     * @param nextSplice first index after splice point for the next marker
     * window.
     *
     * @throws IllegalStateException if {@code this.isClosed()==true}.
     * @throws NullPointerException if {@code hapPairs==null}
     * @throws IllegalArgumentException if
     * {@code this.samples().equals(haps.samples())==false}
     * @throws IllegalArgumentException if
     * {@code gv!=null && this.samples().equals(gv.samples())==false}
     * @throws IllegalArgumentException if
     * {@code gv!=null && haps.markers().equals(gv.markers())==false}
     * @throws IndexOutOfBoundsException if
     * {@code lastSplice<0 || lastSplice>nextOverlap || nextOverlap>nextSplice
     *           || nextSplice>haps.nMarkers()}
     */
    public void print(SampleHapPairs hapPairs,
            GenotypeValues gv, Map<IntPair, List<IbdSegment>> ibdMap,
            int lastSplice, int nextOverlap, int nextSplice) {
        if (isClosed) {
            throw new IllegalStateException("isClosed()==true");
        }
        checkData(hapPairs, gv, lastSplice, nextOverlap, nextSplice);

        if (gv==null) {
            VcfWriter.appendRecords(hapPairs, lastSplice, nextSplice, vcfOut);
        }
        else {
            VcfWriter.appendRecords(hapPairs, gv, lastSplice, nextSplice, vcfOut);
        }
        vcfOut.flush();

        if (ibdMap!=null) {
            printIbd(ibdMap, lastSplice, nextOverlap, nextSplice, hapPairs.nMarkers());
            if (appendIbd==false) {
                appendIbd = true;
            }
        }
    }

    private void checkData(SampleHapPairs haps, GenotypeValues gv,
            int lastSplice, int nextOverlap, int nextSplice) {
        if (samples.equals(haps.samples())==false
                || (gv!=null && samples.equals(gv.samples())==false)  ) {
            throw new IllegalArgumentException("inconsistent samples");
        }
        if (gv!=null && haps.markers().equals(gv.markers())==false) {
            throw new IllegalArgumentException("inconsistent markers");
        }
        if (lastSplice<0 || lastSplice>nextOverlap || nextOverlap>nextSplice
                || nextSplice>haps.nMarkers()) {
            throw new IndexOutOfBoundsException("index error");
        }
    }

    private void printIbd(Map<IntPair, List<IbdSegment>> ibd, int lastSplice,
            int nextOverlap, int nextSplice, int nMarkers) {
        Map<IntPair, IbdSegment> lastBuffer = new HashMap<>(ibdBuffer);
        ibdBuffer.clear();
        try (PrintWriter ibdOut=FileUtil.printWriter(ibdOutFile, appendIbd);
               PrintWriter hbdOut = FileUtil.printWriter(hbdOutFile, appendIbd)) {
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
}
