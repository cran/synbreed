/*
 * Copyright (C) 2015 Brian L. Browning
 *
 * This file is part of Beagle
 *
 * Beagle is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either VERSION 3 of the License, or
 * (at your option) any later VERSION.
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

import blbutil.Const;
import blbutil.FileIt;
import blbutil.FileUtil;
import blbutil.Filter;
import blbutil.InputIt;
import blbutil.IntArray;
import blbutil.SampleFileIt;
import blbutil.Utilities;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.zip.GZIPOutputStream;

/**
 * <p>Class {@code Bref} has methods for reading and writing phased,
 * non-missing genotypes that are stored in a "bref" binary VCF file.
 * </p>
 * <p>Instances of class {@code Bref} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Bref {

    private static final String program = "bref.__REV__.jar";
    private static final int BUFFER_SIZE = 1<<22;
    private static final int SHIFT = 128;
    private static final int MAX_NSEQ = 255; // allow nSeq=0 as sentinal
    private static final int STRING_BUFFER_SIZE = 300;

    private static final String[] bases = new String[] {"A", "C", "G", "T"};
    private static final Set<String> basesSet = basesSet();
    private static final String[][] snvPerms = snvPerms();
    private static final Comparator<String[]> allelesComp = allelesComparator();

    /**
     * The initial long in a bref file created with this bref version.
     */
    public static final int INITIAL_NUMBER = 223579146;

    /**
     * The end of file character for a bref file.
     */
    public static final int EOF = 0;

    private final FileIt<String> it;
    private final DataOutputStream os;
    private final VcfHeader vcfHeader;
    private final BlockingQueue<String[]> stringBuffers;
    private final Function<String, VcfEmission> mapper;
    private final List<VcfEmission> emBuffer;
    private final VcfEmissionCompressor emCompressor;

    /**
     * The {@code main()} method is the entry point to the bref program.
     * See the usage() method for usage instructions.
     *
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        if (args.length != 1) {
            System.out.println(usage());
            System.exit(0);
        }
        String fname = args[0];
        if (fname.endsWith(".vcf") || fname.endsWith(".vcf.gz")) {
            FileIt<String> it = inputIterator(fname);
            File brefFile = brefFile(fname);
            Bref bref = new Bref(it, brefFile);
        }
        else if (fname.endsWith(".bref")) {
            File brefFile = new File(fname);
            try (PrintWriter out = FileUtil.stdOutPrintWriter()) {
                writeVcf(brefFile, out);
            }
        }
        else {
            System.out.println(usage());
            System.out.println("Unrecognized filename extension");
            System.exit(0);
        }
    }

    private static FileIt<String> inputIterator(String fname) {
        if (fname.endsWith(".vcf")) {
            return InputIt.fromTextFile(new File(fname));
        }
        else if (fname.endsWith(".vcf.gz")) {
            return InputIt.fromGzipFile(new File(fname));
        }
        else {
            throw new IllegalArgumentException("invalid filename");
        }
    }

    private static File brefFile(String fname) {
        int x = fname.lastIndexOf(".vcf");
        assert x>=0;
        return new File(fname.substring(0, x) + ".bref");

    }

    /**
     * Returns an array that is obtained by taking the first {@code length}
     * elements of the specified permutation of "A", "C", "G", and "T".
     * The list of 24 permutations of "A", "C", "G", and "T" are sorted
     * in lexicographic order.
     * @param permIndex an index of a permutation of the bases "A",
     * "C", "G", and "T"
     * @param length the number of elements in the returned array
     * @return an array that is obtained by taking the first {@code length}
     * elements of the specified permutation of "A", "C", "G", and "T"
     * @throws IndexOutOfBoundsException if
     * {@code permIndex < 0 || permIndex >= 24}
     * @throws IndexOutOfBoundsException if {@code length < 0 || length >= 4}
     */
    public static String[] alleleString(int permIndex, int length) {
        return Arrays.copyOf(snvPerms[permIndex], length);
    }

    /**
     * Constructs a new {@code BinaryDagWriter} which can write a DAG
     * to the specified file.
     *
     * @param it a file iterator that returns lines of a VCF file.
     * @param brefFile filename for the binary reference files that
     * will be read or written.
     * @throws NullPointerException if {@code file==null}.
     */
    private Bref(FileIt<String> it, File brefFile) {
        this.it = it;
        this.os = dataOutputStream(brefFile);
        this.vcfHeader = vcfHeader(it);
        this.mapper = (String s) -> RefIt.toRef.apply(vcfHeader, s);
        this.stringBuffers = new ArrayBlockingQueue<>(1);
        this.emBuffer = new ArrayList<>(500);
        this.emCompressor = new VcfEmissionCompressor(vcfHeader.samples(),
                MAX_NSEQ);
        try {
            writeHeader(vcfHeader.sampleIds(), os);
            startFileReadingThread();
            writeCompressedRecords();
            os.writeInt(EOF);
            os.close();
        } catch (IOException ex) {
            Utilities.exit("Error writing file", ex);
        }
    }

    private static DataOutputStream dataOutputStream(File file) {
        OutputStream os = null;
        try {
            os = new FileOutputStream(file);
            os = new GZIPOutputStream(os, BUFFER_SIZE);
        } catch (FileNotFoundException ex) {
            Utilities.exit("Error opening: " + file, ex);
        } catch (IOException ex) {
            Utilities.exit("IO error: " + file, ex);
        }
        return new DataOutputStream(os);
    }

    private static VcfHeader vcfHeader(FileIt<String> it) {
        Filter<String> sampleFilter = Filter.acceptAllFilter();
        return new VcfHeader(it, sampleFilter);
    }

    private static void writeHeader(String[] sampleIds, DataOutputStream os)
            throws IOException {
        os.writeInt(Bref.INITIAL_NUMBER);
        os.writeUTF(Bref.program);
        os.writeInt(sampleIds.length);
        for (String id : sampleIds) {
            os.writeUTF(id);
        }
    }

    private void startFileReadingThread() {
        Runnable runnable = () -> {
            String line = readLine(it);
            int bufferSize = STRING_BUFFER_SIZE;
            while (line != null) {
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
            putInBlockingQueue(stringBuffers, new String[0]);    // sentinel
        };
        new Thread(runnable).start();
    }

    private static String readLine(FileIt<String> it) {
        if (it.hasNext()==false) {
            return null;
        }
        String line = it.next();
        if (line.trim().isEmpty()) {
            String s = "Blank line in VCF file: "
                    + (it.file()==null ? "stdin" : it.file());
            throw new IllegalArgumentException(s);
        }
        return line;
    }

    private static String chromFieldPlusTab(String vcfRecord) {
        int tabIndex = vcfRecord.indexOf(Const.tab);
        if (tabIndex == -1) {
            String s = "Missing tab delimiter: " + vcfRecord;
            throw new IllegalArgumentException(s);
        }
        return vcfRecord.substring(0, tabIndex + 1);
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

    private void writeCompressedRecords() throws IOException {
        int lastLength = -1;
        while (lastLength != 0) {
            String[] stringBuffer = takeFromBlockingQueue(stringBuffers);
            lastLength = stringBuffer.length;
            if (stringBuffer.length>0) {
                List<VcfEmission> list = convertStringBuffer(stringBuffer);
                for (VcfEmission em : list) {
                    if (em.storesNonMajorIndices()) {
                        emBuffer.add(em);
                    }
                    else if (em.marker().nAlleles() > MAX_NSEQ) {
                        emBuffer.add(new LowMafRefGT(em.marker(), em.samples(),
                                hapIndices(em)));
                    }
                    else {
                        boolean success = emCompressor.addToCompessedList(em);
                        if (success == false) {
                            writeAndClearVcfEmissions();
                            success = emCompressor.addToCompessedList(em);
                            assert success;
                        }
                        emBuffer.add(null);
                    }
                }
            }
        }
        writeAndClearVcfEmissions();
    }

    private List<VcfEmission> convertStringBuffer(String[] stringBuffer) {
        return Arrays.stream(stringBuffer)
                .parallel()
                .map(mapper)
                .collect(Collectors.toList());
    }

    private static int[][] hapIndices(VcfEmission em) {
        int[] alCnts = alleleCounts(em);
        int majorAllele = majorAllele(alCnts);
        int[][]  hapIndices = new int[alCnts.length][];
        for (int j=0; j<hapIndices.length; ++j) {
            hapIndices[j] = (j == majorAllele) ? null : new int[alCnts[j]];
        }
        int[] indices = new int[alCnts.length];
        for (int h=0, n=em.nHaps(); h<n; ++h) {
            int a = em.allele(h);
            if (a != majorAllele) {
                hapIndices[a][indices[a]++] = h;
            }
        }
        return hapIndices;
    }

    private static int[] alleleCounts(VcfEmission em) {
        int[] cnts = new int[em.marker().nAlleles()];
        for (int h = 0, n = em.nHaps(); h < n; ++h) {
            ++cnts[em.allele(h)];
        }
        return cnts;
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

    private void writeAndClearVcfEmissions() throws IOException {
        if (emBuffer.isEmpty()== false) {
            os.writeInt(emBuffer.size());
            os.writeUTF(chrom(emBuffer, emCompressor));
            os.writeByte(emCompressor.nSeq() - SHIFT);
            IntArray hap2seq = emCompressor.hapToSeq();
            for (int j=0, n=hap2seq.size(); j<n; ++j) {
                os.writeByte(hap2seq.get(j) - SHIFT);
            }
            int index = 0;
            for (VcfEmission ve : emBuffer) {
                if (ve==null) {
                    writeCompressedRecord(emCompressor, index++, os);
                }
                else {
                    writeAlleleIndexRecord(ve, os);
                }
            }
            emBuffer.clear();
            emCompressor.clear();
        }
    }

    private static String chrom(List<VcfEmission> emBuffer,
            VcfEmissionCompressor emCompressor) {
        if (emCompressor.size() > 0) {
            return emCompressor.marker(0).chrom();
        }
        else if (emBuffer.size()>0) {
            return emBuffer.get(0).marker().chrom();
        }
        else {
            throw new IllegalArgumentException();
        }
    }

    private static void writeCompressedRecord(VcfEmissionCompressor emCompressor,
            int index, DataOutputStream os) throws IOException {
        Marker marker = emCompressor.marker(index);
        IntArray seq2Allele = emCompressor.seqToAllele(index);
        writeMarker(marker, os);
        byte codingFlag = 0;
        os.writeByte(codingFlag);
        if (marker.nAlleles() <= 256) {
            for (int j=0, n=seq2Allele.size(); j<n; ++j) {
                os.writeByte(seq2Allele.get(j) - SHIFT);
            }
        }
        else {
            for (int j=0, n=seq2Allele.size(); j<n; ++j) {
                os.writeInt(seq2Allele.get(j));
            }
        }
    }

    private static void writeAlleleIndexRecord(VcfEmission ve, DataOutputStream os)
            throws IOException {
        assert ve.storesNonMajorIndices();
        int nAlleles = ve.nAlleles();
        int majorAllele = ve.majorAllele();
        writeMarker(ve.marker(), os);
        byte codingFlag = 1;
        os.writeByte(codingFlag);
        for (int a=0; a<nAlleles; ++a) {
            if (a == majorAllele) {
                os.writeInt(-1);
            }
            else {
                os.writeInt(ve.alleleCount(a));
                for (int c=0; c<ve.alleleCount(a); ++c) {
                    os.writeInt(ve.hapIndex(a, c));
                }
            }
        }
    }

    private static void writeMarker(Marker marker, DataOutputStream os)
            throws IOException {
        os.writeInt(marker.pos());
        int nIds = Math.min(marker.nIds(), 255);
        os.writeByte(nIds - SHIFT);
        for (int j=0; j<nIds; ++j) {
            os.writeUTF(marker.id(j));
        }
        byte alleleCode = isSNV(marker) ? snvCode(marker.alleles()) : -1;
        os.writeByte(alleleCode);
        if (alleleCode == -1) {
            os.writeInt(marker.nAlleles());
            for (int j=0, n=marker.nAlleles(); j<n; ++j) {
                os.writeUTF(marker.allele(j));
            }
            os.writeInt(marker.end());
        }
    }

    private static byte snvCode(String[] alleles) {
        int x = Arrays.binarySearch(snvPerms, alleles, allelesComp);
        if (x < 0) {
            x = (-x - 1);
        }
        int code = (x << 2) + (alleles.length - 1);
        return (byte) code;
    }

    private static boolean isSNV(Marker marker) {
        for (int j=0, n=marker.nAlleles(); j<n; ++j) {
            if (basesSet.contains(marker.allele(j))==false) {
                return false;
            }
        }
        return true;
    }

    private static void writeVcf(File bref, PrintWriter out) {
        try (SampleFileIt<VcfEmission> brefIt = new BrefIt(bref)) {
            if (brefIt.hasNext()) {
                VcfEmission ve = brefIt.next();
                VcfWriter.writeMetaLinesGT(ve.samples().ids(), program, out);
                out.println(ve.toString());
            }
            while (brefIt.hasNext()) {
                out.println(brefIt.next().toString());
            }
        }
    }

    private static Comparator<String[]> allelesComparator() {
        return (String[] o1, String[] o2) -> {
            int n = Math.min(o1.length, o2.length);
            for (int k=0; k<n; ++k) {
                char c1 = o1[k].charAt(0);
                char c2 = o2[k].charAt(0);
                if (c1 != c2) {
                    return (c1 < c2) ? -1 : 1;
                }
            }
            if (o1.length != o2.length) {
                return o1.length < o2.length ? -1 : 1;
            }
            else {
                return 0;
            }
        };
    }

    private static String[][] snvPerms() {
        List<String[]> perms = new ArrayList<>(24);
        permute(new String[0], bases, perms);
        return perms.toArray(new String[0][]);
    }

    private static Set<String> basesSet() {
        Set<String> set = new HashSet<>(4);
        set.addAll(Arrays.asList(bases));
        return Collections.unmodifiableSet(set);
    }

    private static void permute(String[] start, String[] end, List<String[]> perms) {
        if (end.length==0) {
            perms.add(start);
        }
        else {
            for (int j=0; j<end.length; ++j) {
                String[] newStart = Arrays.copyOf(start, start.length + 1);
                newStart[start.length] = end[j];

                String[] newEnd = new String[end.length - 1];
                if (j > 0) {
                    System.arraycopy(end, 0, newEnd, 0, j);
                }
                if (j < newEnd.length) {
                    System.arraycopy(end, j+1, newEnd, j, (newEnd.length - j));
                }
                permute(newStart, newEnd, perms);
            }
        }
    }

    private static String usage() {
        StringBuilder sb = new StringBuilder(500);
        sb.append("usage: java -jar ");
        sb.append(program);
        sb.append(" [vcf]     (creates a .bref file)");
        sb.append(Const.nl);
        sb.append(" or");
        sb.append(Const.nl);
        sb.append("usage: java -jar ");
        sb.append(program);
        sb.append(" [bref]    (prints a .vcf file to standard out)");
        sb.append(Const.nl);
        sb.append(Const.nl);
        sb.append("where");
        sb.append(Const.nl);
        sb.append("  [vcf]  = A vcf file with phased, non-missing genotype data.  If the VCF");
        sb.append(Const.nl);
        sb.append("           file is a text file, its filename should end in \".vcf\".  If the");
        sb.append(Const.nl);
        sb.append("           VCF file is GZIP-compressed, its filename should end in \".vcf.gz\"");
        sb.append(Const.nl);
        sb.append("  [bref] = A binary reference file.  The filename should end in \".bref\"");
        sb.append(Const.nl);
        return sb.toString();
    }
}
