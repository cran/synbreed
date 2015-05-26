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

import beagleutil.ChromInterval;
import blbutil.SampleFileIterator;
import beagleutil.Samples;
import blbutil.Const;
import blbutil.FileIterator;
import blbutil.Filter;
import blbutil.FilterUtils;
import blbutil.InputIterator;
import java.io.File;
import java.util.HashSet;
import java.util.NoSuchElementException;
import java.util.Set;

/**
 * Class {@code VcfIterator} is an iterator whose {@code next()}
 * method returns VCF records.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class VcfIterator implements SampleFileIterator<VcfRecord> {

    private static final Filter<String> ACCEPT_ALL_FILTER =
            FilterUtils.acceptAllFilter();

    private final VcfHeader vcfHeader;
    private final FileIterator<String> it;
    private final Set<String> chromSet = new HashSet<>();

    private VcfRecord current;
    private VcfRecord next;

    /**
     * Construct and returns a new {@code VcfIterator} instance that
     * reads from standard input.
     * @param vcfFile a file in VCF format.  The VCF records for
     * each chromosome must be contiguous and sorted in order of
     * increasing position.
     * @param sampleFilter a sample filter, or {@code null} if there
     * is no sample filter.
     * @param markerFilter a marker filter, or {@code null} if there
     * is no marker filter.
     * @param chromInterval the the chromosome interval to read, or
     * {@code null} if there is no interval restriction.
     * @return a new {@code SampleFileIterator<VcfRecord>} instance.
     *
     * @throws IllegalArgumentException if any VCF header line
     * does not conform to the VCF specification, or if the first
     * VCF record does not conform to the VCF specification.
     * @throws NullPointerException if {@code vcfFile==null}.
     */
    public static SampleFileIterator<VcfRecord> filteredIterator(File vcfFile,
            Filter<String> sampleFilter, Filter<Marker> markerFilter,
            ChromInterval chromInterval) {
        SampleFileIterator<VcfRecord> it =
                new VcfIterator(vcfFile, sampleFilter);
        if (chromInterval != null) {
            it = new IntervalVcfIterator(it, chromInterval);
        }
        if (markerFilter != null) {
            it = new FilteredVcfIterator(it, markerFilter);
        }
        return it;
    }

    /**
     * Constructs a new {@code VcfIteratorr} instance.
     *
     * @param vcfFile a file in VCF format.  The VCF records for
     * each chromosome must be contiguous and sorted in order of
     * increasing position.
     *
     * @throws IllegalArgumentException if any VCF header line
     * does not conform to the VCF specification, or if the first
     * VCF record does not conform to the VCF specification.
     * @throws NullPointerException if {@code vcfFile==null}.
     */
    public VcfIterator(File vcfFile) {
        this(vcfFile, ACCEPT_ALL_FILTER);
    }

    /**
     * Constructs a new {@code VcfIterator} instance.
     *
     * @param vcfFile a file in VCF format.  The VCF records for
     * each chromosome must be contiguous and sorted in order of
     * increasing position.
     * @param sampleFilter a sample filter, or {@code null} if there
     * is no sample filter.
     *
     * @throws IllegalArgumentException if any VCF header line
     * does not conform to the VCF specification, or if the first
     * VCF record does not conform to the VCF specification.
     * @throws NullPointerException if {@code vcfFile==null}.
     */
    public VcfIterator(File vcfFile, Filter<String> sampleFilter) {
        this(InputIterator.fromGzipFile(vcfFile), sampleFilter);
    }

    /**
     * Constructs a new {@code VcfIterator} instance.
     *
     * @param it an iterator whose {@code next()} method returns
     * lines of a file in VCF format.  The VCF records for each
     * chromosome must be contiguous and sorted in order of increasing
     * position.
     *
     * @throws IllegalArgumentException if any VCF header line
     * does not conform to the VCF specification, or if the first
     * VCF record does not conform to the VCF specification.
     * @throws NullPointerException if {@code it==null}.
     */
    public VcfIterator(FileIterator<String> it) {
        this(it, ACCEPT_ALL_FILTER);
    }

    /**
     * Constructs a new {@code VcfIterator} instance.
     *
     * @param it a {@code FileIterator<String>} whose
     * {@code next()} method returns lines of a file in
     * VCF format.  The VCF records for each chromosome must be
     * contiguous and sorted in order of increasing position.
     * @param sampleFilter a sample filter, or {@code null} if there
     * is no sample filter.
     *
     * @throws IllegalArgumentException if any VCF header line
     * does not conform to the VCF specification, or if the first
     * VCF record does not conform to the VCF specification.
     * @throws NullPointerException if {@code it==null}.
     */
    private VcfIterator(FileIterator<String> it, Filter<String> sampleFilter) {
        if (sampleFilter==null) {
            sampleFilter = ACCEPT_ALL_FILTER;
        }
        this.it = it;
        this.vcfHeader = new VcfHeader(it, sampleFilter);
        this.current = null;
        this.next = readData();
    }

    /**
     * Construct and returns a new {@code VcfIterator} instance that
     * reads from standard input.
     *
     * @return a {@code VcfIterator} instance that reads from standard
     * input.
     * @throws IllegalArgumentException if any VCF header line
     * does not conform to the VCF specification, or if the first
     * VCF record does not conform to the VCF specification.
     */
    public static VcfIterator fromStdin() {
        FileIterator<String> it = new InputIterator(System.in);
        return new VcfIterator(it, ACCEPT_ALL_FILTER);
    }

    /**
     * Construct and returns a new {@code VcfIterator} instance that
     * reads from standard input.
     * @param sampleFilter a sample filter, or {@code null} if there
     * is no sample filter.
     * @return a {@code VcfIterator} instance that reads from standard
     * input.
     *
     * @throws IllegalArgumentException if any VCF header line
     * does not conform to the VCF specification, or if the first
     * VCF record does not conform to the VCF specification.
     */
    public static VcfIterator fromStdin(Filter<String> sampleFilter) {
        FileIterator<String> it = new InputIterator(System.in);
        return new VcfIterator(it, sampleFilter);
    }

    @Override
    public void close() {
        it.close();
        next = null;
    }

    /**
     * Returns {@code true} if the iteration has more elements, and returns
     * {@code false} otherwise.
     * @return {@code true} if the iteration has more elements, and
     * {@code false} otherwise.
     */
    @Override
    public boolean hasNext() {
        return (next != null);
    }

    /**
     * Returns the next element in the iteration.
     * @return the next element in the iteration.
     * @throws NoSuchElementException if the iteration has no more elements.
     */
    @Override
    public VcfRecord next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        current = next;
        next = readData();
        checkMarkerPosOrder(current, next);
        return current;
    }

    /**
     * The {@code remove} method is not supported by this iterator.
     * @throws UnsupportedOperationException if this method is invoked.
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

    private VcfRecord readData() {
        VcfRecord vcfRecord = null;
        while (it.hasNext() && vcfRecord==null) {
            String line = it.next().trim();
            if (line.isEmpty()==false) {
                vcfRecord = new VcfRecord(line, vcfHeader);
            }
        }
        return vcfRecord;
    }

    private void checkMarkerPosOrder(VcfRecord current, VcfRecord next) {
        if (next!=null) {
            Marker m1 = current.marker();
            Marker m2 = next.marker();
            if (m1.chromIndex()==m2.chromIndex()) {
                if (m1.pos() > m2.pos()) {
                    String s = "["
                            + (it.file()==null ? "stdin" : it.file().toString())
                            + "] markers not in chromosomal order: "
                            + Const.nl + m1 + Const.nl + m2;
                    throw new IllegalArgumentException(s);
                }
            }
            else {
                boolean newChrom = chromSet.add(m2.chrom());
                if (newChrom == false) {
                    String s = "["
                            + (it.file()==null ? "stdin" : it.file().toString())
                            + "] non-contiguous markers for chromosome "
                            + m2.chrom();
                    throw new IllegalArgumentException(s);
                }
            }
        }
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("file=");
        sb.append(it.file());
        sb.append(" next=");
        sb.append(next);
        return sb.toString();
    }
}
