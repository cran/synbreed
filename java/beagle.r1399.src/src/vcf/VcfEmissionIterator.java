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
import blbutil.SampleFileIterator;
import java.io.File;
import java.util.NoSuchElementException;
import main.NuclearFamilies;

/**
 * Class {@code VcfEmissionIterator} is an iterator
 * that returns {@code VcfEmission} objects.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class VcfEmissionIterator implements SampleFileIterator<VcfEmission> {

    private final SampleFileIterator<VcfRecord> it;
    private final VcfEmissionFactory emissionFactory;

    private VcfEmissionIterator(SampleFileIterator<VcfRecord> it,
            VcfEmissionFactory emissionFactory) {
        if (it.samples().nSamples()==0) {
            String s = "missing sample data";
            throw new IllegalArgumentException(s);
        }
        if (emissionFactory==null) {
            throw new IllegalArgumentException("emissionFactory==null");
        }
        this.it = it;
        this.emissionFactory = emissionFactory;
    }

    /**
     * Constructs a new {@code SampleFileIterator} instance whose
     * {@code next()} method returns {@code VcfEmission} objects whose
     * emission probabilities are determined by called genotypes.
     * @param it an iterator whose {@code next()} method returns
     * a {@code VcfRecord} objects.
     * @param pedFile a linkage-format pedigree file, or {@code null}
     * if no pedigree relationships are known.  A pedigree file must have
     * at least 4 white-space delimited columns.  The first column of the
     * pedigree file (family ID) is ignored.  The second, third, and fourth
     * columns are the individual ID, father's ID, and mother's ID respectively.
     * @param usePhase {@code true} if phase information in the specified
     * VCF file should be used, and {@code false} otherwise.
     * @return a new {@code SampleFileIterator} instance.
     *
     * @throws IllegalArgumentException if the VCF records contain no samples.
     * @throws IllegalArgumentException if the specified
     * {@code pedFile} is not {@code null} and contains
     * a non-blank line having less than 4 white-space delimited fields
     * or contains duplicate individual identifiers in the second column.
     * @throws NullPointerException if {@code it==null}.
     */
    public static SampleFileIterator<VcfEmission> gt(
            SampleFileIterator<VcfRecord> it, File pedFile, boolean usePhase) {
        NuclearFamilies fam = new NuclearFamilies(it.samples(), pedFile);
        VcfEmissionFactory vef = gt(fam, usePhase);
        return new VcfEmissionIterator(it, vef);
    }

    /**
     * Constructs a new {@code SampleFileIterator} instance whose
     * {@code next()} method returns {@code VcfEmission} objects whose
     * emission probabilities are determined by genotype likelihoods.
     * @param it an iterator whose {@code next()} method returns
     * a {@code VcfRecord} objects. If a VCF record has GL and PL format
     * codes, the PL format field data will be ignored.
     * @param pedFile a linkage-format pedigree file, or {@code null}
     * if no pedigree relationships are known.  A pedigree file must have
     * at least 4 white-space delimited columns.  The first column of the
     * pedigree file (family ID) is ignored.  The second, third, and fourth
     * columns are the individual ID, father's ID, and mother's ID respectively.
     * @param maxLR maximum likelihood ratio.  If the likelihood ratio between
     * two possible genotypes is larger than {@code maxLR}, then the
     * smaller likelihood is set to 0.0, unless the change would create a
     * Mendelian inconsistency in a parent-offspring trio or duo.  In such
     * a case the unmodified likelihoods are used for all members of the
     * inconsistent duo or trio.
     * @return a new {@code SampleFileIterator} instance.
     *
     * @throws IllegalArgumentException if the VCF records contain no samples.
     * @throws IllegalArgumentException if the specified
     * {@code pedFile} is not {@code null} and contains
     * a non-blank line having less than 4 white-space delimited fields
     * or contains duplicate individual identifiers in the second column.
     * @throws IllegalArgumentException if
     * {@code Float.isNaN(maxLR) || maxLR<=1.0f}.
     * @throws NullPointerException if {@code it==null}.
*/
    public static SampleFileIterator<VcfEmission> gl(
            SampleFileIterator<VcfRecord> it, File pedFile, float maxLR) {
        NuclearFamilies fam = new NuclearFamilies(it.samples(), pedFile);
        VcfEmissionFactory vef = gl(fam, maxLR);
        return new VcfEmissionIterator(it, vef);
    }

    /**
     * Constructs a new {@code SampleFileIterator} instance whose
     * {@code next()} method returns {@code VcfEmission} objects whose
     * emission probabilities are determined by called genotypes
     * or genotype likelihoods.
     * @param it an iterator whose {@code next()} method returns
     * a {@code VcfRecord} object.  If a VCF record has GL and PL
     * format codes, the PL format field data will be ignored.
     * @param pedFile a linkage-format pedigree file, or {@code null}
     * if no pedigree relationships are known.  A pedigree file must have
     * at least 4 white-space delimited columns.  The first column of the
     * pedigree file (family ID) is ignored.  The second, third, and fourth
     * columns are the individual ID, father's ID, and mother's ID respectively.
     * @param usePhase {@code true} if phase information in the specified
     * VCF record will be used when a sample's genotype emission probabilities
     * are determined by a called genotype, and {@code false} if
     * phase information in the specified VCF file record will be ignored.
     * @param maxLR maximum likelihood ratio.  If the likelihood ratio between
     * two possible genotypes is larger than {@code maxLR}, then the
     * smaller likelihood is set to 0.0, unless the change would create a
     * Mendelian inconsistency in a parent-offspring trio or duo.  In such
     * a case the unmodified likelihoods are used for all members of the
     * inconsistent duo or trio.
     * @param preferGL {@code true} if genotype emission probabilities
     * should be determined from the VCF record's GL or PL format field if
     * the GL or PL format code is present, and {@code false} if
     * genotype emission probabilities should be determined from the VCF
     * record's a GT format field if the GT format code is present.
     * @return a new {@code SampleFileIterator} instance.
     *
     * @throws IllegalArgumentException if the VCF file contains no samples.
     * @throws IllegalArgumentException if the specified
     * {@code pedFile} is not {@code null} and contains
     * a non-blank line having less than 4 white-space delimited fields
     * or contains duplicate individual identifiers in the second column
     * @throws IllegalArgumentException if
     * {@code Float.isNaN(maxLR) || maxLR<=1.0f}
     * @throws NullPointerException if {@code it==null}
     */
    public static SampleFileIterator<VcfEmission> gtgl(
            SampleFileIterator<VcfRecord> it, File pedFile, boolean usePhase,
            float maxLR, boolean preferGL) {
        NuclearFamilies fam = new NuclearFamilies(it.samples(), pedFile);
        VcfEmissionFactory vef = gtgl(fam, usePhase, maxLR, preferGL);
        return new VcfEmissionIterator(it, vef);
    }

    @Override
    public Samples samples() {
        return it.samples();
    }

    @Override
    public File file() {
        return it.file();
    }

    @Override
    public void close() {
        it.close();
    }

    /**
     * Returns {@code true} if the iteration has more elements.
     * @return {@code true} if the iteration has more elements.
     */
    @Override
    public boolean hasNext() {
        return it.hasNext();
    }

    /**
     * Returns the next element in the iteration.
     * @return the next element in the iteration.
     * @throws NoSuchElementException if the iteration has no more elements.
     */
    @Override
    public VcfEmission next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        return emissionFactory.create(it.next());
    }

    /**
     * The {@code remove} method is not supported by this iterator.
     * @throws UnsupportedOperationException if this method is invoked.
     */
    @Override
    public void remove() {
        throw new UnsupportedOperationException("VcfGtIterator.remove()");
    }

    private static VcfEmissionFactory gt(final NuclearFamilies fam,
            final boolean usePhase) {
        return new VcfEmissionFactory() {
            @Override
            public VcfEmission create(VcfRecord rec) {
                return new BitSetGT(rec, fam, usePhase);
            }
        };
    }

    private static VcfEmissionFactory gl(final NuclearFamilies fam,
            final float maxLR) {
        if (Float.isNaN(maxLR) || maxLR <= 1.0f) {
            throw new IllegalArgumentException("maxLR: " + maxLR);
        }
        return new VcfEmissionFactory() {
            @Override
            public VcfEmission create(VcfRecord rec) {
                return new MedMemGL(rec, fam, maxLR);
            }
        };
    }
    private static VcfEmissionFactory gtgl(final NuclearFamilies fam,
            final boolean usePhase, final float maxLR, final boolean preferGL) {
        if (Float.isNaN(maxLR) || maxLR <= 1.0f) {
            throw new IllegalArgumentException("maxLR: " + maxLR);
        }
        return new VcfEmissionFactory() {
            @Override
            public VcfEmission create(VcfRecord rec) {
                return new MedMemGTGL(rec, fam, usePhase, maxLR, preferGL);
            }
        };
    }
}
