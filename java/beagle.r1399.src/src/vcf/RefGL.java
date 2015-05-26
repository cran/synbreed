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
import blbutil.Const;

/**
 * <p>Class {@code RefGL} represents genotype emission probabilities
 * for a reference panel of phased, non-missing genotypes.
 * </p>
 * Instances of class {@code GL} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class RefGL implements GL {

    private final Samples samples;
    private final Markers markers;
    private final VcfEmission[] vma;

    /**
     * Constructs a {@code PhasedEmissions} instance.
     *
     * @param samples the list of samples with phased genotype data.
     * @param vma genotype emission probabilities.  Each array element
     * stores genotype emission probabilities for a single marker.
     * Array elements corresponding to the same chromosome must be
     * contiguous and sorted in chromosome position order.
     *
     * @throws IllegalArgumentException
     * if elements of {@code vma} corresponding to the same chromosome
     * are not contiguous and sorted in chromosome position order, or if any
     * two {@code vma} elements correspond to the same genetic marker.
     * @throws IllegalArgumentException if
     * {@code vma[j].samples().equals(samples)==false}
     * for some {@code 0<j && j<vma.length}.
     * @throws IllegalArgumentException if
     * {@code vma[j].isRefData()==false}
     * for some {@code 0<j && j<vma.length}.
     *
     * @throws NullPointerException if {@code samples==null}.
     * @throws NullPointerException if {@code vma==null}
     * or if any array element is {@code null}.
     */
    public RefGL(Samples samples, VcfEmission[] vma) {
        checkData(samples, vma);
        this.markers = markers(vma);
        this.samples = samples;
        this.vma = vma.clone();
    }

    private static void checkData(Samples samples, VcfEmission[] pma) {
        for (int j=0; j<pma.length; ++j) {
            if (pma[j].samples().equals(samples)==false) {
                String s = "samples=" + samples
                        + Const.nl + "pma[" + j + "].samples()=" + pma[j].samples();
                throw new IllegalArgumentException(s);
            }
            if (pma[j].isRefData()==false) {
                String s = "non-reference data at marker index " + j;
                throw new IllegalArgumentException(s);
            }
        }
    }

    private static Markers markers(VcfEmission[] vma) {
        Marker[] markers = new Marker[vma.length];
        for (int j=0; j<markers.length; ++j) {
            markers[j] = vma[j].marker();
        }
        return new Markers(markers);
    }

    @Override
    public boolean isRefData() {
        return true;
    }

    @Override
    public float gl(int marker, int sample, byte allele1, byte allele2) {
        byte a1 = vma[marker].allele1(sample);
        byte a2 = vma[marker].allele2(sample);
        return (allele1==a1 && allele2==a2) ? 1.0f : 0.0f;
    }

    @Override
    public byte allele1(int marker, int sample) {
        return vma[marker].allele1(sample);
    }

    @Override
    public byte allele2(int marker, int sample) {
        return vma[marker].allele2(sample);
    }

    @Override
    public int nMarkers() {
        return vma.length;
    }

    @Override
    public Marker marker(int markerIndex) {
        return markers.marker(markerIndex);
    }

    @Override
    public Markers markers() {
        return markers;
    }

    @Override
    public int nSamples() {
        return samples.nSamples();
    }

    @Override
    public Samples samples() {
        return samples;
    }

    @Override
    public String toString() {
        StringBuilder sb  = new StringBuilder();
        sb.append("[RefGL: nMarkers=");
        sb.append(nMarkers());
        sb.append(" nSamples=");
        sb.append(nSamples());
        for (VcfEmission vm : vma) {
            sb.append(Const.nl);
            sb.append(vm);
        }
        sb.append(']');
        return sb.toString();
    }
}
