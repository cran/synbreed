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
 * <p>HClass {@code BasicGL} represents genotype emission probabilities
 * for a set of samples.
 * </p>
 * Instances of class {@code GL} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class BasicGL implements GL {

    private final Samples samples;
    private final Markers markers;
    private final VcfEmission[] vma;
    private final boolean isRefData;

    /**
     * Returns the genotype index corresponding to the
     * specified unordered alleles.
     * @param a1 the first allele index of an unordered genotype.
     * @param a2 the second allele index of an unordered genotype.
     * @return the genotype index corresponding to the
     * specified unordered alleles.
     * @throws IllegalArgumentException if {@code a1 < 0 || a2 < 0}.
     */
    public static int genotype(byte a1, byte a2) {
        if (a1<=a2) {
            if (a1 < 0) {
                String s = "allele < 0: " + a1 + " " + a2;
                throw new IllegalArgumentException(s);
            }
            return (a2*(a2+1))/2 + a1;
        }
        else {
            if (a2<0) {
                String s = "allele < 0: " + a1 + " " + a2;
                throw new IllegalArgumentException(s);
            }
            return (a1*(a1+1))/2 + a2;
        }
    }

    /**
     * Constructs a {@code BasicGL} instance.
     *
     * @param samples the list of samples with genotype data.
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
     * {@code vma[j].samples().equals(samples)==false} for some
     * {@code 0<=j && j<vma.length}.
     *
     * @throws NullPointerException if {@code samples==null}.
     * @throws NullPointerException if {@code vma==null}
     * or if any array element is {@code null}.
     */
    public BasicGL(Samples samples, VcfEmission[] vma) {
        checkSamples(samples, vma);
        this.markers = markers(vma);
        this.samples = samples;
        this.vma = vma.clone();
        this.isRefData = isRefData(vma);
    }

    private static void checkSamples(Samples samples, VcfEmission[] mla) {
        for (int j=0; j<mla.length; ++j) {
            if (mla[j].samples().equals(samples)==false) {
                throw new IllegalArgumentException("inconsistent samples");
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

    private static boolean isRefData(VcfEmission[] vma) {
        boolean isRefData = true;
        for (int j=0; j<vma.length && isRefData==true; ++j) {
            if (vma[j].isRefData()==false) {
                isRefData = false;
            }
        }
        return isRefData;
    }

    @Override
    public boolean isRefData() {
        return isRefData;
    }

    @Override
    public float gl(int marker, int sample, byte allele1, byte allele2) {
        return vma[marker].gl(sample, allele1, allele2);
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
        sb.append("[BasicGL: nMarkers=");
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
