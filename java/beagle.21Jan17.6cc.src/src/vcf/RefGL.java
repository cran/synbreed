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
 * Instances of class {@code RefGL} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class RefGL implements GL {

    private final Samples samples;
    private final Markers markers;
    private final VcfEmission[] vea;

    /**
     * Constructs a {@code RefGL} instance. Each element of the
     * specified array stores genotype emission probabilities for a single
     * marker. Array elements corresponding to the same chromosome must be
     * contiguous and sorted in chromosome position order.
     *
     * @param samples the list of samples with phased genotype data.
     * @param vea genotype emission probabilities.
     *
     * @throws IllegalArgumentException
     * if elements of {@code vea} corresponding to the same chromosome
     * are not contiguous and sorted in chromosome position order
     * @throws IllegalArgumentException if
     * {@code vea[j].marker().equals(vea[k].marker() == true}
     * for any {@code j, k} satisfying {@code 0 <= j && j < k && k < vea.length}
     * @throws IllegalArgumentException if
     * {@code vea[j].samples().equals(samples) == false}
     * for any {@code j} satisfying {@code 0 <= j && j < vea.length}
     * @throws IllegalArgumentException if
     * {@code vea[j].isRefData() == false} for any {@code j} satisfying
     * {@code 0 <= j && j < vea.length}
     *
     * @throws NullPointerException if {@code samples == null}
     * @throws NullPointerException if {@code vea == null}
     * @throws NullPointerException if {@code vea[j] == null} for any
     * {@code j} satisfying {@code 0 <= j && j < vea.length}
     */
    public RefGL(Samples samples, VcfEmission[] vea) {
        checkData(samples, vea);
        this.markers = markers(vea);
        this.samples = samples;
        this.vea = vea.clone();
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

    private static Markers markers(VcfEmission[] vea) {
        Marker[] markers = new Marker[vea.length];
        for (int j=0; j<markers.length; ++j) {
            markers[j] = vea[j].marker();
        }
        return Markers.create(markers);
    }

    @Override
    public boolean isRefData() {
        return true;
    }

    @Override
    public float gl(int marker, int sample, int allele1, int allele2) {
        int a1 = vea[marker].allele1(sample);
        int a2 = vea[marker].allele2(sample);
        return (allele1==a1 && allele2==a2) ? 1.0f : 0.0f;
    }

    @Override
    public boolean isPhased(int marker, int sample) {
        return true;
    }

    @Override
    public int allele1(int marker, int sample) {
        return vea[marker].allele1(sample);
    }

    @Override
    public int allele2(int marker, int sample) {
        return vea[marker].allele2(sample);
    }

    @Override
    public int allele(int marker, int hap) {
        return vea[marker].allele(hap);
    }

    @Override
    public int nMarkers() {
        return vea.length;
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
    public int nHaps() {
        return 2*samples.nSamples();
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
        sb.append('[');
        sb.append(this.getClass().toString());
        sb.append(": nMarkers=");
        sb.append(nMarkers());
        sb.append(" nSamples=");
        sb.append(nSamples());
        for (VcfEmission vm : vea) {
            sb.append(Const.nl);
            sb.append(vm);
        }
        sb.append(']');
        return sb.toString();
    }
}
