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
 * <p>Class {@code HbdAL} represents allele emission probabilities
 * for a set of haplotype pairs under a homozygosity by descent (HBD) model.
 * </p>
 * <p>Instances of class {@code HbdAL} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class HbdAL implements AL {

    private final GL gl;

    /**
     * Constructs an {@code HbdAL} instance.
     *
     * @param gl genotype emission probabilities.
     *
     * @throws NullPointerException if {@code gl==null}.
     */
    public HbdAL(GL gl) {
        if (gl==null) {
            throw new NullPointerException("em==null");
        }
        this.gl = gl;
    }

    @Override
    public float al(int marker, int haplotype, int allele) {
        if (allele<0 || allele >= gl.marker(marker).nAlleles()) {
            String s = "marker=" + marker + " allele: " + allele;
            throw new IllegalArgumentException(s);
        }
        int sample = haplotype/2;
        return gl.gl(marker, sample, allele, allele);
    }

    @Override
    public int allele(int marker, int haplotype) {
        int sample = haplotype/2;
        int a1 = gl.allele1(marker, sample);
        int a2 = gl.allele2(marker, sample);
        return (a1!=a2) ? -1 : a1;
    }

    @Override
    public int nMarkers() {
        return gl.nMarkers();
    }

    @Override
    public Marker marker(int markerIndex) {
        return gl.marker(markerIndex);
    }

    @Override
    public Markers markers() {
        return gl.markers();
    }

    @Override
    public int nSamples() {
        return gl.nSamples();
    }

    @Override
    public Samples samples() {
        return gl.samples();
    }

    @Override
    public int nHaps() {
        return 2*gl.nSamples();
    }

    @Override
    public float errProb() {
        return 0f;
    }

    @Override
    public String toString() {
        StringBuilder sb  = new StringBuilder();
        sb.append("[HbdGL: nMarkers=");
        sb.append(nMarkers());
        sb.append(" nHaps=");
        sb.append(nHaps());
        sb.append(Const.nl);
        sb.append(gl);
        sb.append(Const.nl);
        sb.append(']');
        return sb.toString();
    }
}
