/*
 * Copyright 2013 Brian L. Browning
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package vcf;

import beagleutil.Samples;

/**
 * <p>Class {@code FuzzyGL} is a wrapper for a {@code GL}
 * instance that incorporates a fixed error rate for the
 * observed (emitted) allele to differ from the true allele.  Allele
 * errors are independent.
 * </p>
 * Instances of class {@code FuzzyGL} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class FuzzyGL implements GL {

    private final float ee;
    private final float ef;
    private final float ff;

    private final GL gl;

    /**
     * Constructs a {@code FuzzyGL} instance.
     * @param gl the genotype likelihoods without error
     * @param err the allele error rate
     * @throws IllegalArgumentException if
     * {@code Float.isNaN(err) || err < 0 || err >= 1.0}
     * @throws NullPointerException if {@code gl == null}
     */
    public FuzzyGL(GL gl, float err) {
        if (gl==null) {
            throw new NullPointerException("gl==null");
        }
        if (Double.isNaN(err) || err < 0.0 || err >= 1.0) {
            throw new IllegalArgumentException("err: " + err);
        }
        float e = err;
        float f = 1.0f - err;
        this.ee = e*e;
        this.ef = e*f;
        this.ff = f*f;
        this.gl = gl;
    }

    @Override
    public float gl(int marker, int sample, int a1, int a2) {
        // following algorithm is for both diallelic and multi-allelic markers
        int obs1 = gl.allele1(marker, sample);
        int obs2 = gl.allele2(marker, sample);
        if (obs1>=0 && obs2>=0) {
            if (obs1==obs2 || gl.isPhased(marker, sample)) {
                return phasedGL(obs1, obs2, a1, a2);
            }
            else {
                return phasedGL(obs1, obs2, a1, a2)
                        + phasedGL(obs2, obs1, a1, a2);
            }
        }
        else {
            return gl.gl(marker, sample, a1, a2);
        }
    }

    private float phasedGL(int obs1, int obs2, int a1, int a2) {
        if (obs1==a1) {
            return obs2==a2 ? ff : ef;
        }
        else {
            return obs2==a2 ? ef : ee;
        }
    }

    @Override
    public boolean isRefData() {
        return gl.isRefData();
    }

    @Override
    public boolean isPhased(int marker, int sample) {
        return gl.isPhased(marker, sample);
    }

    @Override
    public int allele1(int marker, int sample) {
        return gl.allele1(marker, sample);
    }

    @Override
    public int allele2(int marker, int sample) {
        return gl.allele2(marker, sample);
    }

    @Override
    public int allele(int marker, int hap) {
        return gl.allele(marker, hap);
    }

    @Override
    public int nMarkers() {
        return gl.nMarkers();
    }

    @Override
    public Marker marker(int marker) {
        return gl.marker(marker);
    }

    @Override
    public Markers markers() {
        return gl.markers();
    }

    @Override
    public int nHaps() {
        return gl.nHaps();
    }

    @Override
    public int nSamples() {
        return gl.nSamples();
    }

    @Override
    public Samples samples() {
        return gl.samples();
    }
}
