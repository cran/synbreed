/*
 * Copyright 2014 Brian L. Browning
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
package haplotype;

import blbutil.Const;
import blbutil.FileUtil;
import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import vcf.GL;
import vcf.Marker;

/**
 * <p>Class {@code GenotypeCorrection} removes any inconsistencies between
 * haplotype pairs and genotypes that determine genotype likelihoods.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class GenotypeCorrection {

    private static final String headerLine = "MARKER" + Const.tab
                        + "SAMPLE" + Const.tab
                        + "REF" + Const.tab + "ALT" + Const.tab
                        + "INPUT_GT" + Const.tab + "ESTIMATED_GT";

    private GenotypeCorrection() {
        // private constructor to prevent instantiation
    }

    /**
     * Removes any inconsistencies between the specified list of
     * haplotype pairs and the genotypes determined by the {@code allele1()}
     * and {@code allele2()} methods of the specified genotype likelihoods.
     * Inconsistencies are resolved by changing the minimum number
     * of alleles in the haplotype pairs.
     *
     * @param hapPairs a list of haplotype pairs
     * @param gl genotype likelihoods
     * @param seed a seed for generating random numbers
     *
     * @throws IllegalArgumentException if
     * {@code hapPairs.get(j).markers().equals(gl.markers()) == false}
     * for any {@code j} satisfying {@code (0 <= j && j < hapPairs.size())}
     * @throws IllegalArgumentException if
     * {@code hapPairs.get(j).samples().equals(gl.samples()) == false}
     * for any {@code j} satisfying {@code (0 <= j && j < hapPairs.size())}
     *
     * @throws NullPointerException if {@code hapPairs == null || gl == null},
     * or if {@code (hapPair.get(j) == null)} for any {@code j} satisfying
     * {@code (0 <= j && j < hapPairs.size())}
     */
    public static void run(List<HapPair> hapPairs, GL gl, long seed) {
        Random random = new Random(seed);
        int[] alleles1 = new int[gl.nMarkers()];
        int[] alleles2 = new int[gl.nMarkers()];
        for (int j=0, n=hapPairs.size(); j<n; ++j) {
            HapPair hapPair = hapPairs.get(j);
            checkMarkersAndSamples(hapPair, gl);
            List<Edit> edits = getEdits(hapPair, gl, random);
            HapPair revHapPair = updatedHapPair(hapPair, edits, alleles1,
                    alleles2);
            hapPairs.set(j, revHapPair);
        }
    }

    /**
     * Removes any inconsistencies between the specified list of
     * haplotype pairs and the genotypes determined by the {@code allele1()}
     * and {@code allele2()} methods of the specified genotype likelihoods.
     * Inconsistencies are resolved by changing the minimum number
     * of alleles in the haplotype pairs.
     *
     * @param hapPairs a list of haplotype pairs
     * @param gl genotype likelihoods
     * @param seed a seed for generating random numbers
     * @param outFile an output file to which a record of the
     * genotype changes will be written
     * @param append {@code true} if the genotype changes should be
     * written to the end of the specified output file
     *
     * @throws IllegalArgumentException if
     * {@code hapPairs.get(j).markers().equals(gl.markers()) == false}
     * for any {@code j} satisfying {@code (0 <= j && j < hapPairs.size())}
     * @throws IllegalArgumentException if
     * {@code hapPairs.get(j).samples().equals(gl.samples()) == false}
     * for any {@code j} satisfying {@code (0 <= j && j < hapPairs.size())}
     *
     * @throws NullPointerException if
     * {@code (hapPairs == null || gl == null || outFile == null)},
     * or if {@code hapPair.get(j) == null} for any {@code j} satisfying
     * {@code (0 <= j && j < hapPairs.size())}
     */
    public static void run(List<HapPair> hapPairs, GL gl, long seed,
            File outFile, boolean append) {
        Random random = new Random(seed);
        int[] alleles1 = new int[gl.nMarkers()];
        int[] alleles2 = new int[gl.nMarkers()];
        try (PrintWriter out = FileUtil.printWriter(outFile, append)) {
            if (append==false) {
                out.println(headerLine);
            }
            for (int j=0, n=hapPairs.size(); j<n; ++j) {
                HapPair hapPair = hapPairs.get(j);
                checkMarkersAndSamples(hapPair, gl);
                List<Edit> edits = getEdits(hapPair, gl, random);
                for (Edit edit : edits) {
                    out.println(edit);
                }
                HapPair revHapPair = updatedHapPair(hapPair, edits, alleles1,
                        alleles2);
                hapPairs.set(j, revHapPair);
            }
        }
    }

    private static void checkMarkersAndSamples(HapPair hapPair, GL gl) {
        if (hapPair.markers().equals(gl.markers())==false) {
            throw new IllegalArgumentException("inconsistent markers");
        }
        if (hapPair.samples().equals(gl.samples())==false) {
            throw new IllegalArgumentException("inconsistent samples");
        }
    }

    private static List<Edit> getEdits(HapPair hapPair, GL gl, Random random) {
        List<Edit> corrections = new ArrayList<>();
        int sample = hapPair.sampleIndex();
        for (int marker=0, n=gl.nMarkers(); marker<n; ++marker) {
            int hapPairA1 = hapPair.allele1(marker);
            int hapPairA2 = hapPair.allele2(marker);
            if (gl.gl(marker, sample, hapPairA1, hapPairA2) <= 0f) {
                int glA1 = gl.allele1(marker, sample);
                int glA2 = gl.allele2(marker, sample);
                if (glA1>=0 && glA2>=0) {
                    if (gl.isPhased(marker, sample)==false && random.nextBoolean()) {
                        int tmp = glA1;
                        glA1 = glA2;
                        glA2 = tmp;
                    }
                    corrections.add(new Edit(hapPair, marker, glA1, glA2));
                }
            }
        }
        return corrections;
    }

    private static HapPair updatedHapPair(HapPair hapPair,
            List<Edit> edits, int[] alleles1, int[] alleles2) {
        if (edits.isEmpty()) {
            return hapPair;
        }
        else {
            copyAlleles(hapPair, alleles1, alleles2);
            for (int j=0, n=edits.size(); j<n; ++j) {
                Edit edit = edits.get(j);
                alleles1[edit.marker()] = edit.newAllele1();
                alleles2[edit.marker()] = edit.newAllele2();
            }
            return new BitHapPair(hapPair.markers(), hapPair.samples(),
                    hapPair.sampleIndex(), alleles1, alleles2);
        }
    }

    private static void copyAlleles(HapPair hapPair, int[] alleles1,
            int[] alleles2) {
        assert hapPair.nMarkers()==alleles1.length;
        assert alleles1.length==alleles2.length;
        for (int m=0, n=hapPair.nMarkers(); m<n; ++m) {
            alleles1[m] = hapPair.allele1(m);
            alleles2[m] = hapPair.allele2(m);
        }
    }

    private static class Edit {
        private final HapPair hapPair;
        private final int marker;
        private final int newAllele1;
        private final int newAllele2;

        /**
         * Constructs a new {@code Edit} instance.
         *
         * @param hapPair a haplotype pair
         * @param marker the marker index
         * @param newAllele1 the post-edit first allele
         * @param newAllele2 the post-edit second allele
         * @throws NullPointerException if {@code hapPair == null}
         * @throws IndexOutOfBoundsException if
         * {@code marker < 0 || marker > hapPair.nMarkers()}
         * @throws IndexOutOfBoundsException if
         * {@code newAllele1 < 0 || newAllele1 >= hapPair.marker(marker).nAlleles()}
         * @throws IndexOutOfBoundsException if
         * {@code newAllele2 < 0 || newAllele2 >= hapPair.marker(marker).nAlleles()}
         */
        public Edit(HapPair hapPair, int marker, int newAllele1, int newAllele2) {
            if (marker<0 || marker>hapPair.nMarkers()) {
                throw new IndexOutOfBoundsException("marker=" + marker);
            }
            if (newAllele1<0 || newAllele1>=hapPair.marker(marker).nAlleles()) {
                throw new IndexOutOfBoundsException("newAllele1=" + newAllele1);
            }
            if (newAllele2<0 || newAllele2>=hapPair.marker(marker).nAlleles()) {
                throw new IndexOutOfBoundsException("newAllele2=" + newAllele2);
            }
            this.hapPair = hapPair;
            this.marker = marker;
            this.newAllele1 = newAllele1;
            this.newAllele2 = newAllele2;
        }

        public int marker() {
            return marker;
        }

        public HapPair hapPair() {
            return hapPair;
        }

        public int newAllele1() {
            return newAllele1;
        }

        public int newAllele2() {
            return newAllele2;
        }

        /**
         * Returns a string description of {@code this}.  The returned string
         * has five tab-delimited fields: 1) Marker identifier, 2) REF allele,
         * 3) ALT alleles, 4) pre-edit genotype, and 5) post-edit genotype.
         * @return a string description of {@code this}
         */
        @Override
        public String toString() {
            Marker m = hapPair.markers().marker(marker);
            String sampleId = hapPair.samples().id(hapPair.sampleIndex());

            StringBuilder sb = new StringBuilder(100);
            sb.append(m.id());
            sb.append(Const.tab);
            sb.append(sampleId);
            for (int j=0, n=m.nAlleles(); j<n; ++j) {
                sb.append(j<2 ? Const.tab : Const.comma);
                sb.append(m.allele(j));
            }
            sb.append(Const.tab);
            sb.append(newAllele1);
            sb.append(Const.unphasedSep);
            sb.append(newAllele2);
            sb.append(Const.tab);
            sb.append(hapPair.allele1(marker));
            sb.append(Const.unphasedSep);
            sb.append(hapPair.allele2(marker));
            return sb.toString();
        }
    }
}
