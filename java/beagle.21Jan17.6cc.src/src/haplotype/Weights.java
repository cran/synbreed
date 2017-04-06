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
package haplotype;

import beagleutil.Samples;
import java.util.HashMap;
import java.util.Map;
import main.NuclearFamilies;

/**
 * <p>Class {@code Weights} represents per-haplotype weights.
 * </p>
 * Instances of class {@code Weights} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Weights {

    private final NuclearFamilies fam;
    private final float nonRefWt;

    /**
     * Constructs a new {@code Weights} instance with a weight of 1.0f
     * for all samples.
     * @param fam the parent-offspring relationships
     * @throws NullPointerException if {@code fam == null}
     */
    public Weights(NuclearFamilies fam) {
        this(fam, 1.0f);
    }

    /**
     * Constructs a new {@code Weights} instance with a weight of 1.0f
     * for reference samples, and a weight of {@code nonRefWt} for
     * non-reference samples.  Non-reference samples are samples
     * which are not present in {@code fam.samples()}.
     * @param fam the parent-offspring data
     * @param nonRefWt the non-reference sample weight
     * @throws IllegalArgumentException if
     * {@code nonRefWt < 0.0f || nonRefWt > 1.0f || Float.isNaN(nonRefWt)}
     * @throws NullPointerException if {@code fam == null}
     */
    public Weights(NuclearFamilies fam, float nonRefWt) {
        if (fam==null) {
            throw new NullPointerException("fam==null");
        }
        if (nonRefWt < 0.0f || nonRefWt > 1.0f || Float.isNaN(nonRefWt)) {
            throw new IllegalArgumentException("nonRefWeight: " + nonRefWt);
        }
        this.fam = fam;
        this.nonRefWt = nonRefWt;
    }

   /**
     * Returns an array of length {@code haps.nHaps()} with
     * per-haplotype weights.  Array elements {@code 2*j} and {@code 2*j + 1}
     * are the weights for the first and second haplotype in the
     * {@code j}-th haplotype pair.  Reference haplotypes are assigned
     * a weight of {@code 1.0f}.  Non-reference haplotypes are assigned
     * a weight of {@code this.nonRefWt()} if the haplotype is not
     * inherited from a parent in the sample, and a weight of {@code 0.01f}
     * if the haplotype is inherited from a parent in the sample.
     * The first haplotype in the offspring is required to be the transmitted
     * transmitted haplotype for a parent-offspring duo.
     *
     * @param haps an array of haplotype pairs
     * @return an array of per-haplotype weights
     *
     * @throws NullPointerException if {@code hapPairs == null}
     */
    public float[] get(HapPairs haps) {
        Samples samples = families().samples();
        float[] fa = new float[haps.nHaps()];
        Map<Integer, Integer> cntMap = cntMap(haps);
        int hapIndex = 0;
        for (int j=0, n=haps.nHapPairs(); j<n; ++j) {
            int idIndex = haps.samples(j).idIndex(haps.sampleIndex(j));
            int sampleIndex = samples.index(idIndex);
            int parentCnt = 0;
            if (sampleIndex != -1) {
                // sample is a non-reference sample
                if (families().father(sampleIndex)>=0) {
                    ++parentCnt;
                }
                if (families().mother(sampleIndex)>=0) {
                    ++parentCnt;
                }
            }
            float sampleWeight = (sampleIndex == -1) ? 1.0f : nonRefWt();
            int cnt = cntMap.get(idIndex);
            float wt = sampleWeight/cnt;
            float MIN_SAMPLE_WEIGHT = 0.01f;
            float minWt = MIN_SAMPLE_WEIGHT/cnt;
            fa[hapIndex++] = parentCnt>0  ? minWt : wt;
            fa[hapIndex++] = parentCnt==2 ? minWt : wt;
        }
        return fa;
    }

    /*
     * Returns a map from the haplotype ID index to the number of
     * haplotype pairs with the ID index.
     */
    private static Map<Integer, Integer> cntMap(HapPairs haps) {
        int nHapPairs = haps.nHapPairs();
        int initCapacity = 1 + (3*nHapPairs + 1)/2;
        Map<Integer, Integer> cntMap = new HashMap<>(initCapacity);
        for (int j=0; j<nHapPairs; ++j) {
            int idIndex = haps.samples(j).idIndex(haps.sampleIndex(j));
            Integer value = cntMap.get(idIndex);
            if (value==null) {
                value = 0;
            }
            cntMap.put(idIndex, (value + 1));
        }
        return cntMap;
    }

    /**
     * Returns the parent-offspring relationships.
     * @return the parent-offspring relationships
     */
    public NuclearFamilies families() {
        return fam;
    }

    /**
     * Returns the non-reference sample weight.
     * @return the non-reference sample weight
     */
    public float nonRefWt() {
        return nonRefWt;
    }
}
