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
package sample;

import dag.Dag;
import haplotype.HapPair;
import java.util.List;
import vcf.GL;

/**
 * <p>Interface {@code SingleBaumInterface} has methods for sampling
 * haplotype pairs.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface SingleBaumInterface {

    /**
     * Returns the directed acyclic graph that determines the transition
     * probabilities.
     * @return the directed acyclic graph that determines the transition
     * probabilities
     */
    Dag dag();

    /**
     * Returns the emission probabilities.
     * @return the emission probabilities
     */
    GL gl();

    /**
     * Returns the number of haplotype pairs that are sampled for each
     * individual.
     * @return the number of haplotype pairs that are sampled for each
     * individual
     */
    int nSamplesPerIndividual();

    /**
     * Returns the initial random seed.
     * @return the initial random seed
     */
    long seed();

    /**
     * <p>Returns a list of {@code this.nSamplesPerIndividual()} sampled
     * haplotype pairs for the specified individual. Haplotype pairs are
     * sampled conditional on the HMM with transition probabilities
     * determined by {@code this.dag()} and emission probabilities
     * determined by {@code this.gl()}.
     * </p>
     * <p>The contract for this method is unspecified if no haplotype pair
     * is consistent with the HMM.
     * </p>
     * @param sample a sample index
     * @return a list of {@code this.nSamplesPerIndividual()} sampled
     * haplotype pairs for the specified individual
     *
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.gl().nSamples()}
     */
    List<HapPair> randomSample(int sample);

    /**
     * <p>Returns a list of {@code this.nSamplesPerIndividual()} sampled
     * haplotype pairs for the specified individual. Haplotype pairs are
     * sampled conditional on the HMM with transition probabilities determined
     * by {@code this.dag()} and emission probabilities determined by
     * {@code this.gl()}. Posterior genotype probabilities are written to
     * the specified array. The posterior probability of the {@code j}-th
     * genotype for the {@code k}-th marker is stored at index
     * {@code gl.markers().sumGenotypes(k) + j} in the {@code gtProbs} array.
     * </p>
     * <p>The contract for this method is unspecified if no haplotype pair
     * is consistent with the HMM.
     * </p>
     * @param sample the sample index
     * @param gtProbs a array to which posterior genotype probabilities
     * for the sample will be written
     * @return a list of {@code this.nSamplesPerIndividual()} sampled
     * haplotype pairs for the specified individual
     *
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.gl().nSamples()}
     * @throws IllegalArgumentException if
     * {@code gtProbs.length != this.gl().markers().sumGenotypes()}
     * @throws NullPointerException if {@code gtProbs == null}
     */
    List<HapPair> randomSample(int sample, double[] gtProbs);
}
