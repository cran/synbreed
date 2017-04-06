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

import blbutil.Utilities;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import haplotype.HapPair;
import haplotype.RevHapPair;
import main.GenotypeValues;

/**
 * <p>Class {@code ConsumeSingleSamples} samples haplotype pairs conditional
 * on the observed genotype data and a haplotype frequency model.
 * Class {@code ConsumeSingleSamples} is designed for use as a consumer in a
 * producer-consumer design pattern.
 * </p>
 * <p>Instances of class {@code ConsumeSingleSamples} are thread-safe if the
 * synchronization requirements for the constructor are satisfied.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class ConsumeSingleSamples implements Runnable {

    /**
     * A sentinel {@code Integer}.
     */
    public static final Integer POISON = -1;

    private final boolean markersAreReversed;
    private final SingleBaumInterface baum;
    private final BlockingQueue<Integer> qIn;
    private final List<HapPair> sampledHaps;
    private final GenotypeValues gv;
    private final double[] gprobs;

    /**
     * Constructs a new {@code ConsumeSingleSample} instance.
     *
     * @param markersAreReversed {@code true} if the {@code baum} parameter
     * {@code randomSample()} method produces sampled haplotype pairs that have
     * their marker order reversed and {@code false} otherwise
     * @param baum a thread-confined instance of class
     * {@code sample.SingleBaumInterface}
     * @param qIn a thread-safe input work queue
     * @param hapList a thread-safe list for storing sampled haplotype pairs
     *
     * @throws NullPointerException if any parameter is {@code null}
     */
    public ConsumeSingleSamples(boolean markersAreReversed,
            SingleBaumInterface baum, BlockingQueue<Integer> qIn,
            List<HapPair> hapList) {
        if (baum == null) {
            throw new NullPointerException("baum=null");
        }
        if (qIn == null) {
            throw new IllegalArgumentException("qIn==null");
        }
        if (hapList == null) {
            throw new IllegalArgumentException("hapList==null");
        }
        this.markersAreReversed = markersAreReversed;
        this.baum = baum;
        this.qIn = qIn;
        this.sampledHaps = hapList;
        this.gv = null;
        this.gprobs = null;
    }

    /**
     * Constructs a new {@code ConsumeSingleSample} instance.
     *
     * @param markersAreReversed {@code true} if the {@code baum} parameter
     * {@code randomSample()} method produces sampled haplotype pairs that have
     * their marker order reversed and {@code false} otherwise
     * @param baum a thread-confined instance of class
     * {@code sample.SingleBaumInterface}
     * @param qIn a thread-safe input work queue
     * @param hapList a thread-safe list for storing sampled haplotype pairs
     * @param gv a thread-safe object which stores scaled posterior genotype
     * probabilities
     *
     * @throws NullPointerException if any parameter is {@code null}
     */
    public ConsumeSingleSamples(boolean markersAreReversed,
            SingleBaumInterface baum, BlockingQueue<Integer> qIn,
            List<HapPair> hapList, GenotypeValues gv) {
        if (baum == null) {
            throw new NullPointerException("baum=null");
        }
        if (qIn == null) {
            throw new IllegalArgumentException("qIn==null");
        }
        if (hapList == null) {
            throw new IllegalArgumentException("hapList==null");
        }
        if (gv == null) {
            throw new IllegalArgumentException("gv==null");
        }
        this.markersAreReversed = markersAreReversed;
        this.baum = baum;
        this.qIn = qIn;
        this.gv = gv;
        this.sampledHaps = hapList;
        int n = baum.gl().markers().sumGenotypes();
        this.gprobs = new double[n];
    }

    /**
     * Takes sample indices from the thread-safe work-queue specified at time of
     * construction and samples haplotype pairs for each sample. The method
     * exits when {@code ConsumeSingleSamples.POISON} is taken from
     * the work queue.
     */
    @Override
    @SuppressWarnings({"BroadCatchBlock", "TooBroadCatch"})
    public void run() {
        try {
            int single = qIn.take();
            while (single != POISON) {
                if (gv == null) {
                    List<HapPair> newHaps = baum.randomSample(single);
                    storeHaps(newHaps);
                } else {
                    List<HapPair> newHaps = baum.randomSample(single, gprobs);
                    storeHaps(newHaps);
                    gv.add(single, gprobs);

                }
                single = qIn.take();
            }
        } catch (Throwable e) {
            Utilities.exit("ConsumeSingleSamples: ERROR", e);
        }
    }

    private void storeHaps(List<HapPair> newHaps) {
        if (markersAreReversed) {
            newHaps.stream().forEach((hp) -> {
                sampledHaps.add(new RevHapPair(hp));
            });
        } else {
            sampledHaps.addAll(newHaps);
        }
    }
}
