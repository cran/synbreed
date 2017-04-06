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
package main;

import beagleutil.Samples;
import blbutil.FileIt;
import blbutil.InputIt;
import blbutil.StringUtil;
import java.io.File;
import java.util.Arrays;

/**
 * <p>Class {@code NuclearFamilies} stores parent-offspring relationships
 * in a list of samples.  In particular, class {@code NuclearFamilies}
 * stores a list of the single individuals in the list of samples,
 * a list of the parent-offspring duos in the list of samples, and a list of
 * the parent-offspring trios in the list of samples. A single individual is
 * an  individuals without a parent or offspring in the list of samples.
 * </p>
 * <p>Instances of class {@code NuclearFamilies} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class NuclearFamilies {

    private final Samples samples;

    private final int[] single;
    private final int[] duoOffspring;
    private final int[] trioOffspring;

    private final int[] mother;
    private final int[] father;

    /**
     * Constructs a new {@code NuclearFamilies} instance.
     *
     * @param samples the list of samples.
     * @param pedFile a linkage-format pedigree file, or {@code null}
     * if no pedigree relationships are known.  A pedigree file must have
     * at least 4 white-space delimited columns.  The first column of the
     * pedigree file (family ID) is ignored.  The second, third, and fourth
     * columns are the individual's ID, the individual's father's ID, and
     * the individual's mother's ID respectively.
     *
     * @throws IllegalArgumentException if a pedigree file is specified,
     * and if the file has a non-blank line with less than 4 white-space
     * delimited fields
     * @throws IllegalArgumentException if a pedigree file is specified,
     * and if the file has duplicate individual identifiers in the
     * second white-space delimited column
     * @throws NullPointerException if {@code samples == null}
     */
    public NuclearFamilies(Samples samples, File pedFile) {
        this.samples = samples;
        this.father = new int[samples.nSamples()];
        this.mother = new int[samples.nSamples()];
        boolean[] isParent = new boolean[samples.nSamples()];
        Arrays.fill(father, -1);
        Arrays.fill(mother, -1);
        if (pedFile != null) {
            identifyParents(samples, pedFile, isParent, father, mother);
        }
        int[] cnts = counts(isParent, father, mother);
        this.single = new int[cnts[0]];
        this.duoOffspring = new int[cnts[1]];
        this.trioOffspring = new int[cnts[2]];
        fillArrays(samples, isParent, father, mother, single, duoOffspring,
                trioOffspring);
    }

    private static void identifyParents(Samples samples, File pedFile,
            boolean[] isParent, int[] father, int[] mother) {
        String MISSING_PARENT = "0";
        boolean[] idHasBeenProcessed = new boolean[samples.nSamples()];
        try (FileIt<String> pedIt=InputIt.fromGzipFile(pedFile)) {
            while (pedIt.hasNext()) {
                String line = pedIt.next().trim();
                if (line.length() > 0) {
                    String[] fields = getPedFields(line);
                    String offspringId = fields[1];
                    String fatherId = fields[2];
                    String motherId = fields[3];
                    int offspring = samples.index(offspringId);
                    if (offspring != -1) {
                        if (idHasBeenProcessed[offspring]) {
                            String s = "duplicate sample in pedigree file: "
                                    + offspringId;
                            throw new IllegalArgumentException(s);
                        }
                        else {
                            idHasBeenProcessed[offspring] = true;
                        }
                        if (fatherId.equals(MISSING_PARENT)==false) {
                            int sampleIndex = samples.index(fatherId);
                            if (sampleIndex != -1) {
                                isParent[sampleIndex] = true;
                                father[offspring] = sampleIndex;
                            }
                        }
                        if (motherId.equals(MISSING_PARENT)==false) {
                            int sampleIndex = samples.index(motherId);
                            if (sampleIndex != -1) {
                                isParent[sampleIndex] = true;
                                mother[offspring] = sampleIndex;
                            }
                        }
                    }
                }
            }
        }
    }

    private static String[] getPedFields(String line) {
        String[] fields = StringUtil.getFields(line, 5);
        if (fields.length < 4) {
            String s = "invalid line in ped file: " + line;
            throw new IllegalArgumentException(s);
        }
        return fields;
    }

    private int[] counts(boolean[] isParent, int[] fathers, int[] mothers) {
        assert isParent.length==fathers.length;
        assert isParent.length==mothers.length;
        int[] cnts = new int[3];
        for (int j=0; j<isParent.length; ++j) {
            int nParents = 0;
            if (fathers[j] >= 0) {
                ++nParents;
            }
            if (mothers[j] >= 0) {
                ++nParents;
            }
            if (nParents==0) {
                if (isParent[j]==false) {
                    ++cnts[0];  // increment single count, cnts[0]
                }
            }
            else {
                // increment duo count, cnts[1], or trio count, cnt[2]
                ++cnts[nParents];
            }
        }
        return cnts;
    }

    private static void fillArrays(Samples samples, boolean[] isParent,
            int[] father, int[] mother, int[] single,
            int[] duoOffspring, int[] trioOffspring) {
        int singleIndex = 0;
        int duoIndex = 0;
        int trioIndex = 0;
        for (int j=0, n=samples.nSamples(); j<n; ++j) {
            int nParents = nParents(j, father, mother);
            switch (nParents) {
                case 0:
                    if (isParent[j]==false) {
                        single[singleIndex++] = j;
                    }
                    break;
                case 1:
                    duoOffspring[duoIndex++] = j;
                    break;
                case 2:
                    trioOffspring[trioIndex++] = j;
                    break;
                default:
                    assert false;
            }
        }
        assert singleIndex==single.length;
        assert duoIndex==duoOffspring.length;
        assert trioIndex==trioOffspring.length;
    }

    private static int nParents(int index, int[] father, int[] mother) {
        int cnt = 0;
        if (father[index]>=0) {
            ++cnt;
        }
        if (mother[index]>=0) {
            ++cnt;
        }
        return cnt;
    }

    /**
     * Returns the list of samples.
     * @return the list of samples
     */
    public Samples samples() {
        return samples;
    }

    /**
     * Returns the number of samples.
     * @return the number of samples
     */
    public int nSamples() {
        return samples.nSamples();
    }

    /**
     * Returns the number of single individuals in the list of samples.
     * A single individual has no parent or offspring in the list of samples.
     * @return the number of single individuals in the sample
     */
    public int nSingles() {
        return single.length;
    }

    /**
     * Returns the number of parent-offspring duos in the list of samples.
     * The offspring of a parent-offspring duo has only one parent
     * in the sample.
     * @return the number of parent-offspring duos in the list of samples
     */
    public int nDuos() {
        return duoOffspring.length;
    }

    /**
     * Returns the number of parent-offspring trios in the list of samples.
     * The offspring of a parent-offspring trio has two parents
     * in the sample.
     * @return the number of parent-offspring trios in the list of samples
     */
    public int nTrios() {
        return trioOffspring.length;
    }

    /**
     * Returns the sample index of the specified single individual.
     * A single individual has no first-degree relative in the list of
     * samples.
     * @param index the index of a single individual
     * @return the sample index of the specified single individual
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nSingles()}
     */
    public int single(int index) {
        return single[index];
    }

    /**
     * Returns the sample index of the parent of the specified
     * parent-offspring duo.
     * @param index the index of a parent-offspring duo
     * @return the sample index of the parent of the specified
     * parent-offspring duo
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nDuos()}
     */
    public int duoParent(int index) {
        int offspring = duoOffspring[index];
        if (father[offspring]>=0) {
            return father[offspring];
        }
        else {
            assert mother[offspring]>=0;
            return mother[offspring];
        }
    }

    /**
     * Returns the sample index of the offspring of the specified
     * parent-offspring duo.
     * @param index the index of a parent-offspring duo
     * @return the sample index of the offspring of the specified
     * parent-offspring duo
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nDuos()}
     */
    public int duoOffspring(int index) {
        return duoOffspring[index];
    }

    /**
     * Returns the sample index of the father of the specified
     * parent-offspring trio.
     * @param index the index of a parent-offspring trio
     * @return the sample index of the father of the specified
     * parent-offspring trio
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nTrios()}
     */
    public int trioFather(int index) {
        return father[trioOffspring[index]];
    }

    /**
     * Returns the sample index of the mother of the specified
     * parent-offspring trio.
     * @param index the index of a parent-offspring trio
     * @return the sample index of the mother of the specified
     * parent-offspring trio
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nTrios()}
     */
    public int trioMother(int index) {
        return mother[trioOffspring[index]];
    }

    /**
     * Returns the sample index of the offspring of the specified
     * parent-offspring trio.
     * @param index the index of a parent-offspring trio
     * @return the sample index of the offspring of the specified
     * parent-offspring trio
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nTrios()}
     */
    public int trioOffspring(int index) {
        return trioOffspring[index];
    }

    /**
     * Returns the sample index of the father of the specified sample,
     * or returns {@code -1} if the father is unknown or is not present
     * in the list of samples.
     * @param sample a sample index
     * @return the sample index of the father of the specified sample,
     * or {@code -1} if the father is unknown or is not present in
     * the list of samples
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()()}
     */
    public int father(int sample) {
        return father[sample];
    }

    /**
     * Returns the sample index of the mother of the specified sample,
     * or returns {@code -1} if the mother is unknown or is not present
     * in the list of samples.
     * @param sample a sample index
     * @return the sample index of the mother of the specified sample,
     * or {@code -1} if the mother is unknown or is not present
     * in the list of samples
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.nSamples()()}
     */
    public int mother(int sample) {
        return mother[sample];
    }

    /**
     * Returns a string representation of {@code this}.  The exact details of
     * the representation are unspecified and subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        return this.getClass().toString();
    }
}
