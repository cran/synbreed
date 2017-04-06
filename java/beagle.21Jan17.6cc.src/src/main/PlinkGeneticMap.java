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

import vcf.Marker;
import beagleutil.ChromIds;
import blbutil.FileIt;
import blbutil.Filter;
import blbutil.InputIt;
import blbutil.StringUtil;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * <p>Class {@code PlinkGeneticMap} represents a genetic map derived
 * from a PLINK map file with map positions in cM units for one or more
 * chromosomes.
 * </p>
 * <p>Instances of class {@code PlinkGeneticMap} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class PlinkGeneticMap implements GeneticMap {

    private final int[][] basePos;
    private final double[][] genPos;

    private PlinkGeneticMap(List<List<String>> chromList) {
        this.basePos = new int[chromList.size()][];
        this.genPos = new double[chromList.size()][];
        for (int j=0, n=chromList.size(); j<n; ++j) {
            List<String> list = chromList.get(j);
            basePos[j] = new int[list.size()];
            genPos[j] = new double[list.size()];
            fillMapPositions(list, basePos[j], genPos[j]);
        }
    }

    private static void fillMapPositions(List<String> list, int[] basePos,
            double[] genPos) {
        int n = list.size();
        assert basePos.length==n && genPos.length==n;
        for (int j=0; j<n; ++j) {
            String[] fields = StringUtil.getFields(list.get(j));
            if (fields.length!=4) {
                String s = "Map file format error: " + list.get(j);
                throw new IllegalArgumentException(s);
            }
            basePos[j] = Integer.parseInt(fields[3]);
            genPos[j] = Double.parseDouble(fields[2]);
            if (Double.isFinite(genPos[j])==false) {
                String s = "invalid map position: " + genPos[j];
                throw new IllegalArgumentException(s);
            }
            if (j>0) {
                if (basePos[j]==basePos[j-1]) {
                    String s = "duplication position: " + list.get(j);
                    throw new IllegalArgumentException(s);
                }
                if (basePos[j]<basePos[j-1] || genPos[j]<genPos[j-1]) {
                    String s = "map positions not in ascending order: "
                            + list.get(j);
                    throw new IllegalArgumentException(s);
                }
            }
        }
        if (n>0 && genPos[0]==genPos[n-1]) {
            String s = "Genetic map has only one map position: " + list.get(0);
            throw new IllegalArgumentException(s);
        }
    }

    /**
     * Constructs and returns a new {@code PlinkGeneticMap} instance from
     * the data in the specified file.
     *
     * @param mapFile a genetic map file in PLINK format with genetic map
     * positions in cM units
     * @return a new {@code PlinkGeneticMap} instance
     *
     * @throws IllegalArgumentException if any map position is infinite
     * or {@code NaN}
     * @throws NullPointerException if {@code mapFile == null}
     * @throws NumberFormatException if the base position on any line of the map
     * file is not a parsable integer
     * @throws NumberFormatException if the genetic map position on any
     * line of the map file is not a parsable double
     * @throws IllegalArgumentException if a non-empty line of the specified
     * genetic map file does not contain 4 fields
     * @throws IllegalArgumentException if the map positions on each
     * chromosome are not sorted in ascending order
     * @throws IllegalArgumentException if there are duplicate
     * base positions on a chromosome
     * @throws IllegalArgumentException if all base positions on a chromosome
     * have the same genetic map position
     */
    public static PlinkGeneticMap fromPlinkMapFile(File mapFile) {
        Filter<String> chromFilter = Filter.acceptAllFilter();
        return new PlinkGeneticMap(divideByChrom(mapFile, chromFilter));
    }

    /**
     * Constructs and returns a new {@code PlinkGeneticMap} instance from
     * the data in the specified file. The returned genetic map will contain
     * only positions on the specified chromosome
     *
     * @param mapFile a genetic map file in PLINK format with genetic map
     * positions in cM units
     * @param chrom a chromosome
     * @return a new {@code PlinkGeneticMap} instance
     *
     * @throws IllegalArgumentException if any map position is infinite or
     * {@code NaN}.
     * @throws NullPointerException if {@code mapFile == null || chrom == null}
     * @throws NumberFormatException if the base position on a line of the map
     * file that corresponds to the specified chromosome is not a parsable
     * integer
     * @throws NumberFormatException if the genetic map position on a line
     * of the map file that corresponds to the specified chromosome is not
     * a parsable double
     * @throws IllegalArgumentException if a non-empty line of the specified
     * genetic map file does not contain 4 fields
     * @throws IllegalArgumentException if the map positions on the specified
     * chromosome are not sorted in ascending order
     * @throws IllegalArgumentException if there are duplicate base positions
     * on the specified chromosome
     * @throws IllegalArgumentException if all base positions on the
     * specified chromosome have the same genetic map position
     * @throws IllegalArgumentException if the specified chromosome does not
     * have at least two distinct positions in the genetic map
     */
    public static PlinkGeneticMap fromPlinkMapFile(File mapFile, String chrom) {
        chrom = chrom.trim();
        Filter<String> chromFilter = singletonFilter(chrom);
        return new PlinkGeneticMap(divideByChrom(mapFile, chromFilter));
    }

    /**
     * Returns a filter that accepts only objects that are equal
     * to the specified object.
     * @param <E> the type of object that is filtered
     * @param singleton the object that will be accepted
     * @return a filter that accepts only objects that are equal
     * to the specified object
     * @throws NullPointerException if {@code singleton == null}
     */
   private static <E> Filter<E> singletonFilter(final E singleton) {
        if (singleton==null) {
            throw new NullPointerException("singleton==null");
        }
        return (E e) -> {
            if (e==null) {
                throw new NullPointerException("e==null");
            }
            return singleton.equals(e);
        };
    }

    private static List<List<String>> divideByChrom(File mapFile,
            Filter<String> chromFilter) {
        int initialMapSize = 200;
        List<List<String>> chromList = new ArrayList<>(25);
        try (FileIt<String> it = InputIt.fromTextFile(mapFile)) {
            while (it.hasNext()) {
                String line = it.next();
                String[] fields = StringUtil.getFields(line, 4);
                if (fields.length > 0) {
                    if (fields.length < 4) {
                        String s = "Map file format error: " + line;
                        throw new IllegalArgumentException(s);
                    } else {
                        String chrom = fields[0];
                        if (chromFilter.accept(chrom)) {
                            int chromIndex = ChromIds.instance().getIndex(fields[0]);
                            while (chromIndex >= chromList.size()) {
                                chromList.add(new ArrayList<>(initialMapSize));
                            }
                            chromList.get(chromIndex).add(line);
                        }
                    }
                }
            }
        }
        return chromList;
    }

    private void checkChromIndex(int chrom) {
        if (chrom < 0 || chrom >= ChromIds.instance().size()) {
            throw new IndexOutOfBoundsException(String.valueOf(chrom));
        }
        if (chrom>=basePos.length || basePos[chrom].length == 0) {
            String s = "missing genetic map for chromosome "
                    + ChromIds.instance().id(chrom);
            throw new IllegalArgumentException(s);
        }
    }

    /**
     * Returns the number of mapped loci in this genetic map.
     *
     * @param chrom a chromosome index
     * @return the number of mapped loci in this genetic map
     * @throws IllegalArgumentException if this genetic map has no
     * map positions for the specified chromosome
     * @throws IndexOutOfBoundsException if
     * {@code chrom < 0 || chrom >= ChromIds.instance().size()}
     */
    public int nMapPositions(int chrom) {
        checkChromIndex(chrom);
        return basePos[chrom].length;
    }

    /**
     * Returns the specified base position
     *
     * @param chrom a chromosome index
     * @param index a map position index
     * @return the specified base position
     *
     * @throws IllegalArgumentException if this genetic map has no
     * map positions for the specified chromosome
     * @throws IndexOutOfBoundsException if
     * {@code chrom < 0 || chrom >= ChromIds.instance().size()}
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nMapPositions(chrom)}
     */
    public int index2BasePos(int chrom, int index) {
        checkChromIndex(chrom);
        return basePos[chrom][index];
    }

    /**
     * Returns the specified genetic map position
     *
     * @param chrom a chromosome index
     * @param index a map position index
     * @return the specified genetic map position
     *
     * @throws IllegalArgumentException if this genetic map has no
     * map positions for the specified chromosome
     * @throws IndexOutOfBoundsException if
     * {@code chrom < 0 || chrom >= ChromIds.instance().size()}
     * @throws IndexOutOfBoundsException if
     * {@code index < 0 || index >= this.nMapPositions(chrom)}
     */
    public double index2GenPos(int chrom, int index) {
        checkChromIndex(chrom);
        return genPos[chrom][index];
    }

    /**
     * Returns the index of the genetic map position that is closest to the
     * specified base position.
     *
     * @param chrom a chromosome index
     * @param basePosition a base position
     * @return the genetic map position index that is closes to the
     * specified base position.
     *
     * @throws IllegalArgumentException if this genetic map has no
     * map positions for the specified chromosome
     * @throws IndexOutOfBoundsException if
     * {@code chrom < 0 || chrom >= ChromIds.instance().size()}
     */
    public int closestIndex(int chrom, int basePosition) {
        checkChromIndex(chrom);
        assert basePos.length>=2;
        int mapIndex = Arrays.binarySearch(basePos[chrom], basePosition);
        if (mapIndex >= 0) {
            return mapIndex;
        } else {
            int insPt = -mapIndex-1;
            if (insPt==0) {
                return 0;
            } else if (insPt==basePos.length) {
                return basePos.length-1;
            } else {
                int distInsPt = basePos[chrom][insPt] - basePosition;
                int distInsPtM1 = basePosition - basePos[chrom][insPt-1];
                return (distInsPt<=distInsPtM1) ? insPt : (insPt-1);
            }
        }
    }

    @Override
    public double genPos(Marker marker) {
        return genPos(marker.chromIndex(), marker.pos());
    }

    @Override
    public double genPos(int chrom, int basePosition) {
        checkChromIndex(chrom);
        assert basePos[chrom].length>=2;
        assert basePos[chrom].length==genPos[chrom].length;
        int index = Arrays.binarySearch(basePos[chrom], basePosition);
        if (index>=0) {
            return genPos[chrom][index];
        } else {
            int insPt = -index-1;
            if (insPt==basePos[chrom].length) {
                --insPt;
            }
            else if (insPt==0) {
                ++insPt;
            }
            int x = basePosition;
            int a = basePos[chrom][insPt-1];
            int b = basePos[chrom][insPt];
            double fa = genPos[chrom][insPt-1];
            double fb = genPos[chrom][insPt];
            return fa + ( ( (double) (x-a)/ (double) (b-a)) * (fb-fa) );
        }
    }

    @Override
    public int basePos(int chrom, double geneticPosition) {
        checkChromIndex(chrom);
        assert basePos[chrom].length>=2;
        assert basePos[chrom].length==genPos[chrom].length;
        int index = Arrays.binarySearch(genPos[chrom], geneticPosition);
        if (index>=0) {
            return basePos[chrom][index];
        } else {
            int insPt = -index-1;
            if (insPt==genPos[chrom].length) {
                --insPt;
                while (genPos[chrom][insPt]==genPos[chrom][insPt-1]) {
                    --insPt;
                }
            }
            else if (insPt==0) {
                ++insPt;
                while (genPos[chrom][insPt]==genPos[chrom][insPt-1]) {
                    ++insPt;
                }
            }
            double x = geneticPosition;
            double a = genPos[chrom][insPt-1];
            double b = genPos[chrom][insPt];
            int fa = basePos[chrom][insPt-1];
            int fb = basePos[chrom][insPt];
            return fa + (int) Math.round( ((x-a)/(b-a)) * (fb-fa) );
        }
    }

    @Override
    public String toString() {
        return this.getClass().toString();
    }
}
