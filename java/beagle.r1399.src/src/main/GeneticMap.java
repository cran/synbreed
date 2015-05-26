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
import blbutil.FileIterator;
import blbutil.Filter;
import blbutil.FilterUtils;
import blbutil.InputIterator;
import blbutil.StringUtil;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * <p>Class {@code GeneticMap} represents a genetic map for one or more
 * chromosomes.
 * </p>
 * Instances of class {@code GeneticMap} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class GeneticMap {

    private final int[][] basePos;
    private final double[][] genPos;

    private GeneticMap(List<List<String>> chromList) {
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
            String[] fields=StringUtil.getFields(list.get(j));
            if (fields.length!=4) {
                String s="Map file format error: "+list.get(j);
                throw new IllegalArgumentException(s);
            }
            basePos[j]=Integer.parseInt(fields[3]);
            genPos[j]=Double.parseDouble(fields[2]);
            if (Double.isInfinite(genPos[j]) || Double.isNaN(genPos[j])) {
                String s="invalid map position: "+genPos[j];
                throw new IllegalArgumentException(s);
            }
            if (j>0) {
                if (basePos[j]==basePos[j-1]) {
                    String s="duplication position: "+list.get(j);
                    throw new IllegalArgumentException(s);
                }
                if (basePos[j]<basePos[j-1] || genPos[j]<genPos[j-1]) {
                    String s="map positions not in ascending order: "
                            +list.get(j);
                    throw new IllegalArgumentException(s);
                }
            }
        }
        if (n>0 && genPos[0]==genPos[n-1]) {
            String s="constant map position: " + list.get(0);
            throw new IllegalArgumentException(s);
        }
    }

    /**
     * Constructs and returns a new {@code GeneticMap} instance.
     *
     * @param mapFile a genetic map file in PLINK format.
     * @return a new {@code GeneticMap} instance.
     *
     * @throws IllegalArgumentException if any map position is infinite
     * or {@code NaN}
     * @throws NullPointerException if {@code mapFile==null}
     * @throws NumberFormatException if the base position on any line of the map
     * file is not a parsable integer, or if the genetic map position on any
     * line of the map file is not a parsable double.
     * @throws IllegalArgumentException if a lines of the specified genetic map
     * file does not contain 4 fields, if the map positions on each
     * chromosome are not sorted in ascending order, if there are duplicate
     * base positions, or if all base positions on a chromosome correspond
     * to the same genetic map position.
     */
    public static GeneticMap fromPlinkMapFile(File mapFile) {
        Filter<String> chromFilter = FilterUtils.acceptAllFilter();
        return new GeneticMap(divideByChrom(mapFile, chromFilter));
    }

    /**
     * Constructs and returns a new {@code GeneticMap} instance.
     *
     * @param mapFile a genetic map file in PLINK format.
     * @param chrom a chromosome. The returned genetic map will contain only
     * positions on the specified chromosome.
     * @return a new {@code GeneticMap} instance.
     *
     * @throws IllegalArgumentException if any map position is infinite or
     * {@code NaN}.
     * @throws NullPointerException if {@code mapFile==null || chrom==null}
     * @throws NumberFormatException if the base position on any line of the map
     * file is not a parsable integer, or if the genetic map position on any
     * line of the map file is not a parsable number.
     * @throws IllegalArgumentException if a lines of the specified genetic map
     * file does not contain 4 fields, if the map positions on each
     * chromosome are not sorted in ascending order, if there are duplicate
     * base positions, or if all base positions on the chromosome correspond
     * to the same genetic map position.
     */
    public static GeneticMap fromPlinkMapFile(File mapFile, String chrom) {
        Filter<String> chromFilter = FilterUtils.singletonFilter(chrom);
        return new GeneticMap(divideByChrom(mapFile, chromFilter));
    }

    private static List<List<String>> divideByChrom(File mapFile,
            Filter<String> chromFilter) {
        int initialMapSize = 200;
        List<List<String>> chromList=new ArrayList<>(25);
        try (FileIterator<String> it=InputIterator.fromTextFile(mapFile)) {
            while (it.hasNext()) {
                String line=it.next();
                String[] fields=StringUtil.getFields(line, 2);
                if (fields.length>0) {
                    if (fields.length<2) {
                        String s="Map file format error: "+line;
                        throw new IllegalArgumentException(s);
                    } else {
                        String chrom = fields[0];
                        if (chromFilter.accept(chrom)) {
                            int chromIndex=ChromIds.instance().indexOf(fields[0]);
                            while (chromIndex>=chromList.size()) {
                                chromList.add(new ArrayList<String>(initialMapSize));
                            }
                            chromList.get(chromIndex).add(line);
                        }
                    }
                }
            }
        }
        return chromList;
    }

    /**
     * Returns {@code true} if there is a genetic map for the
     * specified chromosome and returns {@code false} otherwise.
     *
     * @param chrom a chromosome index.
     * @return {@code true} if there is a genetic map for the
     * specified chromosome and {@code false} otherwise.
     */
    public boolean hasMap(int chrom) {
        if (chrom<0 || chrom>=basePos.length) {
            return false;
        }
        else {
            return basePos[chrom].length>0;
        }
    }

    /**
     * Returns the number of genetic map positions.
     *
     * @param chrom a chromosome index.
     * @return the number of genetic map positions.
     */
    public int nMapPositions(int chrom) {
        return hasMap(chrom) ?  basePos[chrom].length : 0;
    }

    /**
     * Returns the specified base position
     *
     * @param chrom a chromosome index.
     * @param index a map position index.
     * @return the specified base position.
     *
     * @throws IllegalArgumentException if {@code this.hasMap(chrom)==false}.
     * @throws IndexOutOfBoundsException if
     * {@code index<0 || index>=this.nMapPositions(chrom)}.
     */
    public int index2BasePos(int chrom, int index) {
        if (hasMap(chrom)==false) {
            throw new IllegalArgumentException("hasMap(chrom)==false");
        }
        return basePos[chrom][index];
    }

    /**
     * Returns the specified genetic position
     *
     * @param chrom a chromosome index.
     * @param index a map position index.
     * @return the specified genetic position
     *
     * @throws IllegalArgumentException if {@code this.hasMap(chrom)==false}
     * @throws IndexOutOfBoundsException if
     * {@code index<0 || index>=this.nMapPositions(chrom)}
     */
    public double index2GenPos(int chrom, int index) {
        if (hasMap(chrom)==false) {
            throw new IllegalArgumentException("hasMap(chrom)==false");
        }
        return genPos[chrom][index];
    }

    /**
     * Returns the map position index that is closes to the specified base
     * position.
     *
     * @param chrom a chromosome index.
     * @param basePosition a base position
     * @return the map position index that is closes to the specified base
     * position.
     *
     * @throws IllegalArgumentException if {@code this.hasMap(chrom)==false}
     */
    public int closestIndex(int chrom, int basePosition) {
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

    /**
     * Returns the genetic position of the specified marker. The map position is
     * estimated using linear interpolation from the nearest genetic map
     * position.
     *
     * @param marker a genetic marker.
     * @return the genetic position of the specified marker.
     * @throws IllegalArgumentException if {@code this.hasMap(chrom)==false}
     * @throws NullPointerException if {@code marker==null}
     */
    public double genPos(Marker marker) {
        return genPos(marker.chromIndex(), marker.pos());
    }

    /**
     * Returns the genetic position of the specified genome coordinate.
     * The map position is estimated using linear interpolation.
     *
     * @param chrom the chromosome index.
     * @param basePosition the base position on the chromosome.
     * @return the genetic map position of the specified base position.
     * @throws IllegalArgumentException if {@code this.hasMap(chrom)==false}
     * @throws NullPointerException if {@code chrom==null}
     */
    public double genPos(int chrom, int basePosition) {
        if (hasMap(chrom)==false) {
            throw new IllegalArgumentException("hasMap(chrom)==false");
        }
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
            return fa + ( ((x-a)/(b-a)) * (fb-fa) );
        }
    }

    /**
     * Returns the base position corresponding to the specified genetic
     * position. If the genetic position is not a map position then the base
     * position is estimated using linear interpolation.
     *
     * @param chrom the chromosome index.
     * @param geneticPosition the genetic position on the chromosome.
     * @return the base position corresponding to the specified genetic
     * position.
     * @throws IllegalArgumentException if {@code this.hasMap(chrom)==false}
     */
    public int basePos(int chrom, double geneticPosition) {
        if (hasMap(chrom)==false) {
            throw new IllegalArgumentException("hasMap(chrom)==false");
        }
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

    /**
     * Returns a string representation of this genetic map. The exact details
     * of the representation are unspecified and subject to change.
     *
     * @return a string representation of this genetic map.
     */
    @Override
    public String toString() {
        return "main.GeneticMap";
    }
}
