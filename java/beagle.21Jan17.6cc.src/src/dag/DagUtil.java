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
package dag;

import blbutil.Const;
import java.text.DecimalFormat;

/**
 * <p>Class {@code DagUtil} contains static, thread-safe methods for
 * removing elements of an array that have a specified value and
 * for creating a string representation of a DAG.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class DagUtil {

    private DagUtil() {
        // private constructor to prevent instantiation.
    }

    /**
     * Return a string description of the specified {@code DAG} object.  The
     * exact details of the description are unspecified and subject
     * to change.
     *
     * @param dag a directed acyclic graph
     * @return a string description of the specified {@code DAG}
     * @throws NullPointerException if {@code dag == null}
     */
    public static String dagStats(Dag dag) {
        DecimalFormat df2 = new DecimalFormat("0.000");
        int fieldWidth = 7;
        double nEdges = dag.nEdges();
        double nNodes = dag.nNodes();
        double nMarkers = dag.nLevels();

        long nEdgesPerLevel = (int) Math.round(nEdges/nMarkers);
        double nEdgesPerNode = nEdges / nNodes;
        int nHapsPerEdge = (int) Math.round(dag.parentWeight(0, 0) / nEdgesPerLevel);

        String edgesPerLevel = String.valueOf(nEdgesPerLevel);
        String maxEdgesPerLevel = String.valueOf(dag.maxEdges());
        String edgesPerNode = df2.format(nEdgesPerNode);
        String hapsPerNode = String.valueOf(nHapsPerEdge);

        StringBuilder sb = new StringBuilder(100);
        sb.append("mean edges/level: ");
        sb.append(edgesPerLevel);
        padField(sb, fieldWidth - edgesPerLevel.length());

        sb.append("max edges/level: ");
        sb.append(maxEdgesPerLevel);
        sb.append(Const.nl);

        sb.append("mean edges/node:  ");
        sb.append(edgesPerNode);
        padField(sb, fieldWidth - edgesPerNode.length());

        sb.append("mean count/edge: ");
        sb.append(hapsPerNode);
        sb.append(Const.nl);

        return sb.toString();
    }

    private static void padField(StringBuilder sb, int nSpaces) {
        for (int j=0; j<nSpaces; ++j) {
            sb.append(" ");
        }
    }

    /**
     * Returns the number of the elements in the specified array that
     * equal the specified value.
     * @param array an array of integers
     * @param value an integer value
     * @return the number of the elements in the specified array that
     * equal the specified value
     *
     * @throws NullPointerException if {@code array == null}
     */
    private static int count(int[] array, int value) {
        int cnt = 0;
        for (int i : array) {
            if (i == value) {
                ++cnt;
            }
        }
        return cnt;
    }

    /**
     * Returns the number of the elements in the specified array that
     * equal the specified value.
     * @param array an array of float values
     * @param value a float value
     * @return the number of the elements in the specified array that
     * equal the specified value
     *
     * @throws IllegalArgumentException if {@code Float.isNaN(value) == true}
     * @throws NullPointerException if {@code array == null}
     */
    private static int count(float[] array, float value) {
        if (Float.isNaN(value)) {
            throw new IllegalArgumentException(String.valueOf(value));
        }
        int cnt = 0;
        for (float f : array) {
            if (f == value) {
                ++cnt;
            }
        }
        return cnt;
    }

    /**
     * Returns an array obtained by removing all elements in the
     * specified array that equal the specified value.
     * @param array an array of integers
     * @param value an integer value
     * @return an array obtained by removing all elements in the
     * specified array that equal the specified value
     *
     * @throws NullPointerException if {@code array == null}
     */
    public static int[] removeValues(int[] array, int value) {
        int cnt = DagUtil.count(array, value);
        int[] reducedArray = new int[array.length - cnt];
        int index=0;
        for (int j=0; j<array.length; ++j) {
            if (array[j] != value) {
                reducedArray[index++] = array[j];
            }
        }
        assert index==reducedArray.length;
        return reducedArray;
    }

    /**
     * Returns an array obtained by removing all elements in the
     * specified array that equal the specified value.
     * @param array an array of float values
     * @param value a float value
     * @return an array obtained by removing all elements in the
     * specified array that equal the specified value
     *
     * @throws IllegalArgumentException if {@code Float.isNaN(value) == true}
     * @throws NullPointerException if {@code array == null}
     */
    public static float[] removeValues(float[] array, float value) {
        if (Float.isNaN(value)) {
            throw new IllegalArgumentException(String.valueOf(value));
        }
        int cnt = DagUtil.count(array, value);
        float[] reducedArray = new float[array.length - cnt];
        int index=0;
        for (int j=0; j<array.length; ++j) {
            if (array[j] != value) {
                reducedArray[index++] = array[j];
            }
        }
        assert index==reducedArray.length;
        return reducedArray;
    }
}
