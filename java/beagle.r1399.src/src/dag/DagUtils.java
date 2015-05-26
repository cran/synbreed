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
import java.util.concurrent.atomic.AtomicIntegerArray;

/**
 * Class {@code DagUtils} contains static methods for counting and
 * removing elements of an array that have a specified value.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class DagUtils {

    private DagUtils() {
        // private constructor to prevent instantiation.
    }

    /**
     * Return a string description of the specified {@code DAG} object.  The
     * exact details of the description are unspecified and subject
     * to change.
     *
     * @param dag a directed acyclic graph.
     * @return a string description of the specified {@code DAG}.  The
     * exact details of the description are unspecified and subject
     * to change.
     * @throws NullPointerException if {@code dag==null}
     */
    public static String dagStats(Dag dag) {
        DecimalFormat df2 = new DecimalFormat("0.000");
        int fieldWidth = 7;
        double nEdges = dag.nEdges();
        double nNodes = dag.nNodes();
        double nMarkers = dag.nMarkers();

        long nEdgesPerLevel = (int) Math.round(nEdges/nMarkers);
        double nEdgesPerNode = nEdges / nNodes;
        int nHapsPerEdge = (int) Math.round(dag.nodeCnt(0, 0) / nEdgesPerLevel);

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
     * @param array an array of bytes.
     * @param value a byte value.
     * @return the number of the elements in the specified array that
     * equal the specified value.
     *
     * @throws NullPointerException if {@code array==null}.
     */
    public static int count(byte[] array, byte value) {
        int cnt = 0;
        for (byte i : array) {
            if (i == value) {
                ++cnt;
            }
        }
        return cnt;
    }

    /**
     * Returns the number of the elements in the specified array that
     * equal the specified value.
     * @param array an array of integers.
     * @param value an integer value.
     * @return the number of the elements in the specified array that
     * equal the specified value.
     *
     * @throws NullPointerException if {@code array==null}.
     */
    public static int count(int[] array, int value) {
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
     * @param array an array of floats.
     * @param value a float value.
     * @return the number of the elements in the specified array that
     * equal the specified value.
     *
     * @throws IllegalArgumentException if {@code Float.isNaN(value)==true}
     * @throws NullPointerException if {@code array==null}
     */
    public static int count(float[] array, float value) {
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
     * @param array an array of bytes.
     * @param value a byte value.
     * @return an array obtained by removing all elements in the
     * specified array that equal the specified value.
     *
     * @throws NullPointerException if {@code array==null}
     */
    public static byte[] removeValues(byte[] array, byte value) {
        int cnt = DagUtils.count(array, value);
        byte[] reducedArray = new byte[array.length - cnt];
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
     * @param array an array of integers.
     * @param value an integer value.
     * @return an array obtained by removing all elements in the
     * specified array that equal the specified value.
     *
     * @throws NullPointerException if {@code array==null}
     */
    public static int[] removeValues(int[] array, int value) {
        int cnt = DagUtils.count(array, value);
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
     * @param array an array of floats.
     * @param value a float value.
     * @return an array obtained by removing all elements in the
     * specified array that equal the specified value.
     *
     * @throws IllegalArgumentException if {@code Float.isNaN(value)==true}
     * @throws NullPointerException if {@code array==null}
     */
    public static float[] removeValues(float[] array, float value) {
        if (Float.isNaN(value)) {
            throw new IllegalArgumentException(String.valueOf(value));
        }
        int cnt = DagUtils.count(array, value);
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

    /**
     * Returns the number of the elements in the specified array that
     * equal the specified value.
     * @param array an array of integers.
     * @param value an integer value.
     * @return the number of the elements in the specified array that
     * equal the specified value.
     *
     * @throws NullPointerException if {@code array==null}
     */
    public static int count(AtomicIntegerArray array, int value) {
        int cnt = 0;
        for (int j=0, n=array.length(); j<n; ++j) {
            if (array.get(j) == value) {
                ++cnt;
            }
        }
        return cnt;
    }

    /**
     * Returns an array obtained by removing all elements in the
     * specified array that equal the specified value.
     * @param array an array of integers.
     * @param value an integer value.
     * @return an array obtained by removing all elements in the
     * specified array that equal the specified value.
     *
     * @throws NullPointerException if {@code array==null}
     */
    public static int[] removeValues(AtomicIntegerArray array, int value) {
        int cnt = DagUtils.count(array, value);
        int[] reducedArray = new int[array.length() - cnt];
        int index=0;
        for (int j=0, n=array.length(); j<n; ++j) {
            int element = array.get(j);
            if (element != value) {
                reducedArray[index++] = element;
            }
        }
        assert index==reducedArray.length;
        return reducedArray;
    }

    /**
     * Returns an array obtained by removing all elements in the
     * specified array that equal the specified value.
     * @param array an array of bytes.
     * @param value a byte value.
     * @return an array obtained by removing all elements in the
     * specified array that equal the specified value.
     *
     * @throws NullPointerException if {@code array==null}
     */
    public static byte[] removeValues(AtomicIntegerArray array, byte value) {
        int cnt = DagUtils.count(array, value);
        byte[] reducedArray = new byte[array.length() - cnt];
        int index=0;
        for (int j=0, n=array.length(); j<n; ++j) {
            byte element = (byte) array.get(j);
            if (element != value) {
                reducedArray[index++] = element;
            }
        }
        assert index==reducedArray.length;
        return reducedArray;
    }


    /**
     * Returns an array obtained by removing all elements in the
     * specified array that equal the specified value.
     * @param array an array of floats.  Each float is encoded using the
     * the {@code Float.intBitsToFloat()} method.
     * @param value an float value.
     * @return an array obtained by removing all elements in the
     * specified array that equal the specified value.
     *
     * @throws IllegalArgumentException if {@code Float.isNaN(value)==true}
     * @throws NullPointerException if {@code array==null}
     */
    public static float[] removeValues(AtomicIntegerArray array, float value) {
        if (Float.isNaN(value)) {
            throw new IllegalArgumentException(String.valueOf(value));
        }
        int intBits = Float.floatToIntBits(value);
        int cnt = DagUtils.count(array, intBits);
        float[] reducedArray = new float[array.length() - cnt];
        int index=0;
        for (int j=0, n=array.length(); j<n; ++j) {
            int element = array.get(j);
            if (element != intBits) {
                reducedArray[index++] = Float.intBitsToFloat(element);
            }
        }
        assert index==reducedArray.length;
        return reducedArray;
    }


}
