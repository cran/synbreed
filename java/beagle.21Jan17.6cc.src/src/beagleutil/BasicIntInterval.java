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
package beagleutil;

/**
 * <p>Class {@code BasicIntInterval} represents an interval of
 * consecutive integers.
 * </p>
 * Instances of class {@code BasicIntInterval} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}

} */
public final class BasicIntInterval implements IntInterval,
        Comparable<BasicIntInterval> {

    private final int start;
    private final int end;

    /**
     * Constructs an {@code SimpleIntInterval} instance.
     * @param start the starting integer (inclusive).
     * @param end the ending integer (inclusive).
     * @throws IllegalArgumentException if {@code start>end}.
     */
    public BasicIntInterval(int start, int end) {
        if (start > end) {
            String s = "start=" + start + " end=" + end;
            throw new IllegalArgumentException(s);
        }
        this.start = start;
        this.end = end;
    }

    @Override
    public int start() {
        return start;
    }

    @Override
    public int end() {
        return end;
    }

    /**
     * <p>Returns a hash code value for the object.
     * </p>
     * <p>The hash code is defined by the following calculation:
     * </p>
     * <pre>
        int hash = 3;
        hash += 59 * hash + this.start();
        hash += 59 * hash + this.end();
     * </pre>
     * @return a hash code value for the object.
     */
    @Override
    public int hashCode() {
        int hash = 3;
        hash = 59 * hash + this.start;
        hash = 59 * hash + this.end;
        return hash;
    }

    /**
     * Returns {@code true} if the specified object is an
     * {@code BasicIntInterval} instance and
     * {@code this.start() == ((BasicIntInterval) obj).start()}, and
     * {@code this.end() == ((BasicIntInterval) obj).end()},
     * and returns {@code false} otherwise.
     *
     * @param obj the object to be compared with {@code this} for equality.
     * @return {@code true} if the specified object is equal to
     * {@code this}, and returns false otherwise.
     */
    @Override
    public boolean equals(Object obj) {
        if (this==obj) {
            return true;
        }
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final BasicIntInterval other = (BasicIntInterval) obj;
        if (this.start != other.start) {
            return false;
        }
        return this.end==other.end;
    }

    /**
     * Compares the specified {@code BasicIntInterval} with this for order,
     * and returns a negative integer, zero, or a positive integer as
     * {@code this} is less than, equal to, or greater than the specified
     * {@code BasicIntInterval} object.
     * {@code BasicIntInterval} objects are
     * ordered by their start and their end values in that order.
     *
     * @param o the {@code BasicIntInterval} to be compared with this.
     *
     * @return a negative integer, zero, or a positive integer as this
     * object is less than, equal to, or greater than the specified object.
     */
    @Override
    public int compareTo(BasicIntInterval o) {
        if (this.start != o.start) {
            return this.start < o.start ? -1 : 1;
        }
        else if (this.end != o.end) {
            return this.end < o.end ? -1 : 1;
        }
        return 0;
    }

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     *
     * @return a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append('[');
        sb.append(start);
        sb.append(", ");
        sb.append(end);
        sb.append(']');
        return sb.toString();
    }
}
