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
package blbutil;

/**
 * Class {@code Pair} represents a pair of ordered objects.
 *
 * @param <F> the type of the first object in the pair.
 * @param <S> the type of the second object in the pair.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
*/
public final class Pair<F,S> {

    private final F first;
    private final S second;

    /**
     * Constructs a {@code Pair} instance representing the specified
     * ordered pair of objects.
     * @param first the first object in the ordered pair of objects.
     * @param second the second object in the ordered pair of objects.
     * @throws NullPointerException if {@code first==null || second==null}.
     */
    public Pair(F first, S second) {
        if (first==null) {
            throw new NullPointerException("first==null");
        }
        if (second==null) {
            throw new NullPointerException("second==null");
        }
        this.first = first;
        this.second = second;
    }

    /**
     * Returns the first object in the ordered pair of objects.
     * @return the first object in the ordered pair of objects.
     */
    public F first() {
        return first;
    }

    /**
     * Returns the second object in the ordered pair of objects.
     * @return the second object in the ordered pair of objects.
     */
    public S second() {
        return second;
    }

    /**
     * Returns a hash code value for the object.
     * @return a hash code value for the object.
     */
    @Override
    public int hashCode() {
        int hash = 17;
        hash = 37*hash + first.hashCode();
        hash = 37*hash + second.hashCode();
        return hash;
    }

    /**
     * Returns {@code true} if the specified object is a
     * {@code Pair} instance representing the same
     * ordered pair of objects as {@code this}, and
     * returns {@code false} otherwise.  Two ordered pairs,
     * {@code p1} and {@code p2}, are equal if
     * {@code p1.first().equals(p2.first())
     * && p1.second().equals(p2.second())}.
     *
     * @param obj the object to be compared with {@code this} for
     * equality.
     * @return  {@code true} if the specified object is a
     * {@code Pair} instance representing the same
     * ordered pair of objects as {@code this}, and
     * returns {@code false} otherwise.
     */
    @Override
    public boolean equals(Object obj) {
        if (obj==this) {
            return true;
        }
        if ((obj instanceof Pair)==false) {
            return false;
        }
        Pair p = (Pair) obj;
        return (first.equals(p.first) && second.equals(p.second));
    }

    /**
     * Returns a string representation of {@code this}. The
     * exact details of the representation are unspecified and
     * subject to change.
     *
     * @return a string representation of {@code this}.
     */
    @Override
    public String toString() {
        return "[" + first + ", " + second + "]";
    }
}
