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

/**
 * <p>Class {@code Score} represents a similarity score for a pair
 * of trees.
 * </p>
 * Instances of class {@code Score} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class Score implements Comparable<Score> {

    private final int nodeA;
    private final int nodeB;
    private final float score;

    /**
     * Constructs a new {@code Score} instance.  Smaller similarity scores
     * correspond to greater similarity.
     * @param nodeA root node index for the first tree
     * @param nodeB root node index for the second tree
     * @param score the a non-negative similarity score for the two specified
     * trees
     * @param isMergeable {@code true} if the two trees may be
     * merged, and {@code false} otherwise
     * @throws IllegalArgumentException if
     * {@code score < 0 || (score==0 && isMergeable==false)}
     * @throws IllegalArgumentException if {@code Float.isNaN(score)}
     */
    public Score(int nodeA, int nodeB, float score, boolean isMergeable) {
        if (score < 0 || (score==0 && isMergeable==false) || Float.isNaN(score)) {
            throw new IllegalArgumentException(String.valueOf(score));
        }
        this.nodeA = nodeA;
        this.nodeB = nodeB;
        this.score = isMergeable ? score : -score;
    }

    /**
     * Returns the root node index for the first tree.
     *
     * @return the root node index for the first tree
     */
    public int nodeA() {
        return nodeA;
    }

    /**
     * Returns the root node index for the second tree.
     *
     * @return the root node index for the second tree
     */
    public int nodeB() {
        return nodeB;
    }

    /**
     * Returns the similarity score for the two trees.
     *
     * @return the similarity score for the two trees
     */
    public float score() {
        return (score < 0) ? -score : score;
    }

    /**
     * Returns {@code true} if the two trees may be merged, and
     * returns {@code false} otherwise.
     *
     * @return {@code true} if the two trees may be merged, and
     * returns {@code false} otherwise
     */
    public boolean isMergeable() {
        return score>=0;
    }

    /**
     * Compares the specified object with this {@code Score} for
     * equality.  Returns {@code  true} if the specified object
     * is a {@code Score} instance whose {@code nodeA()}, {@code nodeB()},
     * {@code score()}, and {@code isMergeable()} methods return the
     * same values as the corresponding methods for {@code  this}, and
     * returns {@code  false} otherwise.
     * @param obj the object to be compared for equality with this
     * {@code  Score}
     * @return {@code  true} if the specified object is a equal to {@code this}
     */
    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final Score other = (Score) obj;
        if (this.nodeA != other.nodeA) {
            return false;
        }
        if (this.nodeB != other.nodeB) {
            return false;
        }
        return Float.floatToIntBits(this.score)
                == Float.floatToIntBits(other.score);
    }

     /**
     * <p>Returns the hash code value for this object. The hash code is defined
     * by the following calculation:
     * </p>
     * <pre>
     *  int hash = 5;
     *  hash = 53 * hash + this.nodeA();
     *  hash = 53 * hash + this.nodeB();
     *  hash = 53 * hash + Float.floatToIntBits(this.score());
     </pre>
     * @return a hash code value for the object
     */
    @Override
    public int hashCode() {
        int hash = 5;
        hash = 53 * hash + this.nodeA;
        hash = 53 * hash + this.nodeB;
        hash = 53 * hash + Float.floatToIntBits(this.score);
        return hash;
    }

    /**
     * Returns -1, 0, or 1 depending on whether this {@code Score} is less
     * than, equal to, or greater than the specified {@code Score}.  The
     * two scores are ordered first using {@code -Boolean.compare()} on
     * the value returned by {@code isMergeable()}, then by
     * {@code Float.compare()} on the value returned by {@code score()}, then
     * by the value returned by {@code nodeA()}, and then by the value returned
     * by {@code nodeB()}.
     * @param other a {@code Score} element to be compared
     * @return a negative integer, zero, or a positive integer depending
     * on whether this {@code Score} is less than, equal to, or greater
     * than the specified {@code Score}
     *
     * @throws NullPointerException if {@code other == null}
     */
    @Override
    public int compareTo(Score other) {
        int x = -Boolean.compare(this.isMergeable(), other.isMergeable());
        if (x!=0) {
            return x;
        }
        x = Float.compare(this.score(), other.score());
        if (x!=0) {
            return x;
        }
        if (this.nodeA() != other.nodeA()) {
            return (this.nodeA() < other.nodeA()) ? -1 : 1;
        }
        if (this.nodeB() != other.nodeB()) {
            return (this.nodeB() < other.nodeB()) ? -1 : 1;
        }
        return 0;
    }

    /**
     * Returns a string representation of {@code this}. The exact details of
     * the representation are unspecified and subject to change.
     *
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("[nodeA=");
        sb.append(nodeA());
        sb.append(", nodeB=");
        sb.append(nodeB());
        sb.append(", score=");
        sb.append(score());
        sb.append(" isMergeable=");
        sb.append(isMergeable());
        sb.append("]");
        return sb.toString();
    }
}
