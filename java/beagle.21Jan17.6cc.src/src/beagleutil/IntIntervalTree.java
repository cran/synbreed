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

import java.util.Collection;

/**
 * Interface {@code IntIntervalTree} represents an interval
 * tree whose elements are {@code IntInterval} objects.
 *
 * @param <E> the type of objected stored in {@code this}
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}

 */
public interface IntIntervalTree<E extends IntInterval> {

    /**
     * Returns the minimum start (inclusive) of an interval
     * that can be stored in this interval tree.
     * @return the minimum start (inclusive) of an interval
     * that can be stored in this interval tree
     */
    public int start();

    /**
     * Returns the maximum end (inclusive) of an interval
     * that can be stored in this interval tree.
     * @return the maximum end (inclusive) of an interval
     * that can be stored in this interval tree
     */
    public int end();

    /**
     * Removes all of the elements from this interval tree.
     */
    public void clear();

    /**
     * Adds the specified element to this interval tree, and returns
     * {@code true} if the interval tree is changed as a result of
     * the call. The method returns {@code false} if
     * {@code this.contains(E) == true} when the method is invoked.
     *
     * @param element the element to be added
     * @return {@code true} if the interval tree changed as
     * a result of the call
     *
     * @throws IllegalArgumentException if
     * {@code element.start() < this.start() || element.end() > this.end()}
     * @throws NullPointerException if {@code element == null}
     */
    public boolean add(E element);

    /**
     * Returns {@code true} if the interval tree contains the specified
     * element, and returns {@code false} otherwise.
     * @param element the element whose presence in the interval tree
     * is to be tested
     * @return {@code true} if the interval tree contains the specified
     * element
     * @throws NullPointerException if {@code element == null}
     */
    public boolean contains(E element);

    /**
     * Removes the specified element from this interval tree if the
     * specified element is found in the interval tree.
     * @param element the element to be removed from this interval tree
     * @return {@code true} if the interval tree is changed as
     * a result of the call
     * @throws NullPointerException if {@code element == null}
     */
    public boolean remove(E element);

    /**
     * Adds the elements in this interval tree that intersect the specified
     * point to the specified collection.
     *
     * @param point a point
     * @param collection a collection to which will be added the elements of
     * this interval tree that intersect the specified point
     *
     * @throws NullPointerException if {@code collection == null}
     */
    public void intersect(int point, Collection<E> collection);

    /**
     * Adds the elements in this interval tree that intersect any part of
     * the specified interval to the specified collection.
     *
     * @param start the start (inclusive) of the specified interval
     * @param end the end (inclusive) of the specified interval
     * @param collection a collection to which will be added the elements of
     * this interval tree that intersect any part of the specified interval
     *
     * @throws NullPointerException if {@code collection == null}.
     */
    public void intersectPart(int start, int end, Collection<E> collection);

    /**
     * Adds the elements in this interval tree that contain
     * the specified interval to the specified collection.
     *
     * @param start the start (inclusive) of the specified interval
     * @param end the end (inclusive) of the specified interval
     * @param collection a collection to which will be added the elements
     * of this interval tree that contain the specified interval
     *
     * @throws NullPointerException if {@code collection == null}
     */
    public void intersectAll(int start, int end, Collection<E> collection);

    /**
     * Returns {@code true} if this interval tree contains no elements.
     * @return {@code true} if this interval tree contains no elements
     */
    public boolean isEmpty();

    /**
     * Returns the number of elements in this interval tree.
     *
     * @return the number of elements in this interval tree
     */
    public int size();

    /**
     * Returns an array containing all of the elements of this interval tree.
     * @return an array containing all of the elements of this interval tree
     */
    public E[] toArray();
}
