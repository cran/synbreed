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
package sample;

import java.util.NoSuchElementException;

/**
 * <p>Class {@code DiploidStates} represents a list of iterators
 * (one iterator for each marker) that iterate over a subset of diploid
 * HMM states at a marker.
 * </p>
 * <p>Instances of class {@code DiploidStates} are not requires to be
 * thread-safe.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface DiploidStates {

    /**
     * Returns the number of markers.
     * @return the number of markers
     */
    public int nMarkers();

   /**
    * Initializes the iteration of permitted ordered edge pairs for the
    * specified marker.
    * @param marker a marker index
    * @throws IllegalArgumentException if
    * {@code marker < 0 || marker >= this.nMarkers()}
    */
    public void setMarker(int marker);

    /**
     * Returns the current marker index.
     * @return the current marker index
     */
    public int marker();

    /**
     * Returns {@code true} if the iteration of the ordered edge pairs has
     * more elements, and returns {@code false} otherwise.
     * @return {@code true} if the iteration of the ordered edge pairs has
     * more elements
     */
    public boolean hasNext();

    /**
     * Advances the iteration of ordered edge pairs to the next element.
     *
     * @throws NoSuchElementException if {@code this.hasNext() == false}
     */
    public void next();

    /**
     * Returns the first edge of the edge pair that is the current element
     * in the iterations, or returns {@code -1} if {@code this.next()}
     * has not been invoked since the most recent invocation of
     * {@code this.setMarker()}.
     * @return the first edge of the edge pair that is the current element
     * in the iterations
     */
    public int edge1();

    /**
     * Returns the second edge of the edge pair that is the current element
     * in the iteration, or returns {@code -1} if {@code this.next()}
     * has not been invoked since the most recent invocation of
     * {@code this.setMarker()}.
     * @return the second edge of the edge pair that is
     * the current element in the iteration
     */
    public int edge2();
}