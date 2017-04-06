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

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

/**
 * <p>A filter for accepting or rejecting objects.
 * </p>
 * Instances of class {@code Filter} are required to be immutable.
 *
 * @param <E> the type of object that is filtered.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface Filter<E> {

    /**
     * Returns a filter that accepts all non-null objects.
     * @param <E> the type of object that is filtered
     * @return a filter that accepts all non-null objects
     */
    static <E> Filter<E> acceptAllFilter() {
        return (E e) -> {
            if (e == null) {
                throw new NullPointerException("e==null");
            }
            return true;
        };
    }

    /**
     * Returns a filter that accepts all non-null objects that are
     * contained in the specified collection.
     * @param <E> the type of object that is filtered
     * @param include the collection of objects that will be accepted by
     * the filter
     * @return a filter that accepts all non-null objects that are
     * contained in the specified collection
     * @throws NullPointerException if {@code include == null}
     */
    static <E> Filter<E> includeFilter(Collection<E> include) {
        final Set<E> includeSet = new HashSet<>(include);
        return (E e) -> {
            if (e == null) {
                throw new NullPointerException("e==null");
            }
            return includeSet.contains(e);
        };
    }

    /**
     * Returns a filter that accepts all non-null objects that are not
     * contained in the specified collection.
     * @param <E> the type of object that is filtered
     * @param exclude the collection of objects that will be rejected
     * by the filter
     * @return a filter that accepts all non-null objects that are not
     * contained in the specified collection
     * @throws NullPointerException if {@code exclude == null}
     */
    static <E> Filter<E> excludeFilter(Collection<E> exclude) {
        final Set<E> includeSet = new HashSet<>(exclude);
        return (E e) -> {
            if (e == null) {
                throw new NullPointerException("e==null");
            }
            return !includeSet.contains(e);
        };
    }

    /**
     * Returns {@code true} if the specified object is
     * accepted and returns {@code false} if the specified object
     * is rejected.
     * @param e the object to be filtered
     * @return {@code true} if the specified object is
     * accepted
     * @throws NullPointerException if {@code e==null}
     */
    boolean accept(E e);
}
