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
 * Class {@code FilterUtils} contains static methods for constructing
 * filters.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class FilterUtils {

    private FilterUtils() {
        // private constructor to prevent instantiation
    }

    /**
     * Returns a filter that accepts all non-null objects.
     * @param <E> the type of object that is filtered.
     * @return a filter that accepts all non-null objects.
     */
    public static <E> Filter<E> acceptAllFilter() {
        return new Filter<E>() {
            @Override
            public boolean accept(E e) {
                if (e==null) {
                    throw new NullPointerException("e==null");
                }
                return true;
            }
        };
    }

    /**
     * Returns a filter that accepts only objects that are equal
     * to the specified object.
     * @param <E> the type of object that is filtered.
     * @param singleton the object that will be accepted.
     * @return a filter that accepts only objects that are equal
     * to the specified object.
     * @throws NullPointerException if {@code singleton==null}
     */
    public static <E> Filter<E> singletonFilter(final E singleton) {
        if (singleton==null) {
            throw new NullPointerException("singleton==null");
        }
        return new Filter<E>() {
            @Override
            public boolean accept(E e) {
                if (e==null) {
                    throw new NullPointerException("e==null");
                }
                return singleton.equals(e);
            }
        };
    }

    /**
     * Returns a filter that accepts all non-null objects in the
     * specified collection.
     * @param <E> the type of object that is filtered.
     * @param include the collection of elements that will be accepted by
     * the filter.
     * @return a filter that accepts all non-null objects in the
     * specified collection.
     * @throws NullPointerException if {@code include==null}
     */
    public static <E> Filter<E> includeFilter(Collection<E> include) {
        final Set<E> includeSet = new HashSet<>(include);
        return new Filter<E>() {
            @Override
            public boolean accept(E e) {
                if (e==null) {
                    throw new NullPointerException("e==null");
                }
                return includeSet.contains(e);
            }
        };
    }

    /**
     * Returns a filter that accepts all non-null objects that are not
     * contained in the specified collection.
     * @param <E> the type of object that is filtered.
     * @param exclude the collection of elements that will be rejected
     * by the filter.
     * @return a filter that accepts all non-null objects that are not
     * contained in the specified collection.
     * @throws NullPointerException if {@code exclude==null}
     */
    public static <E> Filter<E> excludeFilter(Collection<E> exclude) {
        final Set<E> includeSet = new HashSet<>(exclude);
        return new Filter<E>() {
            @Override
            public boolean accept(E e) {
                if (e==null) {
                    throw new NullPointerException("e==null");
                }
                return !includeSet.contains(e);
            }
        };
    }

    /**
     * Returns a filter that accepts all non-null objects that are contained
     * in the specified {@code include} collection and that are not in the
     * contained in the specified {@code exclude} collection.
     * @param <E> the type of object that is filtered.
     * @param include the collection of elements that may be accepted
     * by the filter.
     * @param exclude the collection of elements that will be rejected
     * by the filter.
     * @return a filter that accepts all non-null objects that are contained
     * in the specified {@code include} collection and that are not in the
     * contained in the specified {@code exclude} collection.
     * @throws NullPointerException if {@code include==null || exclude==null}
     */
    public static <E> Filter<E> includeExcludeFilter(final Collection<E> include,
            final Collection<E> exclude) {
        final Set<E> includeSet = new HashSet<>(include);
        final Set<E> excludeSet = new HashSet<>(exclude);
        return new Filter<E>() {
            @Override
            public boolean accept(E e) {
                if (e==null) {
                    throw new NullPointerException("e==null");
                }
                return includeSet.contains(e) && !excludeSet.contains(e);
            }
        };
    }
}
