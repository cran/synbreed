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
 * A filter for accepting or rejecting objects.
 *
  @param <E> the type of object that is filtered.
 *
 * @author Brian L. Browning
 */
public interface Filter<E> {

    /**
     * Returns {@code true} if the specified object is
     * accepted and returns {@code false} otherwise.
     * @param e the object to be filtered.
     * @return {@code true} if the specified object is
     * accepted and {@code false} if the specified object is
     * rejected.
     * @throws NullPointerException if {@code e==null}.
     */
    public abstract boolean accept(E e);
}
