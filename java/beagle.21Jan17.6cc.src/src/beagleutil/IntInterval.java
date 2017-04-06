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
package beagleutil;

/**
 * <p>Interface {@code IntInterval} represents an interval of
 * consecutive integers.
 * </p>
 * Instances of class {@code IntInterval} are immutable.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface IntInterval {

    /**
     * Returns the start of the interval (inclusive).
     * @return the start of the interval (inclusive).
     */
    public int start();

    /**
     * Returns the end of the interval (inclusive).
     * @return the end of the interval (inclusive).
     */
    public int end();
}
