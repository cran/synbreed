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
 * <p>Class {@code Phase} represents the equivalence of two phased genotypes
 * for a marker or for a set of markers.  Genotype equivalence is defined
 * in terms of allele equivalence.  Two alleles are equivalent if
 * either allele is missing or if both alleles are non-missing and equal.
 * </p>
 * <p>
 * For the case of a single marker with phased (i.e. ordered) genotypes
 * ({@code a1}, {@code a2})
 * and ({@code b1}, {@code b2}), then
 * <br>
 * 1) the genotypes have IDENTICAL phase if a) alleles {@code a1}
 * and {@code b1} are equivalent, b) alleles {@code a2} and
 * {@code b2} are equivalent, and c) either alleles {@code a1}
 * and {@code b2} are not equivalent or alleles
 * {@code a2} and {@code b1} are not equivalent.
 * <br>
 * 2) the genotypes have OPPOSITE phase if a) alleles {@code a1}
 * and {@code b2} are equivalent, b) alleles {@code a2} and
 * {@code b1} are equivalent, and c) either alleles {@code a1}
 * and {@code b1} are not equivalent or alleles {@code a2} and
 * {@code b2} are not equivalent.
 * <br>
 * 3) the genotypes have UNKOWN phase if a) alleles {@code a1}
 * and {@code b1} are equivalent, b) alleles {@code a2} and
 * {@code b2} are equivalent, c) alleles {@code a1} and
 * {@code b2} are equivalent, and d) alleles {@code a2} and
 * {@code b1} are equivalent.
 * <br>
 * 4) the genotypes have INCONSISTENT phase if a) either alleles
 * {@code a1} and {@code b1} are not equivalent or alleles
 * {@code a2} and {@code b2} are not equivalent, and
 * b) either alleles {@code a1} and {@code b2} are not equivalent
 * or alleles {@code a2} and {@code b1} are not equivalent.
 * </p>
 * For the case of two sets of phased genotypes for the same markers,
 * the two sets have
 * <br>
 * 1) IDENTICAL phase if the phase is IDENTICAL for at least
 * one marker and is either IDENTICAL or UNKNOWN for all markers.
 * <br>
 * 2) OPPOSITE phase if the if the phase is OPPOSITE for at
 * lease one marker and is either OPPOSITE or UNKNOWN for all markers.
 * <br>
 * 3) UNKNOWN phase if the phase is UNKNOWN for all markers.
 * <br>
 * 4) INCONSISTENT phase if a) the phase is INCONSISTENT for at least one
 * marker or if b) the relative phase is IDENTICAL for at least one marker and
 * OPPOSITE for at least one marker.
 *
 * @author Brian L. Browning
 */
public enum Phase {
        IDENTICAL,
        OPPOSITE,
        UNKNOWN,
        INCONSISTENT
}
