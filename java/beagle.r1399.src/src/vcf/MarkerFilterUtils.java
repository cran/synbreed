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
package vcf;

import blbutil.Filter;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

/**
 * Class {@code FilterUtils} contains static methods for constructing
 * marker filters.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class MarkerFilterUtils {

    /**
     * Returns a filter that accepts all non-null markers which do not have
     * an identifier in the specified collection.
     * A marker is excluded if {@code exclude.contains(marker.id(j))==true}
     * for any {@code 0 <= j < marker.nIds()} or if
     * {@code exclude.contains(marker.chrom() + ":" + marker.pos())==true}.
     * @param exclude the marker identifiers that will be rejected
     * by the filter.
     * @return a filter that accepts all non-null markers which do not have
     * an identifier in the specified collection.
     * @throws NullPointerException if {@code exclude==null}
     */
    public static Filter<Marker> excludeIdFilter(Collection<String> exclude) {
        final Set<String> excludeSet = new HashSet<>(exclude);
        return new Filter<Marker>() {
            @Override
            public boolean accept(Marker marker) {
                boolean accept = true;
                int n = marker.nIds();
                for (int j=0; j<n && accept==true; ++j) {
                    if (excludeSet.contains(marker.id(j))) {
                        accept = false;
                    }
                }
                if (accept==true) {
                    String posId = marker.chrom() + ':' + marker.pos();
                    if (excludeSet.contains(posId)) {
                        accept = false;
                    }
                }
                return accept;
            }
        };
    }
}
