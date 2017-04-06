package beagleutil;

import blbutil.Const;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;


/**
 * <p>Class {@code CenteredIntIntervalTree} implements a centered
 * interval tree that stores {@code IntInterval} objects.
 * </p>
 * Instances of class {@code CenteredIntIntervalTree} are not thread-safe.
 * @param <E> the objects stored by {@code this}.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class CenteredIntIntervalTree<E extends IntInterval & Comparable<E>>
        implements IntIntervalTree<E> {

    private final int start;
    private final int end;
    private int size;
    private Node<E> root;

    /**
     * Creates a new {@code CenteredIntIntervalTree} instance for the
     * specified range.
     * @param start the minimum start value (inclusive) for intervals stored in
     * this interval tree
     * @param end the maximum end value (inclusive) for intervals stored
     * in this interval tree
     *
     * @throws IllegalArgumentException if {@code end < start}
     */
    public CenteredIntIntervalTree(int start, int end) {
        if (end < start) {
            String s = "end= " + end + " start=" + start;
            throw new IllegalArgumentException(s);
        }
        int length = (end - start + 1);
        int center = start + (length/2);
        this.start = start;
        this.end = end;
        this.size = 0;
        this.root = new Node<>(center);
    }

    @Override
    public int start() {
        return start;
    }

    @Override
    public int end() {
        return end;
    }

    @Override
    public void clear() {
        clear(root);
        size = 0;
    }

    private void clear(Node<E> tree) {
        if (tree==null) {
            return;
        }
        tree.clear();
        clear(tree.leftChild);
        clear(tree.rightChild);
    }

    @Override
    public boolean add(E element) {
        if (element.start() < start || element.end() > end) {
            String s = "element out of range [" + start + ", " + end + ") : "
                    + element;
            throw new IllegalArgumentException(s);
        }
        boolean added = add(root, element);
        if (added) {
            ++size;
        }
        return added;
    }

    private boolean add(Node<E> tree, E element) {
        if (element.end() < tree.center) {
            if (tree.leftChild==null) {
                int nextOffset = nextOffset(tree);
                tree.leftChild = new Node<>(tree.center - nextOffset);
                tree.leftChild.parent = tree;
            }
            return add(tree.leftChild, element);
        }
        else if (element.start() > tree.center) {
            if (tree.rightChild==null) {
                int nextOffset = nextOffset(tree);
                tree.rightChild = new Node<>(tree.center + nextOffset);
                tree.rightChild.parent = tree;
            }
            return add(tree.rightChild, element);
        }
        else {
            return tree.add(element);
        }
    }

    private int nextOffset(Node<E> node) {
        int lastOffset;
        if (node.parent==null) {
            lastOffset = (end - start + 1)/2;
        }
        else {
            lastOffset = Math.abs(node.center - node.parent.center);
        }
        assert lastOffset > 0;
        int offset = (lastOffset+1)/2;
        return offset;
    }

    @Override
    public boolean contains(E element) {
        return contains(root, element);
    }

    private boolean contains(Node<E> tree, E element) {
        if (tree==null) {
            return false;
        }
        else if (element.end() < tree.center) {
            return contains(tree.leftChild, element);
        }
        else if (element.start() > tree.center) {
            return contains(tree.rightChild, element);
        }
        else {
            return tree.contains(element);
        }
    }

    @Override
    public boolean remove(E element) {
        boolean removed = remove(root, element);
        if (removed) {
            --size;
        }
        return removed;
    }

    private boolean remove(Node<E> tree, E element) {
        if (tree==null) {
            return false;
        }
        else if (element.end() < tree.center) {
            return remove(tree.leftChild, element);
        }
        else if (element.start() > tree.center) {
            return remove(tree.rightChild, element);
        }
        else {
            return tree.remove(element);
        }
    }

    @Override
    public void intersect(final int point, Collection<E> collection) {
        intersect(root, point, collection);
    }

    private void intersect(Node<E> tree, int point, Collection<E> collection) {
        if (tree==null) {
            return;
        }
        tree.intersect(point, collection);
        if (point < tree.center) {
            intersect(tree.leftChild, point, collection);
        }
        else if (point > tree.center) {
            intersect(tree.rightChild, point, collection);
        }
    }

    @Override
    public void intersectPart(int start, int end, Collection<E> collection) {
        intersectPart(root, start, end, collection);
    }

    private void intersectPart(Node<E> tree, int start, int end,
            Collection<E> collection) {
        if (tree==null) {
            return;
        }
        tree.intersectPart(start, end, collection);
        if (start < tree.center) {
            intersectPart(tree.leftChild, start, end, collection);
        }
        if (end > tree.center) {
            intersectPart(tree.rightChild, start, end, collection);
        }
    }

    @Override
    public void intersectAll(int start, int end, Collection<E> collection) {
        intersectAll(root, start, end, collection);
    }

    private void intersectAll(Node<E> tree, int start, int end,
            Collection<E> collection) {
        if (tree==null) {
            return;
        }
        tree.intersectAll(start, end, collection);
        if (end < tree.center) {
            intersectAll(tree.leftChild, start, end, collection);
        }
        if (start > tree.center) {
            intersectAll(tree.rightChild, start, end, collection);
        }
    }

    @Override
    public boolean isEmpty() {
        return size==0;
    }

    @Override
    public int size() {
        return size;
    }

    @Override
    public E[] toArray() {
        List<E> list = new ArrayList<>(size);
        toArray(root, list);
        return (E[]) list.toArray();
    }

    private void toArray(Node<E> tree, List<E> list) {
        if (tree==null) {
            return;
        }
        toArray(tree.leftChild, list);
        list.addAll(tree.sortedStart);
        toArray(tree.rightChild, list);
    }

    /**
     * Returns a string representation of {@code this}.  The
     * exact details of the representation are unspecified and
     * subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("[ CenteredIntIntervalTree: ");
        sb.append(Const.nl);
        sb.append("start=");
        sb.append(start);
        sb.append(" end=");
        sb.append(end);
        sb.append(" size=");
        sb.append(size);
        sb.append(Const.nl);
        toString(root, sb);
        sb.append(']');
        return sb.toString();
    }

    private void toString(Node<E> tree, StringBuilder sb) {
        if (tree==null) {
            return;
        }
        toString(tree.leftChild, sb);
        sb.append(tree);
        toString(tree.rightChild, sb);
    }

    private static final class Node<E extends IntInterval & Comparable<E>> {

        private final int center;
        private final SortedSet<E> sortedStart;
        private final SortedSet<E> sortedEnd;
        private Node<E> parent;
        private Node<E> leftChild;
        private Node<E> rightChild;

        /**
         * Returns a {@code Comparator} that is consistent with equals and
         * orders elements in order of increasing start values.
         * @return a {@code Comparator} that is consistent with equals and
         * orders elements in order of increasing start values
         */
        private static <E extends IntInterval & Comparable<E>> Comparator<E> startComparator() {
            return (E e1, E e2) -> {
                int start1 = e1.start();
                int start2 = e2.start();
                if (start1 == start2) {
                    return e1.compareTo(e2);
                }
                else {
                    return (start1 < start2) ? -1 : 1;
                }
            } ;
        }

        /**
         * Returns a {@code Comparator} that is consistent with equals and
         * orders elements in order of decreasing end values.
         * @return a {@code Comparator} that is consistent with equals and
         * orders elements in order of decreasing end values
         */
        private static <E extends IntInterval & Comparable<E>> Comparator<E> endComparator() {
            return (E e1, E e2) -> {
                int end1 = e1.end();
                int end2 = e2.end();
                if (end1 == end2) {
                    return e1.compareTo(e2);
                }
                else {
                    return (end1 > end2) ? -1 : 1;
                }
            } ;
        }

        Node(int center) {
            this.center = center;
            Comparator<E> startComparator = startComparator();
            Comparator<E> endComparator = endComparator();
            this.sortedStart = new TreeSet<>(startComparator);
            this.sortedEnd = new TreeSet<>(endComparator);
            this.leftChild = null;
            this.rightChild = null;
        }

        boolean add(E element) {
            if (element.start() > center || element.end() < center) {
                String s = "element does not overlap center=" + center + ": "
                        + element;
                throw new IllegalArgumentException(s);
            }
            boolean addedStart = sortedStart.add(element);
            boolean addedEnd = sortedEnd.add(element);
            assert addedStart == addedEnd;
            return addedStart;
        }

        boolean contains(E element) {
            boolean startContains = sortedStart.contains(element);
            assert startContains == sortedEnd.contains(element);
            return startContains;
        }

        boolean remove(E element) {
            boolean removedStart = sortedStart.remove(element);
            boolean removedEnd = sortedEnd.remove(element);
            assert removedStart == removedEnd;
            return removedStart;
        }

        void intersect(int point, Collection<E> collection) {
            if (point <= center) {
                boolean finished = false;
                Iterator<E> it = sortedStart.iterator();
                while (it.hasNext() && finished==false) {
                    E e = it.next();
                    if (e.start() <= point) {
                        collection.add(e);
                    }
                    else {
                        finished = true;
                    }
                }
            }
            else {
                boolean finished = false;
                Iterator<E> it = sortedEnd.iterator();
                while (it.hasNext() && finished==false) {
                    E e = it.next();
                    if (e.end() >= point) {
                        collection.add(e);
                    }
                    else {
                        finished = true;
                    }
                }
            }
        }

        void intersectPart(int start, int end, Collection<E> collection) {
            if (end < center) {
                boolean finished = false;
                Iterator<E> it = sortedStart.iterator();
                while (it.hasNext() && finished==false) {
                    E e = it.next();
                    if (e.start() <= end) {
                        collection.add(e);
                    }
                    else {
                        finished = true;
                    }
                }
            }
            else if (start > center) {
                boolean finished = false;
                Iterator<E> it = sortedEnd.iterator();
                while (it.hasNext() && finished==false) {
                    E e = it.next();
                    if (start <= e.end()) {
                        collection.add(e);
                    }
                    else {
                        finished = true;
                    }
                }
            }
            else {
                collection.addAll(sortedStart);
            }
        }

        void intersectAll(int start, int end, Collection<E> collection) {
            boolean finished = false;
            Iterator<E> it = sortedStart.iterator();
            while (it.hasNext() && finished==false) {
                E e = it.next();
                if (e.start() <= start) {
                    if (e.end() >= end) {
                        collection.add(e);
                    }
                }
                else {
                    finished = true;
                }
            }
        }

        void clear() {
            sortedStart.clear();
            sortedEnd.clear();
        }

        /**
         * Returns a string representation of {@code this}.  The
         * exact details of the representation are unspecified and
         * subject to change.
         * @return a string representation of {@code this}
         */
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append(Const.nl);
            sb.append("[ CenteredIntIntervalTree.Node:");
            sb.append(Const.nl);
            sb.append("center=");
            sb.append(center);
            sb.append(" parent.center=");
            sb.append(parent!=null ? parent.center : null);
            sb.append(" leftchild.center=");
            sb.append(leftChild!=null ? leftChild.center : null);
            sb.append(" rightchild.center=");
            sb.append(rightChild!=null ? rightChild.center : null);
            sb.append(Const.nl);
            sb.append("sortedStart: ");
            sb.append(sortedStart);
            sb.append(Const.nl);
            sb.append("sortedEnd: ");
            sb.append(sortedEnd);
            sb.append(Const.nl);
            sb.append(']');
            return sb.toString();
        }
    }

    //<editor-fold defaultstate="collapsed" desc="code for testing class">
    /*
     * The main() method if for testing the CenteredIntIntervalTree class
     */
//    public static void main(String[] args) {
//        main1(args);
//        main2(args);
//    }

    private static void main1(String[] args) {
        int length=16;
        IntIntervalTree<BasicIntInterval> tree = new CenteredIntIntervalTree<>(0,length);
        assert tree.start()==0;
        assert tree.end()==length;
        assert tree.isEmpty();
        assert tree.isEmpty();
        for (int j=0; j<length; ++j) {
            BasicIntInterval i = new BasicIntInterval(j, j+1);
            assert tree.contains(i)==false;
            boolean added = tree.add(i);
            assert added == true;
            assert tree.contains(i)==true;
            added = tree.add(i);
            assert added == false;
        }
        assert tree.size()==length;
        System.out.println("Initial Tree: " + java.util.Arrays.toString(tree.toArray()));
        for (int j=0; j<length; j+=2) {
            BasicIntInterval i = new BasicIntInterval(j, j+1);
            assert tree.contains(i)==true;
            boolean removed = tree.remove(i);
            assert removed == true;
            assert tree.contains(i)==false;
            removed = tree.remove(i);
            assert removed==false;
        }
        assert tree.size()==(length/2);
        System.out.println("Pruned Tree: " + java.util.Arrays.toString(tree.toArray()));

        List<BasicIntInterval> list = new ArrayList<>(length);
        for (int j=0; j<length; ++j) {
            tree.intersect(j, list);
            System.out.println("point=" + j + ": " + list);
            list.clear();
        }

        int intSize = 3;
        for (int j=0; j<length; ++j) {
            tree.intersectPart(j, j+intSize, list);
            System.out.println("start=" + j + " end=" + (j+intSize) + ": " + list);
            list.clear();
        }

        for (int j=0; j<length; ++j) {
            tree.intersectAll(j, j, list);
            System.out.println("start=" + j + " end=" + j + ": " + list);
            list.clear();
        }

        tree.clear();
        assert tree.isEmpty()==true;
        System.out.println("Empty Tree: " + java.util.Arrays.toString(tree.toArray()));
    }

    private static void main2(String[] args) {
        int length=16;
        int nOverlaps = 4;
        IntIntervalTree<BasicIntInterval> tree = new CenteredIntIntervalTree<>(-length,length);
        assert tree.start()==-length;
        assert tree.end()==length;
        assert tree.isEmpty();
        assert tree.isEmpty();
        for (int j=1; j<=nOverlaps; ++j) {
            BasicIntInterval i = new BasicIntInterval(-j - length/2, -j);
            assert tree.contains(i)==false;
            boolean added = tree.add(i);
            assert added == true;
            assert tree.contains(i)==true;
            added = tree.add(i);
            assert added == false;
        }
        assert tree.size()==nOverlaps;
        for (int j=1; j<=nOverlaps; ++j) {
            BasicIntInterval i = new BasicIntInterval(j, j + length/2);
            assert tree.contains(i)==false;
            boolean added = tree.add(i);
            assert added == true;
            assert tree.contains(i)==true;
            added = tree.add(i);
            assert added == false;
        }
        assert tree.size()==2*nOverlaps;
        System.out.println(java.util.Arrays.toString(tree.toArray()));
        System.out.println(tree);

        for (int j=1; j<=nOverlaps; j+=2) {
            BasicIntInterval i = new BasicIntInterval(-j - length/2, -j);
            assert tree.contains(i)==true;
            boolean removed = tree.remove(i);
            assert removed == true;
            assert tree.contains(i)==false;
            removed = tree.remove(i);
            assert removed==false;
        }
        for (int j=1; j<=nOverlaps; j+=2) {
            BasicIntInterval i = new BasicIntInterval(j, j + length/2);
            assert tree.contains(i)==true;
            boolean removed = tree.remove(i);
            assert removed == true;
            assert tree.contains(i)==false;
            removed = tree.remove(i);
            assert removed==false;
        }
        assert tree.size()==nOverlaps;
        System.out.println(java.util.Arrays.toString(tree.toArray()));

        List<BasicIntInterval> list = new ArrayList<>(length);
        for (int j=-length; j<length; ++j) {
            tree.intersect(j, list);
            System.out.println("point=" + j + ": " + list);
            list.clear();
        }

        int intSize = 3;
        for (int j=-length; j<length; ++j) {
            tree.intersectPart(j, j+intSize, list);
            System.out.println("start=" + j + " end=" + (j+intSize) + ": " + list);
            list.clear();
        }

        tree.clear();
        assert tree.isEmpty()==true;
        System.out.println(java.util.Arrays.toString(tree.toArray()));
    }
    //</editor-fold>
}
