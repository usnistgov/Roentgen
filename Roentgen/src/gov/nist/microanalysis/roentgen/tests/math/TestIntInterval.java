package gov.nist.microanalysis.roentgen.tests.math;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.HashSet;
import java.util.Set;

import org.junit.Test;

import gov.nist.microanalysis.roentgen.math.IntInterval;

/**
 * <p>
 * Description
 * </p>
 *
 * @author Nicholas
 * @version 1.0
 */
public class TestIntInterval {

   /**
    * Constructs a TestIntInterval
    */
   @Test
   public void testIntersects() {
      final IntInterval i1 = new IntInterval(1, 3);
      final IntInterval i2 = new IntInterval(4, 8);
      final IntInterval i3 = new IntInterval(7, 9);
      final IntInterval i4 = new IntInterval(11, 13);
      final IntInterval i5 = new IntInterval(1, 13);
      assertFalse(i1.intersects(i2));
      assertTrue(i2.intersects(i3));
      assertFalse(i2.intersects(i1));
      assertTrue(i3.intersects(i2));

      assertTrue(i1.intersects(i5));
      assertTrue(i3.intersects(i5));
      assertFalse(i3.intersects(i4));
      assertFalse(i3.intersects(IntInterval.NULL_INTERVAL));
      assertFalse(IntInterval.NULL_INTERVAL.intersects(i5));
   }

   @Test
   public void testExtent() {
      final IntInterval i1 = new IntInterval(1, 3);
      final IntInterval i2 = new IntInterval(4, 8);
      final IntInterval i3 = new IntInterval(7, 9);
      final IntInterval i4 = new IntInterval(11, 13);
      final IntInterval i5 = new IntInterval(1, 13);
      IntInterval.extent(i1, i2).equals(new IntInterval(1, 8));
      IntInterval.extent(i1, i4).equals(i5);
      IntInterval.extent(i2, i1).equals(new IntInterval(1, 8));
      IntInterval.extent(i3, i2).equals(new IntInterval(4, 9));
   }

   @Test
   public void testIntersection() {
      final IntInterval i1 = new IntInterval(1, 3);
      final IntInterval i2 = new IntInterval(4, 8);
      final IntInterval i3 = new IntInterval(7, 9);
      final IntInterval i4 = new IntInterval(11, 13);
      final IntInterval i5 = new IntInterval(1, 13);
      IntInterval.intersection(i1, i2).equals(IntInterval.NULL_INTERVAL);
      IntInterval.extent(i2, i3).equals(new IntInterval(7, 8));
      IntInterval.extent(i2, i1).equals(new IntInterval(1, 8));
      IntInterval.extent(i2, i5).equals(i2);
      IntInterval.extent(i4, i3).equals(IntInterval.NULL_INTERVAL);
   }

   @Test
   public void testAdd() {
      final IntInterval i1 = new IntInterval(1, 3);
      final IntInterval i2 = new IntInterval(4, 8);
      final IntInterval i3 = new IntInterval(7, 9);
      final IntInterval i4 = new IntInterval(11, 13);
      final IntInterval i5 = new IntInterval(1, 13);
      Set<IntInterval> res = new HashSet<>();
      res = IntInterval.add(res, IntInterval.NULL_INTERVAL);
      assertEquals(0, res.size());
      res = IntInterval.add(res, i1);
      assertEquals(1, res.size());
      assertTrue(res.contains(i1));
      res = IntInterval.add(res, i2);
      assertEquals(2, res.size());
      assertTrue(res.contains(i1));
      assertTrue(res.contains(i2));
      res = IntInterval.add(res, i3);
      assertEquals(2, res.size());
      assertTrue(res.contains(i1));
      assertTrue(res.contains(new IntInterval(4, 9)));
      assertFalse(res.contains(i2));
      assertFalse(res.contains(i3));
      res = IntInterval.add(res, i4);
      assertEquals(3, res.size());
      assertTrue(res.contains(i1));
      assertTrue(res.contains(new IntInterval(4, 9)));
      assertTrue(res.contains(i4));
      res = IntInterval.add(res, i5);
      assertEquals(1, res.size());
      assertTrue(res.contains(i5));
   }

   @Test
   public void testContains() {
      final IntInterval i1 = new IntInterval(1, 3);
      assertTrue(i1.contains(1));
      assertTrue(i1.contains(2));
      assertTrue(i1.contains(3));
      assertFalse(i1.contains(0));
      assertFalse(i1.contains(-1));
      assertFalse(i1.contains(4));
      assertFalse(i1.contains(5));

      assertFalse(IntInterval.NULL_INTERVAL.contains(0));
      assertFalse(IntInterval.NULL_INTERVAL.contains(2));
      assertFalse(IntInterval.NULL_INTERVAL.contains(Integer.MIN_VALUE));
   }

   @Test
   public void testDistance() {
      final IntInterval i1 = new IntInterval(1, 3);
      assertEquals(0, i1.distance(1));
      assertEquals(0, i1.distance(2));
      assertEquals(0, i1.distance(3));
      assertEquals(1, i1.distance(4));
      assertEquals(5, i1.distance(8));
      assertEquals(1, i1.distance(0));
      assertEquals(5, i1.distance(-4));
      assertEquals(Integer.MAX_VALUE, IntInterval.NULL_INTERVAL.distance(0));
      assertEquals(Integer.MAX_VALUE, IntInterval.NULL_INTERVAL.distance(4));
      assertEquals(Integer.MAX_VALUE, IntInterval.NULL_INTERVAL.distance(-4));
   }
}