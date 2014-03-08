package ru.biocad.ig.primer

import scala.collection.mutable.{HashMap => mHashMap, Set => mSet, MultiMap => mMultiMap, Buffer}
import scala.util.Random
import scala.annotation.tailrec

/**
 * @author kfeodorov
 */
object MiscUtils {
  def list2multimap[A, B](list: Seq[(A, B)]): Map[A, Set[B]] =
    list.foldLeft(new mHashMap[A, mSet[B]] with mMultiMap[A, B]){(acc, pair) => acc.addBinding(pair._1, pair._2)}.map{case (k, v) => (k, v.toSet)}.toMap

  /**
   * Weighted random sampling
   * @param prob object -> weight
   * @tparam T object type
   * @return returns randomly chosen object according to weights
   */
  def sample[T](prob: Map[T, Int]): T = {
    val rand = Random.nextInt(prob.values.sum)
    var sum = 0
    for ((k, v) <- prob) {
      sum += v
      if (rand < sum) {
        return k
      }
    }
    prob.head._1
  }

  def countSubstring(str1: String, str2: String): Int = {
    @tailrec def count(pos: Int, c: Int): Int={
      val idx = str1.indexOf(str2, pos)
      if (idx == -1) c else count(idx + str2.size, c + 1)
    }
    count(0, 0)
  }

  def variance(data: List[Double]) = data.map{x => x * x}.sum / data.length - scala.math.pow(data.sum / data.length, 2)

}
