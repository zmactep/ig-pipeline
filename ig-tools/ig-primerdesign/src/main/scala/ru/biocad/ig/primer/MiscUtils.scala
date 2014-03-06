package ru.biocad.ig.primer

import scala.collection.mutable.{HashMap => mHashMap, Set => mSet, MultiMap => mMultiMap, Buffer}
import scala.util.Random

/**
 * @author kfeodorov
 */
object MiscUtils {
  def list2multimap[A, B](list: Buffer[(A, B)]): Map[A, Set[B]] =
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
}
