package ru.biocad.ig.primer

import scala.collection.mutable.{HashMap => mHashMap, Set => mSet, MultiMap => mMultiMap, Buffer}

/**
 * Created by Kos on 22.02.14.
 */
object MiscUtils {
  def list2multimap[A, B](list: Buffer[(A, B)]): Map[A, Set[B]] =
    list.foldLeft(new mHashMap[A, mSet[B]] with mMultiMap[A, B]){(acc, pair) => acc.addBinding(pair._1, pair._2)}.map{case (k, v) => (k, v.toSet)}.toMap

}
