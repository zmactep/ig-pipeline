package igcont.anno

import scala.collection.mutable.ArrayBuffer

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 17.10.13
 * Time: 16:53
 */
class Annotation(n : Int, l : Int) {
  val _node = n
  val _vec  = ArrayBuffer.fill(l)(-1)

  def set(atype : Int, value : Int) = _vec(atype) = value

  def node() : Int = _node

  def data() : ArrayBuffer[Int] = _vec

  def dataOf(atype : Int) : Int = _vec(atype)
}