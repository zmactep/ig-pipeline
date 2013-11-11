package igcont.kmer.slow

import scala.collection.immutable.TreeMap
import scala.collection.mutable.ArrayBuffer

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 08.10.13
 * Time: 19:48
 */
class Hasher(k : Int, alpha : TreeMap[Char, Int]) {
  private val _ksize = k
  private val _alpha = alpha
  private val _asize = alpha.size
  private val _buffer = ArrayBuffer.empty[Int]
  private val _sbuffer = ArrayBuffer.empty[Char]

  def add(c : Char) : Unit = {
    _sbuffer += c
    if (_sbuffer.size >= _ksize) {
      if (_buffer.size == 0) {
        _buffer += _sbuffer.zipWithIndex.foldRight(0)((cc, prev) =>
          prev + _alpha.get(cc._1).get * math.pow(_asize, cc._2).toInt)
      } else {
        val elem = _alpha.get(_sbuffer(_sbuffer.size - _ksize - 1)).get
        val last = _buffer.last
        val prenew  = (last - elem) / _asize
        _buffer += (prenew + _alpha.get(c).get * math.pow(_asize, _ksize - 1).toInt)
      }
    }
  }

  def get : Iterable[Int] = _buffer
}
