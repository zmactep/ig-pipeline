package igcont.kmer

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
  private val _buffer = ArrayBuffer.empty[Int]
  private val _sbuffer = ArrayBuffer.empty[Char]

  private val _didgitSize = math.ceil(math.log(_alpha.maxBy(_._2)._2)).toInt
  private val _indexMask = (1 to _ksize*_didgitSize).foldLeft(0)((acc,_) => acc << 1 | 1)

  def add(c : Char) : Unit = {
    _sbuffer += c
    if (_sbuffer.size >= _ksize) {
      if (_buffer.size == 0) {
        _buffer += _sbuffer.foldLeft(0)((prev, cc) =>
          prev << _didgitSize | _alpha.get(cc).get)
      } else {
        _buffer += (_buffer.last << _didgitSize | _alpha.get(c).get) & _indexMask
      }
    }
  }

  def get : Iterable[Int] = _buffer
}
