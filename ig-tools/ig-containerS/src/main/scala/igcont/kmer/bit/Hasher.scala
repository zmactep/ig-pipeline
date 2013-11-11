package igcont.kmer.bit

import scala.collection.immutable.TreeMap

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 11.11.13
 * Time: 9:37
 */
class Hasher(k : Int, alpha : TreeMap[Char, Int]) {
  private val _ksize   = k
  private val _alpha   = alpha
  private val _bsize   = math.ceil(math.log(alpha.size)/math.log(2)).toInt
  private val _bmask   = Integer.parseInt("1" * _bsize * _ksize, 2)
  private var _buffer  = 0
  private var _counter = 0

  def add(c : Char) : Int = {
    _buffer = ((_buffer << _bsize) + _alpha(c)) & _bmask
    _counter += 1
    if (_counter < _ksize) -1 else _buffer
  }

  def clear() = {
    _buffer = 0
    _counter = 0
  }
}
