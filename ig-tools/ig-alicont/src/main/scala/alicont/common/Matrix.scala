package alicont.common

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 11.11.13
 * Time: 15:51
 */
class Matrix(width : Int, maxheight : Int) {
  private val _arrays    = Array.fill[Double](maxheight + 1, width + 1)(0)
  private val _width     = width
  private val _maxheight = maxheight
  private var _current   = -1

  def move(i : Int) : Unit = _current += i

  def apply(h : Int) : Array[Double] = _arrays(h)

  def apply(h : Int, w : Int) : Double = _arrays(h)(w)

  def update(h : Int, w : Int, x : Int) : Unit = _arrays(h)(w) = x

  def update(h : Int, x : Array[Double]) : Unit = _arrays(h) = x

  def pred : Array[Double] = _arrays(_current - 1)

  def last : Array[Double] = _arrays(_current)

  def height : Int = _current + 1

  def size : (Int, Int) = {
    (_current + 1, _width)
  }

  def fullsize : (Int, Int) = (_maxheight, _width)
}
