package alicont.slow

import scala.Array
import scala.collection.mutable.ArrayBuffer

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 25.10.13
 * Time: 17:20
 */
class Matrix(width : Int) {
  private val _arrays = ArrayBuffer.empty[Array[Array[Int]]]
  private val _width = width

  def push(array : Array[Array[Int]]) : Unit = {
    _arrays += array
  }

  def push(height : Int) : Unit = {
    _arrays += Array.fill(height, _width)(0)
  }

  def pop() : Unit = {
    _arrays.remove(_arrays.size - 1)
  }

  def apply(h : Int) : Array[Int] = {
    val (a, i) = height2coord(h)
    _arrays(a)(i)
  }

  def apply(h : Int, w : Int) : Int = {
    val (a, i) = height2coord(h)
    _arrays(a)(i)(w)
  }

  def update(h : Int, w : Int, x : Int) : Unit = {
    val (a, i) = height2coord(h)
    _arrays(a)(i)(w) = x
  }

  def update(h : Int, x : Array[Int]) : Unit = {
    val (a, i) = height2coord(h)
    _arrays(a)(i) = x
  }

  def height2coord(h : Int) : (Int, Int) = {
    var tmp = h
    var arr = 0

    while (tmp >= _arrays(arr).size) {
      tmp -= _arrays(arr).size
      arr += 1
    }

    (arr, tmp)
  }

  def last : Array[Int] = _arrays.last.last

  def size : (Int, Int) = {
    (_arrays.foldRight(0)((arr, acc) => acc + arr.size), _width)
  }
}
