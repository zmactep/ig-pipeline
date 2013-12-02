package alicont

import scala.collection.immutable.HashMap
import scala.collection.mutable.ArrayBuffer

import alicont.common.Record

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 28.10.13
 * Time: 17:48
 */
class AlignmentResult(s : Double, q : String, t : String) {
  private val _query = q
  private val _target = t
  private val _result = new ArrayBuffer[(Char, HashMap[String, String])]()
  private val _score = s
  private val _similarity = 1.0 * q.zip(t).foldRight(0)((c, acc) => acc + (if (c._1 == c._2) 1 else 0)) / q.size
  private var _tname : String = null

  def this(s : Double, q : String, target_a : String, target : (String, Record)) = {
    this(s, q, target_a)
    _tname = target._2.name
    var i = target._1.indexOf(target_a.replaceAll("-", ""))

    assert(target._1.size == target._2.size)

    target_a.foreach(c => {
      if (c == '-') {
        _result += (('-', null))
      } else {
        _result += ((c, target._2.annotationOf(i)))
        i += 1
      }
    })
  }

  def get : Iterable[(Char, HashMap[String, String])] = _result

  def score : Double = _score

  def query : String = _query

  def target : String = _target

  def name : String = _tname

  def similarity : Double = _similarity
}
