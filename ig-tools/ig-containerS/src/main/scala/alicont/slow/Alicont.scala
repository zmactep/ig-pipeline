package alicont.slow

import scala.collection.mutable

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 28.10.13
 * Time: 15:03
 */
class Alicont(query : String, gap : Int, score_matrix : Array[Array[Int]]) {
  private val _query   = query
  private val _gap     = gap
  private val _score   = score_matrix
  private val _matrix  = new Matrix(query.size)
  private val _strings = new mutable.Stack[String]()

  // Just for capability with fast version
  def this(maxsize : Int, query : String, gap : Int, score_matrix : Array[Array[Int]]) = {
    this(query, gap, score_matrix)
  }

  def push(s : String) : Unit = {
    _strings.push(s)
    if (_strings.size == 1) {
      _matrix.push(Algorithm.needleman_wunsch(s, _query, _gap, _score))
    }
    else {
      _matrix.push(Algorithm.needleman_wunsch(s, _query, _gap, _score, _matrix.last))
    }
  }

  def pop() : Unit = {
    _strings.pop()
    _matrix.pop()
  }

  def target : String = _strings.reverse.mkString("")

  def alignment() : (Int, (String, String)) =
    (_matrix.last.last, Algorithm.traceback(target, _query, _gap, _score, _matrix))
}
