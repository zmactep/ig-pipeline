package alicont.fast

import scala.collection.mutable

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 11.11.13
 * Time: 16:47
 */
class Alicont(maxheight : Int, query : String, gap : Int, score_matrix : Array[Array[Int]]) {
  private val _query   = query
  private val _gap     = gap
  private val _score   = score_matrix
  private val _matrix  = new Matrix(query.size, maxheight)
  private val _strings = new mutable.Stack[String]()

  def push(s : String) : Unit = {
    _strings.push(s)
    Algorithm.needleman_wunsch(s, _query, _gap, _score, _matrix)
  }

  def pop() : Unit = {
    val ls = _strings.pop()
    _matrix.move(-ls.size)
  }

  def target : String = _strings.reverse.mkString("")

  def alignment() : (Int, (String, String)) =
    (_matrix.last.last, Algorithm.traceback(target, _query, _gap, _score, _matrix))
}
