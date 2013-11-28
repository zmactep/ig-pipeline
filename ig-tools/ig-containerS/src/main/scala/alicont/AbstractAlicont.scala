package alicont

import scala.collection.mutable

/**
 * Created with IntelliJ IDEA.
 * User: pavel
 * Date: 27.11.13
 * Time: 22:53
 */
abstract class AbstractAlicont(maxheight : Int, query : String, score_matrix : Array[Array[Int]]) {
  protected val _query   = query
  protected val _score   = score_matrix
  protected val _strings = new mutable.Stack[String]()
  protected val _scoreMatrix  = new Matrix(query.size, maxheight)

  def push(s : String) : Unit
  def pop() : Unit
  def alignment() : (Int, (String, String))
  def target : String = _strings.reverse.mkString("")
}