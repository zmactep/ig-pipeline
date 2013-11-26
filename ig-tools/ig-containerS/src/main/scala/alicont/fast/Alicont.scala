package alicont.fast

import scala.collection.mutable

/**
 * Created with IntelliJ IDEA.
 * User: mactep
 * Date: 11.11.13
 * Time: 16:47
 */

abstract class AlicontBase(maxheight : Int, query : String, score_matrix : Array[Array[Int]]) {

  protected val _query   = query
  protected val _score   = score_matrix
  protected val _strings = new mutable.Stack[String]()
  protected val _scoreMatrix  = new Matrix(query.size, maxheight)

  def push(s : String) : Unit
  def pop() : Unit
  def alignment() : (Int, (String, String))
  def target : String = _strings.reverse.mkString("")
}

abstract class Alicont(maxheight : Int, query : String, gap : Int, score_matrix : Array[Array[Int]])
  extends AlicontBase(maxheight, query, score_matrix){

  protected val _gap     = gap

  def pop() : Unit = {
    val ls = _strings.pop()
    _scoreMatrix.move(-ls.size)
  }
}

abstract class AlicontAffineGap(maxheight : Int, query : String, gapOpen : Int, gapExtend : Int,
                                score_matrix : Array[Array[Int]])
  extends AlicontBase(maxheight, query, score_matrix) {

  protected val _gapOpen = gapOpen
  protected val _gapExtend = gapExtend
  protected val _deletionMatrix = new Matrix(query.size, maxheight)
  protected val _insertionMatrix = new Matrix(query.size, maxheight)

  def pop() : Unit = {
    val ls = _strings.pop()
    _scoreMatrix.move(-ls.size)
    _deletionMatrix.move(-ls.size)
    _insertionMatrix.move(-ls.size)
  }
}

class AlicontGlobal(maxheight : Int, query : String, gap : Int, score_matrix : Array[Array[Int]])
  extends Alicont(maxheight, query, gap, score_matrix) {

  def push(s : String) : Unit = {
    _strings.push(s)
    Algorithm.NeedlemanWunsch.extendAlign(s, _query, _gap, _score, _scoreMatrix)
  }

  def alignment() : (Int, (String, String)) =
    Algorithm.NeedlemanWunsch.traceback(target, _query, _gap, _score, _scoreMatrix)
}

class AlicontLocal(maxheight : Int, query : String, gap : Int, score_matrix : Array[Array[Int]])
  extends Alicont(maxheight, query, gap, score_matrix) {

  def push(s : String) : Unit = {
    _strings.push(s)
    Algorithm.SmithWaterman.extendAlign(s, _query, _gap, _score, _scoreMatrix)
  }

  def alignment() : (Int, (String, String)) =
    Algorithm.SmithWaterman.traceback(target, _query, _gap, _score, _scoreMatrix)
}
