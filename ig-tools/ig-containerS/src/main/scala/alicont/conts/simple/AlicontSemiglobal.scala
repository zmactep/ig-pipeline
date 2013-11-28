package alicont.conts.simple

import alicont.algorithms.simple.SemiglobalAlignment
import alicont.SimpleAlicont

/**
 * Created with IntelliJ IDEA.
 * User: pavel
 * Date: 27.11.13
 * Time: 23:37
 */
class AlicontSemiglobal(maxheight : Int, query : String, gap : Int, score_matrix : Array[Array[Int]])
  extends SimpleAlicont(maxheight : Int, query : String, gap : Int, score_matrix : Array[Array[Int]]) {

  def push(s : String) : Unit = {
    _strings.push(s)
    SemiglobalAlignment.extendMatrix(s, _query, _gap, _score, _scoreMatrix)
  }

  def alignment() : (Int, (String, String)) = {
    SemiglobalAlignment.traceback(target, _query, _gap, _score, _scoreMatrix)
  }
}