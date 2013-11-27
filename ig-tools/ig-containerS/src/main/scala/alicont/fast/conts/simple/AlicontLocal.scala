package alicont.fast.conts.simple

import alicont.fast.SimpleAlicont
import alicont.fast.algorithms.simple.LocalAlignment

/**
 * Created with IntelliJ IDEA.
 * User: pavel
 * Date: 27.11.13
 * Time: 23:36
 */
class AlicontLocal(maxheight : Int, query : String, gap : Int, score_matrix : Array[Array[Int]])
  extends SimpleAlicont(maxheight : Int, query : String, gap : Int, score_matrix : Array[Array[Int]]) {

  def push(s : String) : Unit = {
    _strings.push(s)
    LocalAlignment.extendMatrix(s, _query, _gap, _score, _scoreMatrix)
  }

  def alignment() : (Int, (String, String)) = {
    LocalAlignment.traceback(target, _query, _gap, _score, _scoreMatrix)
  }
}
